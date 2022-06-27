#!/usr/bin R
library(modelr)
library(tidyverse)
library(tidyr)
setwd("~/rds/hpc-work")
source("simulate_MSM.R")
set.seed(20222022)
library(RandomisedTrialsEmulation)
library(MASS)
library(sandwich)
library(foreach)
library(doParallel)
library(parallel)

data.table::setDTthreads(20)

bootstrap_construct <- function(data){
  sampleindex<-sample(unique(data$ID), length(unique(data$ID)),replace=T)
  sampleindex<-sampleindex[order(sampleindex)]
  
  dataindex=newsub=NULL
  
  for (i in 1:length(unique(data$ID)))
  {
    index<-which(data$ID==sampleindex[i])
    dataindex<-c(dataindex,index)
    newsub<-c(newsub, rep(i,length(index)))
  }
  
  bootdata<-data[dataindex,]
  bootdata$ID<-newsub
  return(bootdata)
}

iters <- 1

CI_bootstrap_coefs_PP <- array(,dim = c(10,2,iters))
CI_sandwich_coefs_PP <- array(, dim = c(10,2,iters))

l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
j <- as.numeric(l/10)

cpu_num <- detectCores()
data_direction <- paste0("~/rds/hpc-work/data_",l,"/conf")

for (i in 1:iters){
  ##################### CONFOUNDING STRENGTH #######################################
  
  simdata_censored_conf<-DATA_GEN_censored(1000, 10, conf = j)
  
  ################################### Bootstrap CI ###################################
  bootstrap_iter <- 100
  
  surv_PP_difference_boostrap_estimates_conf <-as.data.frame(matrix(,10,bootstrap_iter))
  
  for (k in 1:bootstrap_iter){
    boot_data_conf <- bootstrap_construct(simdata_censored_conf)
    write.csv(boot_data_conf, paste0(data_direction,"/boot_data_conf.csv"),row.names = F)
    
    #########################################
    
    PP_boot_conf <- RandomisedTrialsEmulation::initiators(paste0(data_direction,"/boot_data_conf.csv"),id='ID', period='t', treatment='A', outcome='Y', eligible ='eligible', cense = 'C',
                                                          model_switchd =c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                                          cov_switchd = c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                                          outcomeCov_var=c( 'X3', 'X4', 'age_s'), outcomeCov =c('X3', 'X4', 'age_s'), model_var = c('assigned_treatment'),
                                                          cov_censed = c( 'X1', 'X2','X3', 'X4', 'age_s'), model_censed =c( 'X1', 'X2','X3', 'X4', 'age_s'), pool_cense=1,
                                                          include_expansion_time_case = 0, include_followup_time_case = c("linear", "quadratic"), include_regime_length = 1,
                                                          use_weight=1, use_censor=1, case_control = 0, data_dir = data_direction, numCores = cpu_num, quiet = TRUE)
    
    switch_data_boot <- read.csv(paste0(data_direction,"/switch_data.csv"))
    design_mat <- expand.grid(id = 1:tail(switch_data_boot$id)[1],
                              for_period = 0:9,
                              followup_time = 0:9) %>% 
      dplyr::mutate(followup_time2 = followup_time^2)
    design_mat <- design_mat[which(10 -design_mat$for_period > design_mat$followup_time),]
    
    fitting_data_treatment_boot <- switch_data_boot %>% 
      dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
      dplyr::select(id,for_period, followup_time, followup_time2, X3, X4, age_s, assigned_treatment) %>% 
      merge(design_mat, by = c("id", "for_period", "followup_time", "followup_time2"), all.y = TRUE) %>% 
      dplyr::group_by(id) %>% 
      tidyr::fill(X3,X4,age_s,assigned_treatment,.direction = "down") %>% 
      dplyr::ungroup() %>% 
      dplyr::select(id, for_period, followup_time, followup_time2, X3, X4, age_s, assigned_treatment) %>% 
      merge(data.frame(id = switch_data_boot$id, for_period = switch_data_boot$for_period), by = c("id", "for_period"), all.y = TRUE) %>% 
      dplyr::arrange(id, for_period, followup_time)
    
    
    fitting_data_treatment_boot <- fitting_data_treatment_boot[!duplicated(fitting_data_treatment_boot),]
    fitting_data_treatment_boot <- fitting_data_treatment_boot[which(!is.na(fitting_data_treatment_boot$X3)),]
    
    fitting_data_control_boot <- fitting_data_treatment_boot %>% 
      dplyr::mutate(assigned_treatment = assigned_treatment*0)
    
    Y_pred_PP_treatment_boot <- predict.glm(PP_boot_conf$model$model, 
                                            fitting_data_treatment_boot, 
                                            type = "response")
    Y_pred_PP_control_boot <- predict.glm(PP_boot_conf$model$model, 
                                          fitting_data_control_boot,
                                          type = "response")
    predicted_probas_PP_boot <- fitting_data_treatment_boot %>% 
      dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment_boot,
                    predicted_proba_control = Y_pred_PP_control_boot) %>% 
      dplyr::group_by(id, for_period) %>% 
      dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                    cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(followup_time) %>% 
      dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                       survival_control = mean(cum_hazard_control))
    
    surv_PP_difference_boostrap_estimates_conf[,k] <- predicted_probas_PP_boot[,2] - predicted_probas_PP_boot[,3]
  }
  
  surv_PP_difference_boostrap_estimates_conf$lb <- apply(surv_PP_difference_boostrap_estimates_conf,
                                                         1,
                                                         quantile,
                                                         probs = c(0.025))
  surv_PP_difference_boostrap_estimates_conf$ub <- apply(surv_PP_difference_boostrap_estimates_conf,
                                                         1,
                                                         quantile,
                                                         probs = c(0.975))
  
  
  
  CI_bootstrap_coefs_PP[,1,i] <- surv_PP_difference_boostrap_estimates_conf$lb
  CI_bootstrap_coefs_PP[,2,i] <- surv_PP_difference_boostrap_estimates_conf$ub
  
  ############################SANDWICH #######################################
  
  write.csv(simdata_censored_conf, file= paste(data_direction,"/MSM_censor_conf.csv",sep = ""), row.names = F) 

  PP <-RandomisedTrialsEmulation::initiators(paste0(data_direction,"/MSM_censor_conf.csv"), id='ID', period='t', treatment='A', outcome='Y', eligible ='eligible', cense = 'C',
                                             model_switchd =c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                             cov_switchd = c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                             outcomeCov_var=c( 'X3', 'X4', 'age_s'), outcomeCov =c('X3', 'X4', 'age_s'), model_var = c('assigned_treatment'),
                                             cov_censed = c( 'X1', 'X2','X3', 'X4', 'age_s'), model_censed =c( 'X1', 'X2','X3', 'X4', 'age_s'), pool_cense=1,
                                             include_expansion_time_case = 0, include_followup_time_case = c("linear", "quadratic"), include_regime_length = 1,
                                             use_weight=1, use_censor=1, case_control = 0, data_dir =data_direction, numCores = cpu_num, quiet = TRUE)
  
  design_mat <- expand.grid(id = 1:1000,
                            for_period = 0:9,
                            followup_time = 0:9) %>% 
    dplyr::mutate(followup_time2 = followup_time^2)
  design_mat <- design_mat[which(10 -design_mat$for_period > design_mat$followup_time),]
  
  switch_data <- read.csv(paste0(data_direction,"/switch_data.csv"))
  
  fitting_data_treatment <-  switch_data %>% 
    dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
    dplyr::select(id,for_period, followup_time, followup_time2, X3, X4, age_s, assigned_treatment) %>% 
    merge(design_mat, by = c("id", "for_period", "followup_time", "followup_time2"), all.y = TRUE) %>% 
    dplyr::group_by(id) %>% 
    tidyr::fill(X3,X4,age_s,assigned_treatment,.direction = "down") %>% 
    dplyr::ungroup() %>% 
    dplyr::select(id, for_period, followup_time, followup_time2, X3, X4, age_s, assigned_treatment) %>% 
    merge(data.frame(id = switch_data$id, for_period = switch_data$for_period), by = c("id", "for_period"), all.y = TRUE) %>% 
    dplyr::arrange(id, for_period, followup_time)
  
  fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]
  fitting_data_treatment <- fitting_data_treatment[which(!is.na(fitting_data_treatment$X3)),]
  
  fitting_data_control <- fitting_data_treatment %>% 
    dplyr::mutate(assigned_treatment = assigned_treatment*0)
  
  covariance_mat_PP <- PP$model$robust$matrix
  
  #Step 1 of algorithm  -- sampling Y_n1, ..., Y_nB ~ MN(coeffs,sandwich covariance)
  sampling_size <- 200
  
  coeffs_sample_PP <- mvrnorm(sampling_size,PP$model$model$coefficients, covariance_mat_PP)
  
  surv_PP_difference_sandwich_estimates <- as.data.frame(matrix(,10,sampling_size))
  
  for (k in 1:sampling_size){
    #Step 1 of algorithm -- same model with new coeffs = one point from MVN sample
    fit_sample <- PP
    fit_sample$model$model$coefficients <- coeffs_sample_PP[k,]
    
    #Step 2 -- calculating survival probas with new model
    Y_pred_sample_treatment <- predict.glm(fit_sample$model$model, 
                                           fitting_data_treatment, 
                                           type = "response")
    Y_pred_sample_control <- predict.glm(fit_sample$model$model, 
                                         fitting_data_control,
                                         type = "response")
    
    predicted_probas_PP_sample <- fitting_data_treatment %>% 
      dplyr::mutate(predicted_proba_treatment = Y_pred_sample_treatment,
                    predicted_proba_control = Y_pred_sample_control) %>% 
      dplyr::group_by(id, for_period) %>% 
      dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                    cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(followup_time) %>% 
      dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                       survival_control = mean(cum_hazard_control))
    surv_PP_difference_sandwich_estimates[,k] <- predicted_probas_PP_sample[,2] - predicted_probas_PP_sample[,3]
    
  }
  
  #Step 3 -- calculating lower and upper bounds by 2.5% and 97.5% quantiles
  
  surv_PP_difference_sandwich_estimates$lb <- apply(surv_PP_difference_sandwich_estimates,
                                                    1,
                                                    quantile,
                                                    probs = c(0.025))
  surv_PP_difference_sandwich_estimates$ub <- apply(surv_PP_difference_sandwich_estimates,
                                                    1,
                                                    quantile,
                                                    probs = c(0.975))
  
  
  
  CI_sandwich_coefs_PP[,1,i] <- surv_PP_difference_sandwich_estimates$lb
  CI_sandwich_coefs_PP[,2,i] <- surv_PP_difference_sandwich_estimates$ub
  
}

save(CI_bootstrap_coefs_PP, file = paste("CI_bootstrap_coefs_PP_",as.character(l),".rda", sep = ""))
save(CI_sandwich_coefs_PP, file = paste("CI_sandwich_coefs_PP_",as.character(l),".rda", sep = ""))


