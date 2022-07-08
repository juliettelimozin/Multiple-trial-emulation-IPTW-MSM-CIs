#!/usr/bin R
library(modelr)
library(tidyverse)
library(tidyr)
setwd("~/rds/hpc-work")
source("simulate_MSM.R")
set.seed(NULL)
library(RandomisedTrialsEmulation)
library(MASS)
library(sandwich)
library(doParallel)
library(doRNG)

conf <- 1:9/10
est_true_value_conf <- array(,dim = c(10,2,9))
est_true_value_treat <- array(,dim = c(10,2,9))
nsample <- 1000000
for (i in conf){
  simdata_censored_conf<-DATA_GEN_censored(nsample, 10, conf = i, censor = F)
  ##### excluding obs before first becoming eligible, with censoring by dropout #######################################
  
  #### PP analysis
  
  PP <-RandomisedTrialsEmulation::initiators( simdata_censored_conf,id='ID', period='t', treatment='A', outcome='Y', eligible ='eligible',
                                              model_switchd =c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                              cov_switchd = c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                              outcomeCov_var=c('X1', 'X2','X3', 'X4', 'age_s'), outcomeCov =c('X1', 'X2','X3', 'X4', 'age_s'), model_var = c('assigned_treatment'),
                                              include_expansion_time_case = 0, include_followup_time_case = c("linear", "quadratic"), include_regime_length = 1,
                                              use_weight=1, use_censor=1, case_control = 0, data_dir =data_direction, numCores = 1, quiet = TRUE)
  
  design_mat <- expand.grid(id = 1:nsample,
                            for_period = 0:9,
                            followup_time = 0:9) %>% 
    dplyr::mutate(followup_time2 = followup_time^2)
  design_mat <- design_mat[which(10 -design_mat$for_period > design_mat$followup_time),]
  
  switch_data <- PP$model$switch_data
  
  fitting_data_treatment <-  switch_data %>% 
    dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
    dplyr::select(id,for_period, followup_time, followup_time2, X1, X2,X3, X4, age_s, assigned_treatment) %>% 
    merge(design_mat, by = c("id", "for_period", "followup_time", "followup_time2"), all.y = TRUE) %>% 
    dplyr::group_by(id) %>% 
    tidyr::fill(X1, X2,X3,X4,age_s,assigned_treatment,.direction = "down") %>% 
    dplyr::ungroup() %>% 
    dplyr::select(id, for_period, followup_time, followup_time2, X1, X2,X3, X4, age_s, assigned_treatment) %>% 
    merge(data.frame(id = switch_data$id, for_period = switch_data$for_period), by = c("id", "for_period"), all.y = TRUE) %>% 
    dplyr::arrange(id, for_period, followup_time) %>% 
    dplyr::filter(for_period == 0)
  
  fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]
  fitting_data_treatment <- fitting_data_treatment[which(!is.na(fitting_data_treatment$X3)),]
  
  fitting_data_control <- fitting_data_treatment %>% 
    dplyr::mutate(assigned_treatment = assigned_treatment*0)
  
  Y_pred_PP_treatment <- predict.glm(PP$model$model, fitting_data_treatment, 
                                     type = "response")
  Y_pred_PP_control <- predict.glm(PP$model$model, fitting_data_control, 
                                   type = "response")
  predicted_probas_PP <- fitting_data_treatment %>% 
    dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment,
                  predicted_proba_control = Y_pred_PP_control) %>% 
    dplyr::group_by(id, for_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  
  est_true_value_conf[,1,10*i] <- predicted_probas_PP$survival_treatment
  est_true_value_conf[,2,10*i] <- predicted_probas_PP$survival_control
  
  ##################################################################################################################
  simdata_censored_treat<-DATA_GEN_censored(nsample, 10, treat_prev = i, censor = F)
  
  
  PP <-RandomisedTrialsEmulation::initiators( simdata_censored_treat,id='ID', period='t', treatment='A', outcome='Y', eligible ='eligible',
                                              model_switchd =c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                              cov_switchd = c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                              outcomeCov_var=c('X1', 'X2','X3', 'X4', 'age_s'), outcomeCov =c('X1', 'X2','X3', 'X4', 'age_s'), model_var = c('assigned_treatment'),
                                              include_expansion_time_case = 0, include_followup_time_case = c("linear", "quadratic"), include_regime_length = 1,
                                              use_weight=1, use_censor=1, case_control = 0, data_dir =data_direction, numCores = 1, quiet = TRUE)
 
  design_mat <- expand.grid(id = 1:nsample,
                            for_period = 0:9,
                            followup_time = 0:9) %>% 
    dplyr::mutate(followup_time2 = followup_time^2)
  design_mat <- design_mat[which(10 -design_mat$for_period > design_mat$followup_time),]
  
  switch_data <- PP$model$switch_data
  
  fitting_data_treatment <-  switch_data %>% 
    dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
    dplyr::select(id,for_period, followup_time, followup_time2, X1, X2,X3, X4, age_s, assigned_treatment) %>% 
    merge(design_mat, by = c("id", "for_period", "followup_time", "followup_time2"), all.y = TRUE) %>% 
    dplyr::group_by(id) %>% 
    tidyr::fill(X1, X2,X3,X4,age_s,assigned_treatment,.direction = "down") %>% 
    dplyr::ungroup() %>% 
    dplyr::select(id, for_period, followup_time, followup_time2, X1, X2,X3, X4, age_s, assigned_treatment) %>% 
    merge(data.frame(id = switch_data$id, for_period = switch_data$for_period), by = c("id", "for_period"), all.y = TRUE) %>% 
    dplyr::arrange(id, for_period, followup_time) %>% 
    dplyr::filter(for_period == 0)
  
  fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]
  fitting_data_treatment <- fitting_data_treatment[which(!is.na(fitting_data_treatment$X3)),]
  
  fitting_data_control <- fitting_data_treatment %>% 
    dplyr::mutate(assigned_treatment = assigned_treatment*0)
  
  Y_pred_PP_treatment <- predict.glm(PP$model$model, fitting_data_treatment, 
                                     type = "response")
  Y_pred_PP_control <- predict.glm(PP$model$model, fitting_data_control, 
                                   type = "response")
  predicted_probas_PP <- fitting_data_treatment %>% 
    dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment,
                  predicted_proba_control = Y_pred_PP_control) %>% 
    dplyr::group_by(id, for_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  
  est_true_value_treat[,1,10*i] <- predicted_probas_PP$survival_treatment
  est_true_value_treat[,2,10*i] <- predicted_probas_PP$survival_control
}

save(est_true_value_conf, file = "est_true_value_conf_noint.rda")
save(est_true_value_treat, file = "est_true_value_treat_noint.rda")
