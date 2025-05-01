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

#Number of MC iterations
iters <- 1000
size_PP <- matrix(,10,iters)

#Fetching array value from HPC parallelisation
l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
j <- 100*l

data_direction <- paste("~/rds/hpc-work/data_",l,sep = "")

# Set number of cores. 67 is sufficient for 200 cores.
registerDoParallel(cores = 67)

#Loop of MC simulations
for (i in 1:iters){
  ##################### CONFOUNDING STRENGTH #######################################
  #Catch error messages from non convergence or non positive definite matrix
  tryCatch({
    print(i)
    #Generate simulated data with specific confounding strength
    simdata_censored_size<-DATA_GEN_censored(j, 10)
    
    ############################SANDWICH CI ######################################
    
    PP_prep <- RandomisedTrialsEmulation::data_preparation(simdata_censored_size, id='ID', period='t', treatment='A', outcome='Y', eligible ='eligible', cense = 'C',
                                                           model_switchd =c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                                           cov_switchd = c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                                           outcomeCov_var=c('X1', 'X2', 'X3', 'X4', 'age_s'), outcomeCov =c('X1', 'X2','X3', 'X4', 'age_s'), model_var = c('assigned_treatment'),
                                                           cov_censed = c( 'X1', 'X2','X3', 'X4', 'age_s'), model_censed =c( 'X1', 'X2','X3', 'X4', 'age_s'), pool_cense=1,
                                                           include_expansion_time_case = 0, include_followup_time_case = c("linear", "quadratic"), include_regime_length = 1,
                                                           use_weight=1, use_censor=1, data_dir = data_direction, numCores = 1, quiet = TRUE)
    switch_data <- PP_prep$switch_data %>% 
      dplyr::mutate(tA = followup_time*assigned_treatment, 
                    tX1 = followup_time*X1,
                    tX2 = followup_time*X2,
                    tX3 = followup_time*X3,
                    tX4 = followup_time*X4,
                    tage_s = followup_time*age_s) 
    
    PP <- RandomisedTrialsEmulation::data_modelling(switch_data = switch_data,
                                                    outcomeCov_var=c('X1', 'X2', 'X3', 'X4', 'age_s', 'tA', 'tX1', 'tX2', 'tX3', 'tX4', 'tage_s'),
                                                    outcomeCov =c('X1', 'X2', 'X3', 'X4', 'age_s', 'tA', 'tX1', 'tX2', 'tX3', 'tX4', 'tage_s'), model_var = c('assigned_treatment'),
                                                    include_expansion_time_case = 0, include_followup_time_case = c("linear", "quadratic"),
                                                    use_weight=1, use_censor=1, numCores = 1, quiet = TRUE, use_sample_weights =  F)
    
    #### Survival function point estimate for PP ####
    design_mat <- expand.grid(id = 1:j,
                              for_period = 0:9,
                              followup_time = 0:9) %>% 
      dplyr::mutate(followup_time2 = followup_time^2)
    design_mat <- design_mat[which(10 -design_mat$for_period > design_mat$followup_time),]
    
    fitting_data_treatment <-  switch_data %>% 
      dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
      dplyr::select(id,for_period, followup_time, followup_time2, X1, X2, X3, X4, age_s, assigned_treatment) %>% 
      merge(design_mat, by = c("id", "for_period", "followup_time", "followup_time2"), all.y = TRUE) %>% 
      dplyr::group_by(id) %>% 
      tidyr::fill( X1, X2,X3,X4,age_s,assigned_treatment,.direction = "down") %>% 
      dplyr::ungroup() %>% 
      dplyr::select(id, for_period, followup_time, followup_time2, X1, X2,X3, X4, age_s, assigned_treatment) %>% 
      merge(data.frame(id = switch_data$id, for_period = switch_data$for_period), by = c("id", "for_period"), all.y = TRUE) %>% 
      dplyr::arrange(id, for_period, followup_time) %>% 
      dplyr::mutate(tA = followup_time*assigned_treatment, 
                    tX1 = followup_time*X1,
                    tX2 = followup_time*X2,
                    tX3 = followup_time*X3,
                    tX4 = followup_time*X4,
                    tage_s = followup_time*age_s) %>% 
      dplyr::filter(for_period == 0)
    
    fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]
    fitting_data_treatment <- fitting_data_treatment[which(!is.na(fitting_data_treatment$X3)),]
    
    fitting_data_control <- fitting_data_treatment %>% 
      dplyr::mutate(assigned_treatment = assigned_treatment*0,
                    tA = followup_time*assigned_treatment)
    
    #Step 2 -- calculating survival probas with new model
    Y_pred_sample_treatment <- predict.glm(PP$model, 
                                           fitting_data_treatment, 
                                           type = "response")
    Y_pred_sample_control <- predict.glm(PP$model, 
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
                       survival_control = mean(cum_hazard_control),
                       survival_difference = survival_treatment - survival_control)
    
    size_PP[,i] <- pull(predicted_probas_PP_sample, survival_difference)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
save(size_PP, file = paste("size_PP_",as.character(l),".rda", sep = ""))



