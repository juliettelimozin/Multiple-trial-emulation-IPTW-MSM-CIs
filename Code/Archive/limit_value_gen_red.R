library(modelr)
library(tidyverse)
library(tidyr)
source("simulate_MSM_simplified.R")
setwd('~/rds/hpc-work')
set.seed(20222022)
library(MASS)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(pammtools)

conf <- 1:9/10
treat <- c(-1,-0.8,-0.5,-0.2,0,0.2,0.5,0.8,1)

limit_value_conf_red <- array(,dim = c(5,2,9))
limit_value_treat_red <- array(, dim = c(5,2,9))

for (i in conf){
  ##################################################################
  simdata_censored_conf<-DATA_GEN_censored_reduced(1000000, 5, conf = i)
  
  PP_prep <- RandomisedTrialsEmulation::data_preparation(simdata_censored_conf, id='ID', period='t', treatment='A', outcome='Y', 
                                                         eligible ='eligible',cense = 'C',
                                                         switch_d_cov = ~X2 + X4,
                                                         outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                                         cense_d_cov = ~X2 + X4,
                                                         include_expansion_time_case = ~1, include_followup_time_case = ~1, 
                                                         include_regime_length = F,
                                                         use_weight=1, use_censor=1, numCores = 1, quiet = T)
  switch_data <- PP_prep$data %>% 
    dplyr::mutate(t_1 = ifelse(followup_time == 1,1,0),
                  t_2 = ifelse(followup_time == 2,1,0),
                  t_3 = ifelse(followup_time == 3,1,0),
                  t_4 = ifelse(followup_time == 4,1,0),
                  t_1A = t_1*assigned_treatment,
                  t_2A = t_2*assigned_treatment,
                  t_3A = t_3*assigned_treatment,
                  t_4A = t_4*assigned_treatment,
                  t_1X2 = t_1*X2,
                  t_2X2 = t_2*X2,
                  t_3X2 = t_3*X2,
                  t_4X2 = t_4*X2,
                  t_1X4 = t_1*X4,
                  t_2X4 = t_2*X4,
                  t_3X4 = t_3*X4,
                  t_4X4 = t_4*X4)
  
  PP <- RandomisedTrialsEmulation::data_modelling(data = switch_data,
                                                  outcome_cov = ~ X2 + X4+ assigned_treatment+
                                                    t_1 + t_2 + t_3 + t_4 +
                                                    t_1A + t_2A + t_3A + t_4A + 
                                                    t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
                                                    t_1X4 + t_2X4 + t_3X4 + t_4X4,
                                                  model_var = c('assigned_treatment'),
                                                  glm_function = 'glm',
                                                  include_expansion_time_case = ~1, include_followup_time_case = ~1,
                                                  use_weight=1, use_censor=1, numCores = 1, quiet = T, use_sample_weights =  F)
  
  #### Survival function point estimate for PP ####
  design_mat <- expand.grid(id = 1:1000000,
                            for_period = 0:4,
                            followup_time = 0:4) 
  design_mat <- design_mat[which(5 -design_mat$for_period > design_mat$followup_time),]
  
  fitting_data_treatment <-  switch_data %>% 
    dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
    dplyr::select(id,for_period, followup_time, X2,  X4, assigned_treatment) %>% 
    merge(design_mat, by = c("id", "for_period", "followup_time"), all.y = TRUE) %>% 
    dplyr::group_by(id) %>% 
    tidyr::fill( X2,X4,assigned_treatment,.direction = "down") %>% 
    dplyr::ungroup() %>% 
    dplyr::select(id, for_period, followup_time, X2, X4, assigned_treatment) %>% 
    merge(data.frame(id = switch_data$id, for_period = switch_data$for_period), by = c("id", "for_period"), all.y = TRUE) %>% 
    dplyr::arrange(id, for_period, followup_time) %>% 
    dplyr::mutate(t_1 = ifelse(followup_time == 1,1,0),
                  t_2 = ifelse(followup_time == 2,1,0),
                  t_3 = ifelse(followup_time == 3,1,0),
                  t_4 = ifelse(followup_time == 4,1,0),
                  t_1A = t_1*assigned_treatment,
                  t_2A = t_2*assigned_treatment,
                  t_3A = t_3*assigned_treatment,
                  t_4A = t_4*assigned_treatment,
                  t_1X2 = t_1*X2,
                  t_2X2 = t_2*X2,
                  t_3X2 = t_3*X2,
                  t_4X2 = t_4*X2,
                  t_1X4 = t_1*X4,
                  t_2X4 = t_2*X4,
                  t_3X4 = t_3*X4,
                  t_4X4 = t_4*X4) %>% 
    dplyr::filter(for_period == 0)
  
  
  fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]
  
  fitting_data_control <- fitting_data_treatment %>% 
    dplyr::mutate(assigned_treatment = assigned_treatment*0,
                  t_1A = t_1*0,
                  t_2A = t_2*0,
                  t_3A = t_3*0,
                  t_4A = t_4*0,)
  
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
  
  limit_value_conf_red[,1,10*i] <- pull(predicted_probas_PP_sample, survival_treatment)
  limit_value_conf_red[,2,10*i] <- pull(predicted_probas_PP_sample, survival_control)
  
}
save(limit_value_conf_red, file = "limit_value_conf_red.rda")

for (i in 1:9){
  ############################SANDWICH CI ######################################
  simdata_censored_treat<-DATA_GEN_censored_reduced(1000000, 5, treat_prev = treat[i])
  
  PP_prep <- RandomisedTrialsEmulation::data_preparation(simdata_censored_treat, id='ID', period='t', treatment='A', outcome='Y', 
                                                         eligible ='eligible',cense = 'C',
                                                         switch_d_cov = ~X2 + X4,
                                                         outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                                         cense_d_cov = ~X2 + X4,
                                                         include_expansion_time_case = ~1, include_followup_time_case = ~1, 
                                                         include_regime_length = F,
                                                         use_weight=1, use_censor=1, numCores = 1, quiet = T)
  switch_data <- PP_prep$data %>% 
    dplyr::mutate(t_1 = ifelse(followup_time == 1,1,0),
                  t_2 = ifelse(followup_time == 2,1,0),
                  t_3 = ifelse(followup_time == 3,1,0),
                  t_4 = ifelse(followup_time == 4,1,0),
                  t_1A = t_1*assigned_treatment,
                  t_2A = t_2*assigned_treatment,
                  t_3A = t_3*assigned_treatment,
                  t_4A = t_4*assigned_treatment,
                  t_1X2 = t_1*X2,
                  t_2X2 = t_2*X2,
                  t_3X2 = t_3*X2,
                  t_4X2 = t_4*X2,
                  t_1X4 = t_1*X4,
                  t_2X4 = t_2*X4,
                  t_3X4 = t_3*X4,
                  t_4X4 = t_4*X4)
  
  PP <- RandomisedTrialsEmulation::data_modelling(data = switch_data,
                                                  outcome_cov = ~ X2 + X4+ assigned_treatment+
                                                    t_1 + t_2 + t_3 + t_4 +
                                                    t_1A + t_2A + t_3A + t_4A + 
                                                    t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
                                                    t_1X4 + t_2X4 + t_3X4 + t_4X4,
                                                  model_var = c('assigned_treatment'),
                                                  glm_function = 'glm',
                                                  include_expansion_time_case = ~1, include_followup_time_case = ~1,
                                                  use_weight=1, use_censor=1, numCores = 1, quiet = T, use_sample_weights =  F)
  
  #### Survival function point estimate for PP ####
  design_mat <- expand.grid(id = 1:1000000,
                            for_period = 0:4,
                            followup_time = 0:4) 
  design_mat <- design_mat[which(5 -design_mat$for_period > design_mat$followup_time),]
  
  fitting_data_treatment <-  switch_data %>% 
    dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
    dplyr::select(id,for_period, followup_time, X2,  X4, assigned_treatment) %>% 
    merge(design_mat, by = c("id", "for_period", "followup_time"), all.y = TRUE) %>% 
    dplyr::group_by(id) %>% 
    tidyr::fill( X2,X4,assigned_treatment,.direction = "down") %>% 
    dplyr::ungroup() %>% 
    dplyr::select(id, for_period, followup_time, X2, X4, assigned_treatment) %>% 
    merge(data.frame(id = switch_data$id, for_period = switch_data$for_period), by = c("id", "for_period"), all.y = TRUE) %>% 
    dplyr::arrange(id, for_period, followup_time) %>% 
    dplyr::mutate(t_1 = ifelse(followup_time == 1,1,0),
                  t_2 = ifelse(followup_time == 2,1,0),
                  t_3 = ifelse(followup_time == 3,1,0),
                  t_4 = ifelse(followup_time == 4,1,0),
                  t_1A = t_1*assigned_treatment,
                  t_2A = t_2*assigned_treatment,
                  t_3A = t_3*assigned_treatment,
                  t_4A = t_4*assigned_treatment,
                  t_1X2 = t_1*X2,
                  t_2X2 = t_2*X2,
                  t_3X2 = t_3*X2,
                  t_4X2 = t_4*X2,
                  t_1X4 = t_1*X4,
                  t_2X4 = t_2*X4,
                  t_3X4 = t_3*X4,
                  t_4X4 = t_4*X4) %>% 
    dplyr::filter(for_period == 0)
  
  
  fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]
  
  fitting_data_control <- fitting_data_treatment %>% 
    dplyr::mutate(assigned_treatment = assigned_treatment*0,
                  t_1A = t_1*0,
                  t_2A = t_2*0,
                  t_3A = t_3*0,
                  t_4A = t_4*0,)
  
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
  
  limit_value_treat_red[,1,i] <- pull(predicted_probas_PP_sample, survival_treatment)
  limit_value_treat_red[,2,i] <- pull(predicted_probas_PP_sample, survival_control)
}

save(limit_value_treat_red, file = 'limit_value_treat_red.rda')