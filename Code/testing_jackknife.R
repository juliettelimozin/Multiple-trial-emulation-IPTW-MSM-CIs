#!/usr/bin R
library(modelr)
library(tidyverse)
library(tidyr)
source("simulate_MSM_simplified.R")
source("weight_func.R")
set.seed(NULL)
library(TrialEmulation)
library(MASS)
library(sandwich)
library(doParallel)
library(doRNG)

iters <- 1000
bootstrap_iter <- 500

size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)

scenarios <- tidyr::crossing(size,conf, treat)

   

simdata_censored<-DATA_GEN_censored_reduced(200, 5, 
                                            treat_prev = as.numeric(0),
                                            outcome_prev = -3.8,
                                            censor = F)
PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                            eligible ='eligible',
                                            switch_d_cov = ~X2 + X4,
                                            outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                            use_weight=T, use_censor=T, quiet = T,
                                            save_weight_models = T,
                                            data_dir = data_direction)
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

PP <- TrialEmulation::trial_msm(data = switch_data,
                                outcome_cov = ~ X2 + X4+ assigned_treatment+
                                  t_1 + t_2 + t_3 + t_4 +
                                  t_1A + t_2A + t_3A + t_4A + 
                                  t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
                                  t_1X4 + t_2X4 + t_3X4 + t_4X4,
                                model_var = c('assigned_treatment'),
                                glm_function = 'glm',
                                include_trial_period = ~1, include_followup_time = ~1,
                                use_weight=T, use_censor = T, quiet = T, use_sample_weights =  F)
switch_data$p_i <- predict.glm(PP$model, switch_data,type = 'response')

switch_d0 <- readRDS(paste(data_direction,'/weight_model_switch_d0.rds', sep = ""))
switch_n0 <- readRDS(paste(data_direction,'/weight_model_switch_n0.rds', sep = ""))
switch_d1 <- readRDS(paste(data_direction,'/weight_model_switch_d1.rds', sep = ""))
switch_n1 <- readRDS(paste(data_direction,'/weight_model_switch_n1.rds', sep = ""))

design_mat <- expand.grid(id = 1:200,
                          trial_period = 0:4,
                          followup_time = 0:4) 
design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]

fitting_data_treatment <-  switch_data %>% 
  dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
  dplyr::select(id,trial_period, followup_time, X2,  X4, assigned_treatment) %>% 
  merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>% 
  dplyr::group_by(id) %>% 
  tidyr::fill( X2,X4,assigned_treatment,.direction = "down") %>% 
  dplyr::ungroup() %>% 
  dplyr::select(id, trial_period, followup_time, X2, X4, assigned_treatment) %>% 
  merge(data.frame(id = switch_data$id, trial_period = switch_data$trial_period), by = c("id", "trial_period"), all.y = TRUE) %>% 
  dplyr::arrange(id, trial_period, followup_time) %>% 
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
  dplyr::filter(trial_period == 0)

fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]

fitting_data_control <- fitting_data_treatment %>% 
  dplyr::mutate(assigned_treatment = assigned_treatment*0,
                t_1A = t_1*0,
                t_2A = t_2*0,
                t_3A = t_3*0,
                t_4A = t_4*0)

Y_pred_PP_treatment <- predict.glm(PP$model, 
                                   fitting_data_treatment, 
                                   type = "response")
Y_pred_PP_control <- predict.glm(PP$model, 
                                 fitting_data_control,
                                 type = "response")
predicted_probas_PP <- fitting_data_treatment %>% 
  dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment,
                predicted_proba_control = Y_pred_PP_control) %>% 
  dplyr::group_by(id, trial_period) %>% 
  dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(followup_time) %>% 
  dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                   survival_control = mean(cum_hazard_control),
                   survival_difference = survival_treatment - survival_control,
                   mrd = -survival_difference)

################################### Bootstrap C I ###################################

# Move the call to bootstrap_construct outside the for loop to avoid problems with random numbers 
boot_data <- list()
for (k in 1:bootstrap_iter) {
  boot_data[[k]] <- sort(sample(unique(switch_data$id), length(unique(switch_data$id)), replace = TRUE))
}
time <- proc.time()
############################# DIRECT BOOTSTRAP #############################
surv_PP_difference_boostrap_estimates <-as.data.frame(matrix(,5,bootstrap_iter))
surv_PP_difference_boostrap_estimates <- foreach(k = 1:bootstrap_iter, .combine=cbind) %dopar% {
  
  weights_table_boot <- data.frame(id = 1:200) %>% 
    rowwise() %>% 
    dplyr::mutate(weight_boot = length(boot_data[[k]][boot_data[[k]] == id])) #bootstrap weight is number of times they were sampled
  
  IP_model <- weight_func_bootstrap(data = simdata_censored, expanded_data = switch_data, 
                                    switch_d_cov = ~ X2 + X4,
                                    weight_model_d0 = switch_d0,
                                    weight_model_n0 = switch_n0,
                                    weight_model_d1 = switch_d1,
                                    weight_model_n1 = switch_n1,
                                    boot_idx = boot_data[[k]], remodel = TRUE, quiet = TRUE)
  
  #calculate IP weights from bootstrap sample
  
  boot_design_data <- IP_model$data %>%
    merge(weights_table_boot, by = 'id', all.y = TRUE) %>% 
    dplyr::mutate(weight = ifelse(weight_boot !=0,weight*weight_boot,0))
  
  #Direct bootstrap
  
  PP_boot <- TrialEmulation::trial_msm(data = boot_design_data,
                                            outcome_cov = ~ X2 + X4+ assigned_treatment+
                                              t_1 + t_2 + t_3 + t_4 +
                                              t_1A + t_2A + t_3A + t_4A + 
                                              t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
                                              t_1X4 + t_2X4 + t_3X4 + t_4X4,
                                            model_var = c('assigned_treatment'),
                                            glm_function = 'glm',
                                            include_trial_period  = ~1, include_followup_time = ~1,
                                            use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
  
  design_mat <- expand.grid(id = 1:tail(boot_design_data$id, n = 1), 
                            trial_period = 0:4,
                            followup_time = 0:4)
  design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]
  
  fitting_data_treatment_boot <- boot_design_data %>% 
    dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
    dplyr::select(id,trial_period, followup_time, X2,  X4, assigned_treatment) %>% 
    merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>% 
    dplyr::group_by(id) %>% 
    tidyr::fill( X2,X4,assigned_treatment,.direction = "down") %>% 
    dplyr::ungroup() %>% 
    dplyr::select(id, trial_period, followup_time, X2, X4, assigned_treatment) %>% 
    merge(data.frame(id = boot_design_data$id, trial_period = boot_design_data$trial_period), by = c("id", "trial_period"), all.y = TRUE) %>% 
    dplyr::arrange(id, trial_period, followup_time) %>% 
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
    dplyr::filter(trial_period == 0)
  
  
  fitting_data_treatment_boot <- fitting_data_treatment_boot[!duplicated(fitting_data_treatment_boot),]
  
  fitting_data_control_boot <- fitting_data_treatment_boot %>% 
    dplyr::mutate(assigned_treatment = assigned_treatment*0,
                  t_1A = t_1*0,
                  t_2A = t_2*0,
                  t_3A = t_3*0,
                  t_4A = t_4*0)
  
  Y_pred_PP_treatment_boot <- predict.glm(PP_boot$model, 
                                          fitting_data_treatment_boot, 
                                          type = "response")
  Y_pred_PP_control_boot <- predict.glm(PP_boot$model, 
                                        fitting_data_control_boot,
                                        type = "response")
  predicted_probas_PP_boot <- fitting_data_treatment_boot %>% 
    dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment_boot,
                  predicted_proba_control = Y_pred_PP_control_boot) %>% 
    dplyr::group_by(id, trial_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  
  predicted_probas_PP_boot[,3] - predicted_probas_PP_boot[,2]
  
}

bootstrap_CI <- cbind(apply(surv_PP_difference_boostrap_estimates,
                                                  1,
                                                  quantile,
                                                  probs = c(0.025)),
                      apply(surv_PP_difference_boostrap_estimates,
                                                  1,
                                                  quantile,
                                                  probs = c(0.975)))

################### LEF OUTCOME ONLY ##################################
X <- model.matrix(PP$model)
e <- PP$model$y - PP$model$fitted.values

surv_PP_difference_LEF_outcome_estimates <- as.data.frame(matrix(,5,bootstrap_iter))
surv_PP_difference_LEF_outcome_estimates <- foreach(k = 1:bootstrap_iter, .combine=cbind) %dopar% {
  
  weights_table_boot <- data.frame(id = 1:200) %>% 
    rowwise() %>% 
    dplyr::mutate(weight_boot = length(boot_data[[k]][boot_data[[k]] == id])) #bootstrap weight is number of times they were sampled
  
  IP_model <- weight_func_bootstrap(data = simdata_censored, expanded_data = switch_data, 
                                    switch_d_cov = ~ X2 + X4, 
                                    weight_model_d0 = switch_d0,
                                    weight_model_n0 = switch_n0,
                                    weight_model_d1 = switch_d1,
                                    weight_model_n1 = switch_n1,
                                    boot_idx = boot_data[[k]], remodel = TRUE, quiet = TRUE)
  
  #calculate IP weights from bootstrap sample
  
  boot_design_data <- IP_model$data %>%
    merge(weights_table_boot, by = 'id', all.y = TRUE) %>% 
    dplyr::mutate(weight = ifelse(weight_boot !=0,weight*weight_boot,0))
  
  LEFs <- t(X)%*%(boot_design_data$weight*e)
  LEFs[is.na(LEFs)] <- 0
  variance_mat <- vcov(PP$model)
  variance_mat[is.na(variance_mat)] <- 0
  #Calculate \hat \beta(b)
  beta <- PP$model$coefficients + variance_mat%*%LEFs
  PP_boot <- PP
  PP_boot$model$coefficients <- beta
  
  Y_pred_PP_treatment_boot <- predict.glm(PP_boot$model, 
                                          fitting_data_treatment, 
                                          type = "response")
  Y_pred_PP_control_boot <- predict.glm(PP_boot$model, 
                                        fitting_data_control,
                                        type = "response")
  predicted_probas_PP_boot <- fitting_data_treatment %>% 
    dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment_boot,
                  predicted_proba_control = Y_pred_PP_control_boot) %>% 
    dplyr::group_by(id, trial_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  predicted_probas_PP_boot[,3] - predicted_probas_PP_boot[,2]
}

LEF_outcome_CI<-cbind(apply(surv_PP_difference_LEF_outcome_estimates,
                                                    1,
                                                    quantile,
                                                    probs = c(0.025)),
                      apply(surv_PP_difference_LEF_outcome_estimates,
                                                     1,
                                                     quantile,
                                                     probs = c(0.975)))



#################### LEF WEIGHT AND OUTCOME  ##############################
X_sw_d0 <- model.matrix(switch_d0)
e_sw_d0 <- switch_d0$y - switch_d0$fitted.values

X_sw_n0 <- model.matrix(switch_n0)
e_sw_n0 <- switch_n0$y - switch_n0$fitted.values

X_sw_d1 <- model.matrix(switch_d1)
e_sw_d1 <- switch_d1$y - switch_d1$fitted.values

X_sw_n1 <- model.matrix(switch_n1)
e_sw_n1 <- switch_n1$y - switch_n1$fitted.values


surv_PP_difference_LEF_both_estimates <- as.data.frame(matrix(,5,bootstrap_iter))
surv_PP_difference_LEF_both_estimates <- foreach(k = 1:bootstrap_iter, .combine=cbind) %dopar% {
  weights_table_boot <- data.frame(id = 1:200) %>% 
    rowwise() %>% 
    dplyr::mutate(weight_boot = length(boot_data[[k]][boot_data[[k]] == id])) #bootstrap weight is number of times they were sampled
  
  data_0 <- merge(weights_table_boot, switch_d0$data, on = id, all.y = TRUE)
  data_1 <- merge(weights_table_boot, switch_d1$data, on = id, all.y = TRUE)
  
  LEF_sw_d0_boot <- t(X_sw_d0)%*%(data_0$weight_boot*e_sw_d0)
  LEF_sw_n0_boot <- t(X_sw_n0)%*%(data_0$weight_boot*e_sw_n0)
  LEF_sw_d1_boot <- t(X_sw_d1)%*%(data_1$weight_boot*e_sw_d1)
  LEF_sw_n1_boot <- t(X_sw_n1)%*%(data_1$weight_boot*e_sw_n1)
  
  #Calculate \hat \beta(b)
  beta_sw_d0 <- switch_d0$coefficients + vcov(switch_d0)%*%LEF_sw_d0_boot
  beta_sw_n0 <- switch_n0$coefficients + vcov(switch_n0)%*%LEF_sw_n0_boot
  beta_sw_d1 <- switch_d1$coefficients + vcov(switch_d1)%*%LEF_sw_d1_boot
  beta_sw_n1 <- switch_n1$coefficients + vcov(switch_n1)%*%LEF_sw_n1_boot
  
  IP_model <- weight_func_bootstrap(data = simdata_censored, expanded_data = switch_data, 
                                    switch_d_cov = ~ X2 + X4,
                                    weight_model_d0 = switch_d0,
                                    weight_model_n0 = switch_n0,
                                    weight_model_d1 = switch_d1,
                                    weight_model_n1 = switch_n1,
                                    new_coef_sw_d0 = beta_sw_d0,
                                    new_coef_sw_n0 = beta_sw_n0,
                                    new_coef_sw_d1 = beta_sw_d1,
                                    new_coef_sw_n1 = beta_sw_n1,
                                    boot_idx = boot_data[[k]], remodel = FALSE, quiet = TRUE)
  
  #calculate IP weights from bootstrap sample
  
  boot_design_data <- IP_model$data %>%
    merge(weights_table_boot, by = 'id', all.y = TRUE) %>% 
    dplyr::mutate(weight = ifelse(weight_boot !=0,weight*weight_boot,0))
  
  LEFs <- t(X)%*%(boot_design_data$weight*e)
  LEFs[is.na(LEFs)] <- 0
  variance_mat <- vcov(PP$model)
  variance_mat[is.na(variance_mat)] <- 0
  #Calculate \hat \beta(b)
  beta <- PP$model$coefficients + variance_mat%*%LEFs
  
  PP_boot <- PP$model
  PP_boot$coefficients <- beta
  
  Y_pred_PP_treatment_boot <- predict.glm(PP_boot, 
                                          fitting_data_treatment, 
                                          type = "response")
  Y_pred_PP_control_boot <- predict.glm(PP_boot, 
                                        fitting_data_control,
                                        type = "response")
  predicted_probas_PP_boot <- fitting_data_treatment %>% 
    dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment_boot,
                  predicted_proba_control = Y_pred_PP_control_boot) %>% 
    dplyr::group_by(id, trial_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  
  predicted_probas_PP_boot[,3] - predicted_probas_PP_boot[,2]
}
LEF_both_CI <-cbind(apply(surv_PP_difference_LEF_both_estimates,
                                                 1,
                                                 quantile,
                                                 probs = c(0.025)),
                    apply(surv_PP_difference_LEF_both_estimates,
                                                  1,
                                                  quantile,
                                                  probs = c(0.975)))

############################SANDWICH #######################################
covariance_mat <-PP$robust$matrix
if (all(eigen(covariance_mat)$values > 0) == F){
  not_pos_def <- not_pos_def + 1.0
  next
}
#Step 1 of algorithm  -- sampling Y_n1, ..., Y_nB ~ MN(coeffs,sandwich covariance)
sampling_size <- 200
coeffs_sample <- mvrnorm(sampling_size,PP$model$coefficients, covariance_mat)

surv_PP_difference_sandwich_estimates <- foreach(k = 1:sampling_size, .combine=cbind) %dopar% {
  
  #Step 1 of algorithm -- same model with new coeffs = one point from MVN sample
  fit_sample <- PP
  fit_sample$model$coefficients <- coeffs_sample[k,]
  
  #Step 2 -- calculating survival probas with new model
  Y_pred_sample_treatment <- predict.glm(fit_sample$model, 
                                         fitting_data_treatment, 
                                         type = "response")
  Y_pred_sample_control <- predict.glm(fit_sample$model, 
                                       fitting_data_control,
                                       type = "response")
  
  predicted_probas_PP_sample <- fitting_data_treatment %>% 
    dplyr::mutate(predicted_proba_treatment = Y_pred_sample_treatment,
                  predicted_proba_control = Y_pred_sample_control) %>% 
    dplyr::group_by(id, trial_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  
  predicted_probas_PP_sample[,3] - predicted_probas_PP_sample[,2]
}

#Step 3 -- calculating lower and upper bounds by 2.5% and 97.5% quantiles
sandwich_CI <- cbind(apply(surv_PP_difference_sandwich_estimates,
                                                  1,
                                                  quantile,
                                                  probs = c(0.025)),
                     apply(surv_PP_difference_sandwich_estimates,
                                                  1,
                                                  quantile,
                                                  probs = c(0.975)))

############################# JACKKNIFE #############################
# Move the call to bootstrap_construct outside the for loop to avoid problems with random numbers 
boot_data <- list()
lenid <- length(unique(switch_data$id))
for (k in 1:lenid) {
  boot_data[[k]] <- unique(switch_data$id[switch_data$id != k])
}

jackknife_est_theta<- array(, dim = c(length( PP$model$coefficients),lenid))
jackknife_est_mrd<- array(, dim = c(5,lenid))


for (k in 1:lenid){
  print(k)
  weights_table_boot <- data.frame(id = 1:200) %>% 
    rowwise() %>% 
    dplyr::mutate(weight_boot = length(boot_data[[k]][boot_data[[k]] == id])) #bootstrap weight is number of times they were sampled
  
  IP_model <- weight_func_bootstrap(data = simdata_censored, expanded_data = switch_data, 
                                    switch_d_cov = ~ X2 + X4,
                                    weight_model_d0 = switch_d0,
                                    weight_model_n0 = switch_n0,
                                    weight_model_d1 = switch_d1,
                                    weight_model_n1 = switch_n1,
                                    boot_idx = boot_data[[k]], remodel = TRUE, quiet = TRUE)
  
  #calculate IP weights from bootstrap sample
  
  boot_design_data <- IP_model$data %>%
    merge(weights_table_boot, by = 'id', all.y = TRUE) %>% 
    dplyr::mutate(weight = ifelse(weight_boot !=0,weight*weight_boot,0))
  
  #Direct bootstrap
  
  PP_boot <- TrialEmulation::trial_msm(data = boot_design_data,
                                       outcome_cov = ~ X2 + X4+ assigned_treatment+
                                         t_1 + t_2 + t_3 + t_4 +
                                         t_1A + t_2A + t_3A + t_4A + 
                                         t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
                                         t_1X4 + t_2X4 + t_3X4 + t_4X4,
                                       model_var = c('assigned_treatment'),
                                       glm_function = 'glm',
                                       include_trial_period = ~1, include_followup_time = ~1,
                                       use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
  
  jackknife_est_theta[,k] <-PP_boot$model$coefficients
  design_mat <- expand.grid(id = 1:tail(boot_design_data$id, n = 1), 
                              trial_period = 0:4,
                              followup_time = 0:4)
  design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]
  
  fitting_data_treatment_boot <- boot_design_data %>% 
    dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
    dplyr::select(id,trial_period, followup_time, X2,  X4, assigned_treatment) %>% 
    merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>% 
    dplyr::group_by(id) %>% 
    tidyr::fill( X2,X4,assigned_treatment,.direction = "down") %>% 
    dplyr::ungroup() %>% 
    dplyr::select(id, trial_period, followup_time, X2, X4, assigned_treatment) %>% 
    merge(data.frame(id = boot_design_data$id, trial_period = boot_design_data$trial_period), by = c("id", "trial_period"), all.y = TRUE) %>% 
    dplyr::arrange(id, trial_period, followup_time) %>% 
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
    dplyr::filter(trial_period == 0)
  
  
  fitting_data_treatment_boot <- fitting_data_treatment_boot[!duplicated(fitting_data_treatment_boot),]
  
  fitting_data_control_boot <- fitting_data_treatment_boot %>% 
    dplyr::mutate(assigned_treatment = assigned_treatment*0,
                  t_1A = t_1*0,
                  t_2A = t_2*0,
                  t_3A = t_3*0,
                  t_4A = t_4*0)
  
  Y_pred_PP_treatment_boot <- predict.glm(PP_boot$model, 
                                          fitting_data_treatment_boot, 
                                          type = "response")
  Y_pred_PP_control_boot <- predict.glm(PP_boot$model, 
                                        fitting_data_control_boot,
                                        type = "response")
  predicted_probas_PP_boot <- fitting_data_treatment_boot %>% 
    dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment_boot,
                  predicted_proba_control = Y_pred_PP_control_boot) %>% 
    dplyr::group_by(id, trial_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control),
                     mrd = survival_control-survival_treatment)
  jackknife_est_mrd[,k] <- pull(predicted_probas_PP_boot,mrd)
}

jackknife_theta_var <- array(0,dim = c(20,20))
for (k in 1:lenid){
  jackknife_theta_var <- jackknife_theta_var + outer(jackknife_est_theta[,k]- rowMeans(jackknife_est_theta[,]),jackknife_est_theta[,k]- rowMeans(jackknife_est_theta[,]))
}
jackknife_theta_var <- ((200-1)/200)*jackknife_theta_var

sampling_size <- 200
coeffs_sample <- mvrnorm(sampling_size,PP$model$coefficients, jackknife_theta_var)

surv_PP_difference_sandwich_estimates <- foreach(k = 1:sampling_size, .combine=cbind) %dopar% {
  
  #Step 1 of algorithm -- same model with new coeffs = one point from MVN sample
  fit_sample <- PP
  fit_sample$model$coefficients <- coeffs_sample[k,]
  
  #Step 2 -- calculating survival probas with new model
  Y_pred_sample_treatment <- predict.glm(fit_sample$model, 
                                         fitting_data_treatment, 
                                         type = "response")
  Y_pred_sample_control <- predict.glm(fit_sample$model, 
                                       fitting_data_control,
                                       type = "response")
  
  predicted_probas_PP_sample <- fitting_data_treatment %>% 
    dplyr::mutate(predicted_proba_treatment = Y_pred_sample_treatment,
                  predicted_proba_control = Y_pred_sample_control) %>% 
    dplyr::group_by(id, trial_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  
  predicted_probas_PP_sample[,3] - predicted_probas_PP_sample[,2]
}

jackknife_mrd_mvn_CI <- cbind(apply(surv_PP_difference_sandwich_estimates,
                                    1,
                                    quantile,
                                    probs = c(0.025)), apply(surv_PP_difference_sandwich_estimates,
                                                             1,
                                                             quantile,
                                                             probs = c(0.975)))



mrd <-pull(predicted_probas_PP,mrd)
jacknife_mrd_se <- array(,dim = c(5))
for (t in 1:5){
  jacknife_mrd_se[t] <- sqrt(((200-1)/200)*sum((jackknife_est_mrd[t,] - mean(jackknife_est_mrd[t,]))^2))
}

jackknife_mrd_CI <- cbind(mrd - 1.96*jacknife_mrd_se, mrd + 1.96*jacknife_mrd_se)

ggplot(data = predicted_probas_PP,aes(x = 0:4)) +
  geom_step(aes(x = 0:4,y = mrd)) +
  geom_stepribbon(aes(ymin = LEF_outcome_CI[,1],
                      ymax = LEF_outcome_CI[,2], color = "LEF outcome"), alpha = 0.1) +
  geom_stepribbon(aes(ymin = LEF_both_CI[,1],
                      ymax = LEF_both_CI[,2], color = "LEF both"), alpha = 0.1) +
  geom_stepribbon(aes(ymin = sandwich_CI[,1],
                      ymax = sandwich_CI[,2], color = "Sandwich"), alpha = 0.1) +
  geom_stepribbon(aes(ymin = bootstrap_CI[,1], 
                      ymax = bootstrap_CI[,2], color = "Nonparametric Boot."), alpha = 0.1) +
  geom_stepribbon(aes(ymin = jackknife_mrd_mvn_CI[,1], 
                      ymax = jackknife_mrd_mvn_CI[,2], color = "Jackknife MVN sampling"), alpha = 0.1) +
  geom_stepribbon(aes(ymin = jackknife_mrd_CI[,1], 
                      ymax = jackknife_mrd_CI[,2], color = "Jackknife Wald"), alpha = 0.1) +
  scale_color_manual(name = "CI method", values = c("Jackknife MVN sampling" = 'orange',"Jackknife Wald" = 'pink',"LEF outcome"= "green", "Nonparametric Boot." = "red", "LEF both" = 'purple', "Sandwich" = 'blue')) +
  labs(x = 'Follow-up time', 
       y = "Marginal risk difference")


