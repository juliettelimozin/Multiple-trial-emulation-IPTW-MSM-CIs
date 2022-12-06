library(modelr)
library(tidyverse)
library(tidyr)
source("simulate_MSM_simplified.R")
library(MASS)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(matlib)
library(pammtools)
library(RandomisedTrialsEmulation)

#### Generate overall weights and LEF U(\hat beta)
data <- DATA_GEN_censored_reduced(1000,5, treat_prev = 0)
PP_prep <- RandomisedTrialsEmulation::data_preparation(data, id='ID', period='t', treatment='A', outcome='Y', 
                                                       eligible ='eligible',cense = 'C',
                                                       switch_d_cov = ~X2 + X4,
                                                       outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                                       cense_d_cov = ~X2 + X4,
                                                       include_expansion_time_case = ~1, include_followup_time_case = ~1, 
                                                       include_regime_length = F,
                                                       use_weight=1, use_censor=1, quiet = T)
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
                                                use_weight=1, use_censor=1, quiet = T, use_sample_weights =  F)
switch_data$p_i <- predict.glm(PP$model, switch_data,type = 'response')


X <- model.matrix(PP$model)
e <- PP$model$y - PP$model$fitted.values
LEF <- t(X)%*%(switch_data$weight*e)

design_mat <- expand.grid(id = 1:1000,
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
                t_4A = t_4*0)

###### LEF Bootstrap ######
bootstrap_iter <- 200

LEFs <- array(,dim = c(20,bootstrap_iter))
betas <- array(,dim = c(20,bootstrap_iter))
surv_PP_treatment_boostrap_estimates <- as.data.frame(matrix(,5,bootstrap_iter))
surv_PP_control_boostrap_estimates <-as.data.frame(matrix(,5,bootstrap_iter))
surv_PP_difference_boostrap_estimates <-as.data.frame(matrix(,5,bootstrap_iter))

for (i in 1:bootstrap_iter){
  #bootstrap sample of ids
  idx <- sort(sample(unique(switch_data$id), length(unique(switch_data$id)), replace = TRUE))
  weights_table_boot <- data.frame(id = 1:1000) %>% 
    rowwise() %>% 
    mutate(weight_boot = length(idx[idx == id])) #bootstrap weight is number of times they were sampled
  
  boot_data <- data[which(data$ID %in% idx),]
  
  
  #calculate IP weights from bootstrap sample
  PP_prep_boot <- RandomisedTrialsEmulation::data_preparation(boot_data, id='ID', period='t', 
                       treatment='A', outcome='Y', 
                       eligible ='eligible',cense = 'C',
                       switch_d_cov = ~X2 + X4,
                       outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                       cense_d_cov = ~X2 + X4,
                       include_expansion_time_case = ~1, include_followup_time_case = ~1, 
                       include_regime_length = F,
                       use_weight=1, use_censor=1, quiet = T)
 
   boot_design_data <- PP_prep_boot$data %>%
    select(id, for_period, followup_time, weight) %>% 
    merge(switch_data, by = c('id', 'for_period','followup_time'), all.y = T) %>% 
    merge(weights_table_boot, by = 'id', all.y = TRUE) %>% 
    mutate(weight = ifelse(weight_boot !=0,weight.x*weight_boot,0))

  LEFs[,i] <- t(X)%*%(boot_design_data$weight*e)
  
  #Calculate \hat \beta(b)
  betas[,i] <- PP$model$coefficients + vcov(PP$model)%*%LEFs[,i]
  
  PP_boot <- RandomisedTrialsEmulation::data_modelling(data = boot_design_data,
                                                  outcome_cov = ~ X2 + X4+ assigned_treatment+
                                                    t_1 + t_2 + t_3 + t_4 +
                                                    t_1A + t_2A + t_3A + t_4A + 
                                                    t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
                                                    t_1X4 + t_2X4 + t_3X4 + t_4X4,
                                                  model_var = c('assigned_treatment'),
                                                  glm_function = 'glm',
                                                  include_expansion_time_case = ~1, include_followup_time_case = ~1,
                                                  use_weight=1, use_censor=1, quiet = T, use_sample_weights =  F)
  
  design_mat <- expand.grid(id = 1:1000,
                            for_period = 0:4,
                            followup_time = 0:4)
  design_mat <- design_mat[which(5 -design_mat$for_period > design_mat$followup_time),]
  
  fitting_data_treatment_boot <-  boot_design_data %>% 
    dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
    dplyr::select(id,for_period, followup_time, X2,  X4, assigned_treatment) %>% 
    merge(design_mat, by = c("id", "for_period", "followup_time"), all.y = TRUE) %>% 
    dplyr::group_by(id) %>% 
    tidyr::fill( X2,X4,assigned_treatment,.direction = "down") %>% 
    dplyr::ungroup() %>% 
    dplyr::select(id, for_period, followup_time, X2, X4, assigned_treatment) %>% 
    merge(data.frame(id = boot_design_data$id, for_period = boot_design_data$for_period), by = c("id", "for_period"), all.y = TRUE) %>% 
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
  
  fitting_data_treatment_boot <- fitting_data_treatment_boot[!duplicated(fitting_data_treatment_boot),]
  
  fitting_data_control_boot <- fitting_data_treatment_boot %>% 
    dplyr::mutate(assigned_treatment = assigned_treatment*0,
                  t_1A = t_1*0,
                  t_2A = t_2*0,
                  t_3A = t_3*0,
                  t_4A = t_4*0,)
  
  Y_pred_PP_treatment_boot <- predict.glm(PP_boot$model, 
                                          fitting_data_treatment_boot, 
                                          type = "response")
  Y_pred_PP_control_boot <- predict.glm(PP_boot$model, 
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
                     survival_control = mean(cum_hazard_control),
                     survival_difference = survival_treatment - survival_control)
  
  surv_PP_treatment_boostrap_estimates[,i] <- pull(predicted_probas_PP_boot, survival_treatment)
  surv_PP_control_boostrap_estimates[,i] <- pull(predicted_probas_PP_boot,survival_control)
  surv_PP_difference_boostrap_estimates[,i] <- pull(predicted_probas_PP_boot,survival_difference)
}

#Calculate LEF variance estimate 
var_boot <- outer(betas[,1] - PP$model$coefficients, betas[,1] - PP$model$coefficients)
for (i in 2:bootstrap_iter){
  var_boot <- var_boot + outer(betas[,i] - PP$model$coefficients, betas[,i] - PP$model$coefficients)
}
var_boot <- var_boot/bootstrap_iter

#########Sampling Y_n1, ..., Y_nB ~ MN(coeffs,LEF covariance)
sampling_size <- 200
coeffs_sample <- mvrnorm(sampling_size,PP$model$coefficients, var_boot)

surv_PP_treatment_LEF_estimates <- as.data.frame(matrix(,5,sampling_size))
surv_PP_control_LEF_estimates <- as.data.frame(matrix(,5,sampling_size))
surv_PP_difference_LEF_estimates <- as.data.frame(matrix(,5,sampling_size))

for (i in 1:sampling_size){
  #Step 1 of algorithm -- same model with new coeffs = one point from MVN sample
  fit_sample <- PP
  fit_sample$model$coefficients <- coeffs_sample[i,]
  
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
    dplyr::group_by(id, for_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  surv_PP_treatment_LEF_estimates[,i] <- predicted_probas_PP_sample[,2]
  surv_PP_control_LEF_estimates[,i] <- predicted_probas_PP_sample[,3]
  surv_PP_difference_LEF_estimates[,i] <- predicted_probas_PP_sample[,2] - predicted_probas_PP_sample[,3]
}

#Step 3 -- calculating lower and upper bounds by 2.5% and 97.5% quantiles
surv_PP_treatment_LEF_estimates$lb <-apply(surv_PP_treatment_LEF_estimates,
                                                1,
                                                quantile,
                                                probs = c(0.025))
surv_PP_treatment_LEF_estimates$ub <- apply(surv_PP_treatment_LEF_estimates,
                                                 1,
                                                 quantile,
                                                 probs = c(0.975))
surv_PP_control_LEF_estimates$lb <-apply(surv_PP_control_LEF_estimates,
                                              1,
                                              quantile,
                                              probs = c(0.025))
surv_PP_control_LEF_estimates$ub <- apply(surv_PP_control_LEF_estimates,
                                               1,
                                               quantile,
                                               probs = c(0.975))
surv_PP_difference_LEF_estimates$lb <- apply(surv_PP_difference_LEF_estimates,
                                                  1,
                                                  quantile,
                                                  probs = c(0.025))
surv_PP_difference_LEF_estimates$ub <- apply(surv_PP_difference_LEF_estimates,
                                                  1,
                                                  quantile,
                                                  probs = c(0.975))
######## DIrect BOOTSTRAP #################
surv_PP_treatment_boostrap_estimates$lb <- apply(surv_PP_treatment_boostrap_estimates,
                                                 1,
                                                 quantile,
                                                 probs = c(0.025))
surv_PP_treatment_boostrap_estimates$ub <- apply(surv_PP_treatment_boostrap_estimates,
                                                 1,
                                                 quantile,
                                                 probs = c(0.975))
surv_PP_control_boostrap_estimates$lb <- apply(surv_PP_control_boostrap_estimates,
                                               1,
                                               quantile,
                                               probs = c(0.025))
surv_PP_control_boostrap_estimates$ub <- apply(surv_PP_control_boostrap_estimates,
                                               1,
                                               quantile,
                                               probs = c(0.975))
surv_PP_difference_boostrap_estimates$lb <- apply(surv_PP_difference_boostrap_estimates,
                                                  1,
                                                  quantile,
                                                  probs = c(0.025))
surv_PP_difference_boostrap_estimates$ub <- apply(surv_PP_difference_boostrap_estimates,
                                                  1,
                                                  quantile,
                                                  probs = c(0.975))


#### Survival function point estimate for PP ####

#Step 2 -- calculating survival probas with new model
Y_pred_treatment <- predict.glm(PP$model, 
                                       fitting_data_treatment, 
                                       type = "response")
Y_pred_control <- predict.glm(PP$model, 
                                     fitting_data_control,
                                     type = "response")

predicted_probas_PP <- fitting_data_treatment %>% 
  dplyr::mutate(predicted_proba_treatment = Y_pred_treatment,
                predicted_proba_control = Y_pred_control) %>% 
  dplyr::group_by(id, for_period) %>% 
  dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(followup_time) %>% 
  dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                   survival_control = mean(cum_hazard_control),
                   survival_difference = survival_treatment - survival_control)

#Update dataframe with results from bootstrap LEF CI
predicted_probas_PP <- predicted_probas_PP %>% 
  dplyr::mutate(survival_treatment_lb_LEF = surv_PP_treatment_LEF_estimates$lb,
                survival_treatment_ub_LEF = surv_PP_treatment_LEF_estimates$ub,
                survival_control_lb_LEF = surv_PP_control_LEF_estimates$lb,
                survival_control_ub_LEF = surv_PP_control_LEF_estimates$ub,
                survival_difference_lb_LEF = surv_PP_difference_LEF_estimates$lb,
                survival_difference_ub_LEF = surv_PP_difference_LEF_estimates$ub,
                survival_treatment_lb_boostrap = surv_PP_treatment_boostrap_estimates$lb,
                survival_treatment_ub_boostrap = surv_PP_treatment_boostrap_estimates$ub,
                survival_control_lb_boostrap = surv_PP_control_boostrap_estimates$lb,
                survival_control_ub_boostrap = surv_PP_control_boostrap_estimates$ub,
                survival_difference_lb_boostrap = surv_PP_difference_boostrap_estimates$lb,
                survival_difference_ub_boostrap = surv_PP_difference_boostrap_estimates$ub)


##### Plot survival curves #####
library(ggplot2)
library(pammtools)

ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_treatment, color = "Treatment")) +
  geom_stepribbon(aes(ymin = survival_treatment_lb_LEF, ymax = survival_treatment_ub_LEF), alpha = 0.1) +
  geom_step(aes(y = survival_control, color = "Control")) + 
  geom_stepribbon(aes(ymin = survival_control_lb_LEF, ymax = survival_control_ub_LEF), alpha = 0.3) +
  scale_color_manual(name = "Treatment assignment", values = c("Treatment"= "red", "Control" = "blue")) +
  labs(x = 'Follow-up time', y = "Survival function", color = "Assigned treatment", 
       title = "PP analysis - LEF CIs")

ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_difference)) +
  geom_stepribbon(aes(ymin = survival_difference_lb_LEF, 
                      ymax = survival_difference_ub_LEF, color = "LEF variance"), alpha = 0.1) +
  geom_stepribbon(aes(ymin = survival_difference_lb_boostrap, 
                      ymax = survival_difference_ub_boostrap, color = "Direct Boot."), alpha = 0.3) +
  scale_color_manual(name = "CI method", values = c("LEF variance"= "red", "Direct Boot." = "blue")) +
  labs(x = 'Follow-up time', 
       y = "Marginal risk difference", title = "PP analysis- CIs")

################### LEF Percentile method ####################

surv_PP_treatment_LEF_point_estimates <- as.data.frame(matrix(,5,sampling_size))
surv_PP_control_LEF_point_estimates <- as.data.frame(matrix(,5,sampling_size))
surv_PP_difference_LEF_point_estimates <- as.data.frame(matrix(,5,sampling_size))

for (i in 1:bootstrap_iter){
  #Step 1 of algorithm -- same model with new coeffs = one point from MVN sample
  fit_sample <- PP
  fit_sample$model$coefficients <- betas[,i]
  
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
    dplyr::group_by(id, for_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  surv_PP_treatment_LEF_point_estimates[,i] <- predicted_probas_PP_sample[,2]
  surv_PP_control_LEF_point_estimates[,i] <- predicted_probas_PP_sample[,3]
  surv_PP_difference_LEF_point_estimates[,i] <- predicted_probas_PP_sample[,2] - predicted_probas_PP_sample[,3]
}

surv_PP_treatment_LEF_point_estimates$lb <-apply(surv_PP_treatment_LEF_point_estimates,
                                           1,
                                           quantile,
                                           probs = c(0.025))
surv_PP_treatment_LEF_point_estimates$ub <- apply(surv_PP_treatment_LEF_point_estimates,
                                            1,
                                            quantile,
                                            probs = c(0.975))
surv_PP_control_LEF_point_estimates$lb <-apply(surv_PP_control_LEF_point_estimates,
                                         1,
                                         quantile,
                                         probs = c(0.025))
surv_PP_control_LEF_point_estimates$ub <- apply(surv_PP_control_LEF_point_estimates,
                                          1,
                                          quantile,
                                          probs = c(0.975))
surv_PP_difference_LEF_point_estimates$lb <- apply(surv_PP_difference_LEF_point_estimates,
                                             1,
                                             quantile,
                                             probs = c(0.025))
surv_PP_difference_LEF_point_estimates$ub <- apply(surv_PP_difference_LEF_point_estimates,
                                             1,
                                             quantile,
                                             probs = c(0.975))

predicted_probas_PP <- predicted_probas_PP %>% 
  dplyr::mutate(survival_treatment_lb_LEF_point = surv_PP_treatment_LEF_point_estimates$lb,
                survival_treatment_ub_LEF_point = surv_PP_treatment_LEF_point_estimates$ub,
                survival_control_lb_LEF_point = surv_PP_control_LEF_point_estimates$lb,
                survival_control_ub_LEF_point = surv_PP_control_LEF_point_estimates$ub,
                survival_difference_lb_LEF_point = surv_PP_difference_LEF_point_estimates$lb,
                survival_difference_ub_LEF_point = surv_PP_difference_LEF_point_estimates$ub)

ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_treatment, color = "Treatment")) +
  geom_stepribbon(aes(ymin = survival_treatment_lb_LEF_point, ymax = survival_treatment_ub_LEF_point), alpha = 0.1) +
  geom_step(aes(y = survival_control, color = "Control")) + 
  geom_stepribbon(aes(ymin = survival_control_lb_LEF_point, ymax = survival_control_ub_LEF_point), alpha = 0.3) +
  scale_color_manual(name = "Treatment assignment", values = c("Treatment"= "red", "Control" = "blue")) +
  labs(x = 'Follow-up time', y = "Survival function", color = "Assigned treatment", 
       title = "PP analysis - LEF percentile CIs")

ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_difference)) +
  geom_stepribbon(aes(ymin = survival_difference_lb_LEF_point, 
                      ymax = survival_difference_ub_LEF_point, color = "LEF percentile"), alpha = 0.1) +
  geom_stepribbon(aes(ymin = survival_difference_lb_boostrap, 
                      ymax = survival_difference_ub_boostrap, color = "Direct Boot."), alpha = 0.3) +
  scale_color_manual(name = "CI method", values = c("LEF percentile"= "red", "Direct Boot." = "blue")) +
  labs(x = 'Follow-up time', 
       y = "Marginal risk difference", title = "PP analysis- CIs")
