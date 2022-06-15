## simulate data for testing TrialEmulation package, using the algorithm in Young and Tchetgen Tchetgen (2014) 
library(modelr)
library(reshape2)
library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/MPhil-dissertation/Code")
source("simulate_MSM.R")
set.seed(20222022)

#RUN ITT FILE FIRST TO ENSURE SIMULATED DATA ARE THE SAME

library(RandomisedTrialsEmulation, lib.loc='/home/li/lib/R/R_LIBS/')

###################account for censoring##############################################

data_path = "MSM_censor.csv"

##### excluding obs before first becoming eligible, with censoring by dropout #######################################

#### PP analysis

PP <- RandomisedTrialsEmulation::initiators(data_path, id='ID', period='t', treatment='A', outcome='Y', eligible ='eligible', cense = 'C',
                                            model_switchd =c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                            cov_switchd = c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                            outcomeCov_var=c( 'X3', 'X4', 'age_s'), outcomeCov =c('X3', 'X4', 'age_s'), model_var = c('assigned_treatment'),
                                            cov_censed = c( 'X1', 'X2','X3', 'X4', 'age_s'), model_censed =c( 'X1', 'X2','X3', 'X4', 'age_s'), pool_cense=1,
                                            include_expansion_time_case = 0, include_followup_time_case = c("linear", "quadratic"), include_regime_length = 1,
                                            use_weight=1, use_censor=1, case_control = 0, data_dir =getwd(), numCores = 1, quiet = FALSE)

#### Survival function point estimate for PP ####

fitting_data_treatment <- read.csv("fitting_data_treatment.csv")
fitting_data_control <- read.csv("fitting_data_control.csv")

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
                   survival_control = mean(cum_hazard_control),
                   survival_difference = survival_treatment - survival_control,
                   survival_ratio = survival_treatment/survival_control)

###### Bootstrap CI ######
bootstrap_iter <- 100

surv_PP_treatment_boostrap_estimates <- as.data.frame(matrix(,10,bootstrap_iter))
surv_PP_control_boostrap_estimates <-as.data.frame(matrix(,10,bootstrap_iter))
surv_PP_difference_boostrap_estimates <-as.data.frame(matrix(,10,bootstrap_iter))
surv_PP_ratio_boostrap_estimates <-as.data.frame(matrix(,10,bootstrap_iter))

bootstrap_construct <- function(data){
  idx <- sort(sample(unique(data$ID), length(unique(data$ID)), replace = TRUE))
  data_boot <- data[which(data$ID == idx[1]),]
  for (i in 2:length(idx)){
    data_duplicate <- data[which(data$ID == idx[i]),]
    if (idx[i] == idx[i-1]){
      data_duplicate$ID <- replace(data_duplicate$ID, 
                                   data_duplicate$ID == idx[i], 
                                   tail(idx, n = 1) + i)
    }
    data_boot <- rbind(data_boot,data_duplicate)
  }
  return(data_boot)
}

for (i in 1:bootstrap_iter){
  boot_data <- bootstrap_construct(simdata_censored)
  write.csv(boot_data, "boot_data.csv",row.names = F)
  PP_boot <- RandomisedTrialsEmulation::initiators("boot_data.csv", id='ID', period='t', treatment='A', outcome='Y', eligible ='eligible', cense = 'C',
                                                   model_switchd =c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                                   cov_switchd = c( 'X1', 'X2', 'X3', 'X4', 'age_s'),
                                                   outcomeCov_var=c( 'X3', 'X4', 'age_s'), outcomeCov =c('X3', 'X4', 'age_s'), model_var = c('assigned_treatment'),
                                                   cov_censed = c( 'X1', 'X2','X3', 'X4', 'age_s'), model_censed =c( 'X1', 'X2','X3', 'X4', 'age_s'), pool_cense=1,
                                                   include_expansion_time_case = 0, include_followup_time_case = c("linear", "quadratic"), include_regime_length = 1,
                                                   use_weight=1, use_censor=1, case_control = 0, data_dir =getwd(), numCores = 1, quiet = TRUE)
  switch_data_boot <- read.csv("switch_data.csv")
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
  
  Y_pred_PP_treatment_boot <- predict.glm(PP_boot$model$model, 
                                           fitting_data_treatment_boot, 
                                           type = "response")
  Y_pred_PP_control_boot <- predict.glm(PP_boot$model$model, 
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
  
  surv_PP_treatment_boostrap_estimates[,i] <- predicted_probas_PP_boot[,2]
  surv_PP_control_boostrap_estimates[,i] <- predicted_probas_PP_boot[,3]
  surv_PP_difference_boostrap_estimates[,i] <- predicted_probas_PP_boot[,2] - predicted_probas_PP_boot[,3]
  surv_PP_ratio_boostrap_estimates[,i] <- predicted_probas_PP_boot[,2]/predicted_probas_PP_boot[,3]
}
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
surv_PP_ratio_boostrap_estimates$lb <- apply(surv_PP_ratio_boostrap_estimates,
                                              1,
                                              quantile,
                                              probs = c(0.025))
surv_PP_ratio_boostrap_estimates$ub <- apply(surv_PP_ratio_boostrap_estimates,
                                              1,
                                              quantile,
                                              probs = c(0.975))

##### Survival function values and CI #####
predicted_probas_PP <- predicted_probas_PP %>% 
  dplyr::mutate(survival_treatment_lb = surv_PP_treatment_boostrap_estimates$lb,
                survival_treatment_ub = surv_PP_treatment_boostrap_estimates$ub,
                survival_control_lb = surv_PP_control_boostrap_estimates$lb,
                survival_control_ub = surv_PP_control_boostrap_estimates$ub,
                survival_difference_lb = surv_PP_difference_boostrap_estimates$lb,
                survival_difference_ub = surv_PP_difference_boostrap_estimates$ub,
                survival_ratio_lb = surv_PP_ratio_boostrap_estimates$lb,
                survival_ratio_ub = surv_PP_ratio_boostrap_estimates$ub)

write.csv(predicted_probas_PP, "predicted_probas_PP.csv")
save(PP,file = "PP_fit.rda")

##### Plot survival curves #####
library(ggplot2)
library(pammtools)
ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_treatment, color = "Treatment")) +
  geom_stepribbon(aes(ymin = survival_treatment_lb, ymax = survival_treatment_ub), alpha = 0.1) +
  geom_step(aes(y = survival_control, color = "Control")) + 
  geom_stepribbon(aes(ymin = survival_control_lb, ymax = survival_control_ub), alpha = 0.3) +
  scale_color_manual(name = "Treatment assignment", values = c("Treatment"= "red", "Control" = "blue")) +
  labs(x = 'Follow-up time', y = "Survival function", color = "Assigned treatment", title = "PP Analysis")

ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_ratio)) +
  geom_stepribbon(aes(ymin = survival_ratio_lb, 
                      ymax = survival_ratio_ub), alpha = 0.1) +
  labs(x = 'Follow-up time', 
       y = "Ratio of survival functions for each treatment type", title = "PP analysis")

ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_difference)) +
  geom_stepribbon(aes(ymin = survival_difference_lb, 
                      ymax = survival_difference_ub), alpha = 0.1) +
  labs(x = 'Follow-up time', 
       y = "Difference between treatment and control survival functions", title = "PP analysis")

