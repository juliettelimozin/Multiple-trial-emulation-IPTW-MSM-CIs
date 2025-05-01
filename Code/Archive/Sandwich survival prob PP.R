## simulate data for testing TrialEmulation package, using the algorithm in Young and Tchetgen Tchetgen (2014) 
library(modelr)
library(reshape2)
library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/MPhil-dissertation/Code")
source("simulate_MSM.R")
set.seed(20222022)
library(MASS)
library(RandomisedTrialsEmulation)
library(matrixStats)
library(Metrics)
library(sandwich)


#RUN BOOTSTRAP FILES TO WORK ON SAME SIMULATED DATA AND MODEL FIT
load("PP_fit.rda")
predicted_probas_PP <- read.csv("predicted_probas_PP.csv")
fitting_data_treatment <- read.csv("fitting_data_treatment.csv")
fitting_data_control <- read.csv("fitting_data_control.csv")

#Calculate robust SE's covariance matrix 

covariance_mat <- PP$model$robust$matrix

#Step 1 of algorithm  -- sampling Y_n1, ..., Y_nB ~ MN(coeffs,sandwich covariance)
sampling_size <- 200
coeffs_sample <- mvrnorm(sampling_size,PP$model$model$coefficients, covariance_mat)

surv_PP_treatment_sandwich_estimates <- as.data.frame(matrix(,10,sampling_size))
surv_PP_control_sandwich_estimates <- as.data.frame(matrix(,10,sampling_size))
surv_PP_difference_sandwich_estimates <- as.data.frame(matrix(,10,sampling_size))
surv_PP_ratio_sandwich_estimates <- as.data.frame(matrix(,10,sampling_size))

for (i in 1:sampling_size){
  #Step 1 of algorithm -- same model with new coeffs = one point from MVN sample
  fit_sample <- PP
  fit_sample$model$model$coefficients <- coeffs_sample[i,]
  
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
  surv_PP_treatment_sandwich_estimates[,i] <- predicted_probas_PP_sample[,2]
  surv_PP_control_sandwich_estimates[,i] <- predicted_probas_PP_sample[,3]
  surv_PP_difference_sandwich_estimates[,i] <- predicted_probas_PP_sample[,2] - predicted_probas_PP_sample[,3]
  surv_PP_ratio_sandwich_estimates[,i] <- predicted_probas_PP_sample[,2]/predicted_probas_PP_sample[,3]
}

#Step 3 -- calculating lower and upper bounds by 2.5% and 97.5% quantiles
surv_PP_treatment_sandwich_estimates$lb <-apply(surv_PP_treatment_sandwich_estimates,
                                                 1,
                                                 quantile,
                                                 probs = c(0.025))
surv_PP_treatment_sandwich_estimates$ub <- apply(surv_PP_treatment_sandwich_estimates,
                                                  1,
                                                  quantile,
                                                  probs = c(0.975))
surv_PP_control_sandwich_estimates$lb <-apply(surv_PP_control_sandwich_estimates,
                                               1,
                                               quantile,
                                               probs = c(0.025))
surv_PP_control_sandwich_estimates$ub <- apply(surv_PP_control_sandwich_estimates,
                                                1,
                                                quantile,
                                                probs = c(0.975))
surv_PP_difference_sandwich_estimates$lb <- apply(surv_PP_difference_sandwich_estimates,
                                                   1,
                                                   quantile,
                                                   probs = c(0.025))
surv_PP_difference_sandwich_estimates$ub <- apply(surv_PP_difference_sandwich_estimates,
                                                   1,
                                                   quantile,
                                                   probs = c(0.975))

surv_PP_ratio_sandwich_estimates$lb <- apply(surv_PP_ratio_sandwich_estimates,
                                              1,
                                              quantile,
                                              probs = c(0.025))
surv_PP_ratio_sandwich_estimates$ub <- apply(surv_PP_ratio_sandwich_estimates,
                                              1,
                                              quantile,
                                              probs = c(0.975))


#Update dataframe with results
predicted_probas_PP <- predicted_probas_PP %>% 
  dplyr::mutate(survival_treatment_lb_sandwich = surv_PP_treatment_sandwich_estimates$lb,
                survival_treatment_ub_sandwich = surv_PP_treatment_sandwich_estimates$ub,
                survival_control_lb_sandwich = surv_PP_control_sandwich_estimates$lb,
                survival_control_ub_sandwich = surv_PP_control_sandwich_estimates$ub,
                survival_difference_lb_sandwich = surv_PP_difference_sandwich_estimates$lb,
                survival_difference_ub_sandwich = surv_PP_difference_sandwich_estimates$ub,
                survival_ratio_lb_sandwich = surv_PP_ratio_sandwich_estimates$lb,
                survival_ratio_ub_sandwich = surv_PP_ratio_sandwich_estimates$ub)


##### Plot survival curves #####
library(ggplot2)
library(pammtools)

ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_treatment, color = "Treatment")) +
  geom_stepribbon(aes(ymin = survival_treatment_lb, ymax = survival_treatment_ub), alpha = 0.1) +
  geom_step(aes(y = survival_control, color = "Control")) + 
  geom_stepribbon(aes(ymin = survival_control_lb, ymax = survival_control_ub), alpha = 0.3) +
  scale_color_manual(name = "Treatment assignment", values = c("Treatment"= "red", "Control" = "blue")) +
  labs(x = 'Follow-up time', y = "Survival function", color = "Assigned treatment", 
       title = "PP analysis - bootstrap CIs") 
ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_treatment, color = "Treatment")) +
  geom_stepribbon(aes(ymin = survival_treatment_lb_sandwich, ymax = survival_treatment_ub_sandwich), alpha = 0.1) +
  geom_step(aes(y = survival_control, color = "Control")) + 
  geom_stepribbon(aes(ymin = survival_control_lb_sandwich, ymax = survival_control_ub_sandwich), alpha = 0.3) +
  scale_color_manual(name = "Treatment assignment", values = c("Treatment"= "red", "Control" = "blue")) +
  labs(x = 'Follow-up time', y = "Survival function", color = "Assigned treatment", 
       title = "PP analysis - sandwich CIs")

ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_ratio)) +
  geom_stepribbon(aes(ymin = survival_ratio_lb_sandwich, 
                      ymax = survival_ratio_ub_sandwich), alpha = 0.1) +
  labs(x = 'Follow-up time', 
       y = "Ratio of survival functions for each treatment type", title = "PP analysis - sandwich CIs")

ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_difference)) +
  geom_stepribbon(aes(ymin = survival_difference_lb_sandwich, 
                      ymax = survival_difference_ub_sandwich), alpha = 0.1) +
  labs(x = 'Follow-up time', 
       y = "Difference between treatment and control survival functions", title = "PP analysis- sandwich CIs")

