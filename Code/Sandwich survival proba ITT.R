## simulate data for testing TrialEmulation package, using the algorithm in Young and Tchetgen Tchetgen (2014) 
library(modelr)
library(reshape2)
library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/MPhil-dissertation/Code")
source("simulate_MSM.R")
set.seed(20222022)
library(MASS)
library(RandomisedTrialsEmulation, lib.loc='/home/li/lib/R/R_LIBS/')
library(matrixStats)
library(Metrics)


#RUN BOOTSTRAP FILES TO WORK ON SAME SIMULATED DATA AND MODEL FIT
load("ITT_fit.rda")
predicted_probas_ITT <- read.csv("predicted_probas_ITT.csv")
fitting_data_treatment <- read.csv("fitting_data_treatment.csv")
fitting_data_control <- read.csv("fitting_data_control.csv")

#Step 1 of algorithm  -- sampling Y_n1, ..., Y_nB ~ MN(0,covmat(coeffs))
sampling_size <- 100
coeffs_sample <- mvrnorm(sampling_size,c(0,0,0,0,0,0,0), vcov(ITT$model$model))

surv_ITT_treatment_sandwich_estimates <- as.data.frame(matrix(,10,sampling_size))
surv_ITT_control_sandwich_estimates <- as.data.frame(matrix(,10,sampling_size))
surv_ITT_difference_sandwich_estimates <- as.data.frame(matrix(,10,sampling_size))
surv_ITT_ratio_sandwich_estimates <- as.data.frame(matrix(,10,sampling_size))

for (i in 1:sampling_size){
  #Step 1 of algorithm -- calculating theta*_nb = coeffs + B^-1/2 * Y_nb
  fit_sample <- ITT
  fit_sample$model$model$coefficients <-ITT$model$model$coefficients + coeffs_sample[i,]*sampling_size^(-1/2)
  
  #Step 2 -- calculating g(theta*)
  Y_pred_sample_treatment <- predict.glm(fit_sample$model$model, 
                                           fitting_data_treatment, 
                                           type = "response")
  Y_pred_sample_control <- predict.glm(fit_sample$model$model, 
                                         fitting_data_control,
                                         type = "response")
  
  predicted_probas_ITT_sample <- fitting_data_treatment %>% 
    dplyr::mutate(predicted_proba_treatment = Y_pred_sample_treatment,
                  predicted_proba_control = Y_pred_sample_control) %>% 
    dplyr::group_by(id, for_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  surv_ITT_treatment_sandwich_estimates[,i] <- predicted_probas_ITT_sample[,2]
  surv_ITT_control_sandwich_estimates[,i] <- predicted_probas_ITT_sample[,3]
  surv_ITT_difference_sandwich_estimates[,i] <- predicted_probas_ITT_sample[,2] - predicted_probas_ITT_sample[,3]
  surv_ITT_ratio_sandwich_estimates[,i] <- predicted_probas_ITT_sample[,2]/predicted_probas_ITT_sample[,3]
}

#Step 3 -- calculating g(theta*) +- 1.96*MSE(g(theta*))
surv_ITT_treatment_sandwich_estimates <- (surv_ITT_treatment_sandwich_estimates - predicted_probas_ITT$survival_treatment)^2
surv_ITT_treatment_sandwich_estimates$MSE <- rowMeans(as.matrix(surv_ITT_treatment_sandwich_estimates))

surv_ITT_control_sandwich_estimates <- (surv_ITT_control_sandwich_estimates - predicted_probas_ITT$survival_control)^2
surv_ITT_control_sandwich_estimates$MSE <- rowMeans(as.matrix(surv_ITT_control_sandwich_estimates))

surv_ITT_difference_sandwich_estimates <- (surv_ITT_difference_sandwich_estimates - predicted_probas_ITT$survival_difference)^2
surv_ITT_difference_sandwich_estimates$MSE <- rowMeans(as.matrix(surv_ITT_difference_sandwich_estimates))

surv_ITT_ratio_sandwich_estimates <- (surv_ITT_ratio_sandwich_estimates - predicted_probas_ITT$survival_ratio)^2
surv_ITT_ratio_sandwich_estimates$MSE <- rowMeans(as.matrix(surv_ITT_ratio_sandwich_estimates))


#Update dataframe with results
predicted_probas_ITT <- predicted_probas_ITT %>% 
  dplyr::mutate(survival_treatment_lb_sandwich = survival_treatment - 1.96 * (surv_ITT_treatment_sandwich_estimates$MSE)^(1/2),
                survival_treatment_ub_sandwich = survival_treatment + 1.96 * (surv_ITT_treatment_sandwich_estimates$MSE)^(1/2),
                survival_control_lb_sandwich = survival_control - 1.96 * (surv_ITT_control_sandwich_estimates$MSE)^(1/2),
                survival_control_ub_sandwich = survival_control + 1.96 * (surv_ITT_control_sandwich_estimates$MSE)^(1/2),
                survival_difference_lb_sandwich = survival_difference - 1.96 * (surv_ITT_difference_sandwich_estimates$MSE)^(1/2),
                survival_difference_ub_sandwich = survival_difference + 1.96 * (surv_ITT_difference_sandwich_estimates$MSE)^(1/2),
                survival_ratio_lb_sandwich = survival_ratio - 1.96 * (surv_ITT_ratio_sandwich_estimates$MSE)^(1/2),
                survival_ratio_lb_sandwich = survival_ratio + 1.96 * (surv_ITT_ratio_sandwich_estimates$MSE)^(1/2)
                )

##### Plot survival curves #####
library(ggplot2)
library(pammtools)

ggplot(data = predicted_probas_ITT, aes(followup_time)) +
  geom_step(aes(y = survival_treatment, color = "Treatment")) +
  geom_stepribbon(aes(ymin = survival_treatment_lb, ymax = survival_treatment_ub), alpha = 0.1) +
  geom_step(aes(y = survival_control, color = "Control")) + 
  geom_stepribbon(aes(ymin = survival_control_lb, ymax = survival_control_ub), alpha = 0.3) +
  scale_color_manual(name = "Treatment assignment", values = c("Treatment"= "red", "Control" = "blue")) +
  labs(x = 'Follow-up time', y = "Survival function", color = "Assigned treatment", 
       title = "ITT analysis - bootstrap CIs") +
  ylim(0.96,1)
ggplot(data = predicted_probas_ITT, aes(followup_time)) +
  geom_step(aes(y = survival_treatment, color = "Treatment")) +
  geom_stepribbon(aes(ymin = survival_treatment_lb_sandwich, ymax = survival_treatment_ub_sandwich), alpha = 0.1) +
  geom_step(aes(y = survival_control, color = "Control")) + 
  geom_stepribbon(aes(ymin = survival_control_lb_sandwich, ymax = survival_control_ub_sandwich), alpha = 0.3) +
  scale_color_manual(name = "Treatment assignment", values = c("Treatment"= "red", "Control" = "blue")) +
  labs(x = 'Follow-up time', y = "Survival function", color = "Assigned treatment", 
       title = "ITT analysis - sandwich CIs") +
  ylim(0.96,1)

ggplot(data = predicted_probas_ITT, aes(followup_time)) +
  geom_step(aes(y = survival_ratio)) +
  geom_stepribbon(aes(ymin = survival_ratio_lb_sandwich, 
                      ymax = survival_ratio_ub_sandwich), alpha = 0.1) +
  labs(x = 'Follow-up time', 
       y = "Ratio of survival functions for each treatment type", title = "ITT analysis")

ggplot(data = predicted_probas_ITT, aes(followup_time)) +
  geom_step(aes(y = survival_difference)) +
  geom_stepribbon(aes(ymin = survival_difference_lb_sandwich, 
                      ymax = survival_difference_ub_sandwich), alpha = 0.1) +
  labs(x = 'Follow-up time', 
       y = "Difference between treatment and control survival functions", title = "ITT analysis")
  
  
  
  