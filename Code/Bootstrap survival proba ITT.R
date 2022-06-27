## simulate data for testing TrialEmulation package, using the algorithm in Young and Tchetgen Tchetgen (2014) 
library(modelr)
library(reshape2)
library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/MPhil-dissertation/Code")
source("simulate_MSM.R")
set.seed(20222022)

simdata_censored<-DATA_GEN_censored(1000, 10) 


write.csv(simdata_censored, file= "MSM_censor.csv", row.names = F) 


library(RandomisedTrialsEmulation)

###################account for censoring##############################################

data_path = "MSM_censor.csv"

##### excluding obs before first becoming eligible, with censoring by dropout #######################################

#### ITT analysis

ITT <- RandomisedTrialsEmulation::initiators(data_path, id='ID', period='t', treatment='A', outcome='Y', eligible ='eligible', cense = 'C',
           outcomeCov_var=c('X3', 'X4', 'age_s'), outcomeCov =c( 'X3', 'X4', 'age_s'), model_var = c('assigned_treatment'),
           cov_censed = c( 'X1', 'X2','X3', 'X4', 'age_s'), model_censed =c( 'X1', 'X2','X3', 'X4', 'age_s'), pool_cense=1,
           include_expansion_time_case = 0,  include_followup_time_case = c("linear", "quadratic"), include_regime_length = 1,
           use_weight=1, use_censor=0, case_control = 0, data_dir =getwd(), numCores = 1, quiet = FALSE)

#### Survival function point estimate for ITT ####

design_mat <- expand.grid(id = 1:1000,
                          for_period = 0:9,
                          followup_time = 0:9) %>% 
  dplyr::mutate(followup_time2 = followup_time^2)
design_mat <- design_mat[which(10 -design_mat$for_period > design_mat$followup_time),]

switch_data <- read.csv("switch_data.csv")
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

write.csv(fitting_data_treatment, "fitting_data_treatment.csv")
write.csv(fitting_data_control, "fitting_data_control.csv")

Y_pred_ITT_treatment <- predict.glm(ITT$model$model, fitting_data_treatment, 
                                    type = "response")
Y_pred_ITT_control <- predict.glm(ITT$model$model, fitting_data_control, 
                                    type = "response")
predicted_probas_ITT <- fitting_data_treatment %>% 
  dplyr::mutate(predicted_proba_treatment = Y_pred_ITT_treatment,
                predicted_proba_control = Y_pred_ITT_control) %>% 
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

surv_ITT_treatment_boostrap_estimates <- as.data.frame(matrix(,10,bootstrap_iter))
surv_ITT_control_boostrap_estimates <-as.data.frame(matrix(,10,bootstrap_iter))
surv_ITT_difference_boostrap_estimates <-as.data.frame(matrix(,10,bootstrap_iter))
surv_ITT_ratio_boostrap_estimates <-as.data.frame(matrix(,10,bootstrap_iter))

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

for (i in 1:bootstrap_iter){
  boot_data <- bootstrap_construct(simdata_censored)
  write.csv(boot_data, "boot_data.csv",row.names = F)
  ITT_boot <- RandomisedTrialsEmulation::initiators("boot_data.csv", id='ID', period='t', treatment='A', outcome='Y', eligible ='eligible', cense = 'C',
                         outcomeCov_var=c('X3', 'X4', 'age_s'), outcomeCov =c( 'X3', 'X4', 'age_s'), model_var = c('assigned_treatment'),
                         cov_censed = c( 'X1', 'X2','X3', 'X4', 'age_s'), model_censed =c( 'X1', 'X2','X3', 'X4', 'age_s'), pool_cense=1,
                         include_expansion_time_case = 0,  include_followup_time_case = c("linear", "quadratic"), include_regime_length = 1,
                         use_weight=1, use_censor=0, case_control = 0, data_dir =getwd(), numCores = 4, quiet = TRUE)
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
  
  Y_pred_ITT_treatment_boot <- predict.glm(ITT_boot$model$model, 
                                                             fitting_data_treatment_boot, 
                                                            type = "response")
  Y_pred_ITT_control_boot <- predict.glm(ITT_boot$model$model, 
                                                           fitting_data_control_boot,
                                                           type = "response")
  predicted_probas_ITT_boot <- fitting_data_treatment_boot %>% 
    dplyr::mutate(predicted_proba_treatment = Y_pred_ITT_treatment_boot,
                  predicted_proba_control = Y_pred_ITT_control_boot) %>% 
    dplyr::group_by(id, for_period) %>% 
    dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                  cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(followup_time) %>% 
    dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                     survival_control = mean(cum_hazard_control))
  surv_ITT_treatment_boostrap_estimates[,i] <- predicted_probas_ITT_boot[,2]
  surv_ITT_control_boostrap_estimates[,i] <- predicted_probas_ITT_boot[,3]
  surv_ITT_difference_boostrap_estimates[,i] <- predicted_probas_ITT_boot[,2] - predicted_probas_ITT_boot[,3]
  surv_ITT_ratio_boostrap_estimates[,i] <- predicted_probas_ITT_boot[,2]/predicted_probas_ITT_boot[,3]
  
}
surv_ITT_treatment_boostrap_estimates$lb <- apply(surv_ITT_treatment_boostrap_estimates,
                                                    1,
                                                    quantile,
                                                    probs = c(0.025))
surv_ITT_treatment_boostrap_estimates$ub <- apply(surv_ITT_treatment_boostrap_estimates,
                                                    1,
                                                    quantile,
                                                    probs = c(0.975))
surv_ITT_control_boostrap_estimates$lb <- apply(surv_ITT_control_boostrap_estimates,
                                                    1,
                                                    quantile,
                                                    probs = c(0.025))
surv_ITT_control_boostrap_estimates$ub <- apply(surv_ITT_control_boostrap_estimates,
                                                    1,
                                                    quantile,
                                                    probs = c(0.975))
surv_ITT_difference_boostrap_estimates$lb <- apply(surv_ITT_difference_boostrap_estimates,
                                                1,
                                                quantile,
                                                probs = c(0.025))
surv_ITT_difference_boostrap_estimates$ub <- apply(surv_ITT_difference_boostrap_estimates,
                                                1,
                                                quantile,
                                                probs = c(0.975))
surv_ITT_ratio_boostrap_estimates$lb <- apply(surv_ITT_ratio_boostrap_estimates,
                                                   1,
                                                   quantile,
                                                   probs = c(0.025))
surv_ITT_ratio_boostrap_estimates$ub <- apply(surv_ITT_ratio_boostrap_estimates,
                                                   1,
                                                   quantile,
                                                   probs = c(0.975))

##### Survival function values and CI #####
predicted_probas_ITT <- predicted_probas_ITT %>% 
  dplyr::mutate(survival_treatment_lb = surv_ITT_treatment_boostrap_estimates$lb,
                survival_treatment_ub = surv_ITT_treatment_boostrap_estimates$ub,
                survival_control_lb = surv_ITT_control_boostrap_estimates$lb,
                survival_control_ub = surv_ITT_control_boostrap_estimates$ub,
                survival_difference_lb = surv_ITT_difference_boostrap_estimates$lb,
                survival_difference_ub = surv_ITT_difference_boostrap_estimates$ub,
                survival_ratio_lb = surv_ITT_ratio_boostrap_estimates$lb,
                survival_ratio_ub = surv_ITT_ratio_boostrap_estimates$ub)

write.csv(predicted_probas_ITT, "predicted_probas_ITT.csv")
save(ITT,file = "ITT_fit.rda")

##### Plot survival curves #####
library(ggplot2)
library(pammtools)
ggplot(data = predicted_probas_ITT, aes(followup_time)) +
  geom_step(aes(y = survival_treatment, color = "Treatment")) +
  geom_stepribbon(aes(ymin = survival_treatment_lb, ymax = survival_treatment_ub), alpha = 0.1) +
  geom_step(aes(y = survival_control, color = "Control")) + 
  geom_stepribbon(aes(ymin = survival_control_lb, ymax = survival_control_ub), alpha = 0.3) +
  scale_color_manual(name = "Treatment assignment", values = c("Treatment"= "red", "Control" = "blue")) +
  labs(x = 'Follow-up time', y = "Survival function", color = "Assigned treatment", title = "ITT analysis") 

ggplot(data = predicted_probas_ITT, aes(followup_time)) +
  geom_step(aes(y = survival_ratio)) +
  geom_stepribbon(aes(ymin = survival_ratio_lb, 
                      ymax = survival_ratio_ub), alpha = 0.1) +
  labs(x = 'Follow-up time', 
       y = "Ratio of survival functions for each treatment type", title = "ITT analysis")

ggplot(data = predicted_probas_ITT, aes(followup_time)) +
  geom_step(aes(y = survival_difference)) +
  geom_stepribbon(aes(ymin = survival_difference_lb, 
                      ymax = survival_difference_ub), alpha = 0.1) +
  labs(x = 'Follow-up time', 
       y = "Difference between treatment and control survival functions", title = "ITT analysis")

