#!/usr/bin R

library(modelr)
library(tidyverse)
library(tidyr)
setwd("~/rds/hpc-work/Project1")
source("~/rds/hpc-work/Multiple-trial-emulation-IPTW-MSM-CIs/Code/simulate_MSM_simplified.R")
source("~/rds/hpc-work/Multiple-trial-emulation-IPTW-MSM-CIs/Code/weight_func.R")
set.seed(20250228)
library(TrialEmulation, lib.loc = '/home/jml219/R/x86_64-redhat-linux-gnu-library/4.3')
library(MASS)
library(sandwich)
library(doParallel)
library(doRNG)
library(rlist)
iters <- 50
bootstrap_iter <- 500
sampling_size <- 500

size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)
scenarios <- tidyr::crossing(size,conf, treat)

estimates <- array(, dim = c(5,iters))

l <- 1

coeff_dim <- array(,dim = c(iters))

data_direction <- paste("~/rds/hpc-work/Project1/models_scenario_low_",l,sep = "")
# Set number of cores. 67 is sufficient for 200 cores.
registerDoParallel(cores = 67)

for (i in 1:iters){
  tryCatch({
    print(i)
    simdata_censored<-DATA_GEN_censored_reduced(as.numeric(scenarios[l,1]), 5, 
                                                conf = as.numeric(scenarios[l,2]), 
                                                treat_prev = as.numeric(scenarios[l,3]),
                                                outcome_prev = -4.7,
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
    
    my_covariates <- ~ X2 + X4+ assigned_treatment+
      t_1 + t_2 + t_3 + t_4 +
      t_1A + t_2A + t_3A + t_4A + 
      t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
      t_1X4 + t_2X4 + t_3X4 + t_4X4
    
    PP <- TrialEmulation::trial_msm(data = switch_data,
                                    outcome_cov = my_covariates,
                                    model_var = c('assigned_treatment'),
                                    glm_function = 'glm',
                                    include_trial_period = ~1, include_followup_time = ~1,
                                    use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
    coeff_dim[i] <- dim(PP$robust$matrix)[1]
    
    if(is.na(PP$model$coefficients['t_4A']) == T){
      PP$model <- update(PP$model, . ~ . - t_4A, data = switch_data)
      my_covariates <- ~ X2 + X4+ assigned_treatment+
        t_1 + t_2 + t_3 + t_4 +
        t_1A + t_2A + t_3A + 
        t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
        t_1X4 + t_2X4 + t_3X4 + t_4X4
      if(is.na(PP$model$coefficients['t_3A']) == T){
        PP$model <- update(PP$model, . ~ . - t_3A, data = switch_data)
        my_covariates <- ~ X2 + X4+ assigned_treatment+
          t_1 + t_2 + t_3 + t_4 +
          t_1A + t_2A + 
          t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
          t_1X4 + t_2X4 + t_3X4 + t_4X4
      }
    }
    
    switch_data$p_i <- predict.glm(PP$model, switch_data,type = 'response')
    
    switch_d0 <- readRDS(paste(data_direction,'/weight_model_switch_d0.rds', sep = ""))
    switch_n0 <- readRDS(paste(data_direction,'/weight_model_switch_n0.rds', sep = ""))
    switch_d1 <- readRDS(paste(data_direction,'/weight_model_switch_d1.rds', sep = ""))
    switch_n1 <- readRDS(paste(data_direction,'/weight_model_switch_n1.rds', sep = ""))
    
    design_mat <- expand.grid(id = 1:as.numeric(scenarios[l,1]),
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
    
    fitting_data_treatment <- fitting_data_treatment %>% distinct()
    
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
                       mrd = - survival_difference )
    
    estimates[,i] <- pull(predicted_probas_PP,mrd)
  })
}

save(estimates, file = paste("~/rds/hpc-work/Multiple-trial-emulation-IPTW-MSM-CIs/Code/esimates_test_2_low_", as.character(l), ".rda", sep = ""))
  