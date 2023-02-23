#!/usr/bin R
library(modelr)
library(tidyverse)
library(tidyr)
setwd("~/rds/hpc-work")
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
l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)

scenarios <- tidyr::crossing(size,conf, treat)

CI_bootstrap_PP_red <- array(,dim = c(5,2,iters))
CI_sandwich_PP_red <- array(, dim = c(5,2,iters))
CI_LEF_outcome_PP_red <- array(, dim = c(5,2,iters))
CI_LEF_both_PP_red <- array(, dim = c(5,2,iters))

computation_time <- array(,dim = c(4,iters))

estimates <- array(, dim = c(5,iters))

l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

not_pos_def <- 0.0

data_direction <- paste("~/rds/hpc-work/models_scenario_",l,sep = "")
# Set number of cores. 67 is sufficient for 200 cores.
registerDoParallel(cores = 67)

for (i in 1:iters){
  tryCatch({
    simdata_censored<-DATA_GEN_censored_reduced(as.numeric(scenarios[l,1]), 5, 
                                                conf = as.numeric(scenarios[l,2]), 
                                                treat_prev = as.numeric(scenarios[l,3]),
                                                outcome_prev = -3,
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
    
    PP <- TrialEmulation::data_modelling(data = switch_data,
                                         outcome_cov = ~ X2 + X4+ assigned_treatment+
                                           t_1 + t_2 + t_3 + t_4 +
                                           t_1A + t_2A + t_3A + t_4A + 
                                           t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
                                           t_1X4 + t_2X4 + t_3X4 + t_4X4,
                                         model_var = c('assigned_treatment'),
                                         glm_function = 'glm',
                                         include_expansion_time = ~1, include_followup_time = ~1,
                                         use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
    switch_data$p_i <- predict.glm(PP$model, switch_data,type = 'response')
    
    switch_d0 <- readRDS(paste(data_direction,'/weight_model_switch_d0.rds', sep = ""))
    switch_n0 <- readRDS(paste(data_direction,'/weight_model_switch_n0.rds', sep = ""))
    switch_d1 <- readRDS(paste(data_direction,'/weight_model_switch_d1.rds', sep = ""))
    switch_n1 <- readRDS(paste(data_direction,'/weight_model_switch_n1.rds', sep = ""))
    
    design_mat <- expand.grid(id = 1:as.numeric(scenarios[l,1]),
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
    
    ########## ESTIMATES ####################
    Y_pred_PP_treatment <- predict.glm(PP$model, 
                                       fitting_data_treatment, 
                                       type = "response")
    Y_pred_PP_control <- predict.glm(PP$model, 
                                     fitting_data_control,
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
                       survival_difference = survival_treatment - survival_control)
    
    estimates[,i] <- pull(predicted_probas_PP,survival_difference)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
save(estimates, file = paste("estimates_red_high_", as.character(l), ".rda", sep = ""))

