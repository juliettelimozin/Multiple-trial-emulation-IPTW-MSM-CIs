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

data_direction <- paste("~/rds/hpc-work/models_scenario_med_",l,sep = "")
# Set number of cores. 67 is sufficient for 200 cores.
registerDoParallel(cores = 67)

for (i in 1:iters){
  tryCatch({
    simdata_censored<-DATA_GEN_censored_reduced(as.numeric(scenarios[l,1]), 5, 
                                                conf = as.numeric(scenarios[l,2]), 
                                                treat_prev = as.numeric(scenarios[l,3]),
                                                outcome_prev = -3.8,
                                                censor = F)
    PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                                eligible ='eligible',
                                                switch_d_cov = ~X2 + X4,
                                                outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                                include_regime_length = F,
                                                use_weight=1, use_censor=1, quiet = T,
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
                                         use_weight=1, use_censor=1, quiet = T, use_sample_weights =  F)
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
      
      weights_table_boot <- data.frame(id = 1:as.numeric(scenarios[l,1])) %>% 
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
      
      PP_boot <- TrialEmulation::data_modelling(data = boot_design_data,
                                                outcome_cov = ~ X2 + X4+ assigned_treatment+
                                                  t_1 + t_2 + t_3 + t_4 +
                                                  t_1A + t_2A + t_3A + t_4A + 
                                                  t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
                                                  t_1X4 + t_2X4 + t_3X4 + t_4X4,
                                                model_var = c('assigned_treatment'),
                                                glm_function = 'glm',
                                                include_expansion_time = ~1, include_followup_time = ~1,
                                                use_weight=1, use_censor=1, quiet = T, use_sample_weights =  F)
      
      design_mat <- expand.grid(id = 1:tail(boot_design_data$id, n = 1), 
                                for_period = 0:4,
                                followup_time = 0:4)
      design_mat <- design_mat[which(5 -design_mat$for_period > design_mat$followup_time),]
      
      fitting_data_treatment_boot <- boot_design_data %>% 
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
        dplyr::group_by(id, for_period) %>% 
        dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                      cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
        dplyr::ungroup() %>% 
        dplyr::group_by(followup_time) %>% 
        dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                         survival_control = mean(cum_hazard_control))
      
      predicted_probas_PP_boot[,2] - predicted_probas_PP_boot[,3]
      
    }
    
    surv_PP_difference_boostrap_estimates$lb <- apply(surv_PP_difference_boostrap_estimates,
                                                      1,
                                                      quantile,
                                                      probs = c(0.025))
    surv_PP_difference_boostrap_estimates$ub <- apply(surv_PP_difference_boostrap_estimates,
                                                      1,
                                                      quantile,
                                                      probs = c(0.975))
    
    computation_time[1,i] <- (proc.time() - time)[[3]]
    CI_bootstrap_PP_red[,1,i] <- surv_PP_difference_boostrap_estimates$lb
    CI_bootstrap_PP_red[,2,i] <- surv_PP_difference_boostrap_estimates$ub
    
    time <- proc.time()
    ################### LEF OUTCOME ONLY ##################################
    X <- model.matrix(PP$model)
    e <- PP$model$y - PP$model$fitted.values
    
    surv_PP_difference_LEF_outcome_estimates <- as.data.frame(matrix(,5,bootstrap_iter))
    surv_PP_difference_LEF_outcome_estimates <- foreach(k = 1:bootstrap_iter, .combine=cbind) %dopar% {
      
      weights_table_boot <- data.frame(id = 1:as.numeric(scenarios[l,1])) %>% 
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
        dplyr::group_by(id, for_period) %>% 
        dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                      cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
        dplyr::ungroup() %>% 
        dplyr::group_by(followup_time) %>% 
        dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                         survival_control = mean(cum_hazard_control))
      predicted_probas_PP_boot[,2] - predicted_probas_PP_boot[,3]
    }
    surv_PP_difference_LEF_outcome_estimates$lb <-apply(surv_PP_difference_LEF_outcome_estimates,
                                                        1,
                                                        quantile,
                                                        probs = c(0.025))
    surv_PP_difference_LEF_outcome_estimates$ub <- apply(surv_PP_difference_LEF_outcome_estimates,
                                                         1,
                                                         quantile,
                                                         probs = c(0.975))
    
    
    computation_time[2,i] <- (proc.time() - time)[[3]]
    CI_LEF_outcome_PP_red[,1,i] <- surv_PP_difference_LEF_outcome_estimates$lb
    CI_LEF_outcome_PP_red[,2,i] <- surv_PP_difference_LEF_outcome_estimates$ub
    
    time <- proc.time()
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
      weights_table_boot <- data.frame(id = 1:as.numeric(scenarios[l,1])) %>% 
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
        dplyr::group_by(id, for_period) %>% 
        dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                      cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
        dplyr::ungroup() %>% 
        dplyr::group_by(followup_time) %>% 
        dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                         survival_control = mean(cum_hazard_control))
      
      predicted_probas_PP_boot[,2] - predicted_probas_PP_boot[,3]
    }
    surv_PP_difference_LEF_both_estimates$lb <-apply(surv_PP_difference_LEF_both_estimates,
                                                     1,
                                                     quantile,
                                                     probs = c(0.025))
    surv_PP_difference_LEF_both_estimates$ub <- apply(surv_PP_difference_LEF_both_estimates,
                                                      1,
                                                      quantile,
                                                      probs = c(0.975))
    
    
    
    computation_time[3,i] <- (proc.time() - time)[[3]]
    CI_LEF_both_PP_red[,1,i] <- surv_PP_difference_LEF_both_estimates$lb
    CI_LEF_both_PP_red[,2,i] <- surv_PP_difference_LEF_both_estimates$ub
    
    
    time <- proc.time()
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
        dplyr::group_by(id, for_period) %>% 
        dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                      cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
        dplyr::ungroup() %>% 
        dplyr::group_by(followup_time) %>% 
        dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                         survival_control = mean(cum_hazard_control))
      
      predicted_probas_PP_sample[,2] - predicted_probas_PP_sample[,3]
    }
    
    #Step 3 -- calculating lower and upper bounds by 2.5% and 97.5% quantiles
    surv_PP_difference_sandwich_estimates$lb <- apply(surv_PP_difference_sandwich_estimates,
                                                      1,
                                                      quantile,
                                                      probs = c(0.025))
    surv_PP_difference_sandwich_estimates$ub <- apply(surv_PP_difference_sandwich_estimates,
                                                      1,
                                                      quantile,
                                                      probs = c(0.975))
    
    
    computation_time[4,i] <- (proc.time() - time)[[3]]
    CI_sandwich_PP_red[,1,i] <- surv_PP_difference_sandwich_estimates$lb
    CI_sandwich_PP_red[,2,i] <- surv_PP_difference_sandwich_estimates$ub
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
print(paste0("% not pos def: ", not_pos_def*100/iters))
save(CI_bootstrap_PP_red, file = paste("CI_bootstrap_PP_red_med_",as.character(l),".rda", sep = ""))
save(CI_sandwich_PP_red, file = paste("CI_sandwich_PP_red_med_",as.character(l),".rda", sep = ""))
save(CI_LEF_outcome_PP_red, file = paste("CI_LEF_outcome_PP_red_med_",as.character(l),".rda", sep = ""))
save(CI_LEF_both_PP_red, file = paste("CI_LEF_both_PP_red_med_",as.character(l),".rda", sep = ""))
save(computation_time, file = paste("computation_time_med_",as.character(l),".rda", sep = ""))
save(estimates, file = paste("estimates_red_med_", as.character(l), ".rda", sep = ""))

