#!/usr/bin R
library(modelr)
library(tidyverse)
library(tidyr)
setwd("~/rds/hpc-work/Project1")
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
sampling_size <- 500

size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)
scenarios <- tidyr::crossing(size,conf, treat)

CI_bootstrap_PP_red <- array(,dim = c(5,2,iters))
CI_sandwich_PP_red <- array(, dim = c(5,2,iters))
CI_LEF_outcome_PP_red <- array(, dim = c(5,2,iters))
CI_LEF_both_PP_red <- array(, dim = c(5,2,iters))
CI_jackknife_mvn_PP_red <- array(, dim = c(5,2,iters))
CI_jackknife_wald_PP_red <- array(, dim = c(5,2,iters))

computation_time <- array(,dim = c(6,iters))

estimates <- array(, dim = c(5,iters))
survival_treatment_estimates <- array(, dim = c(5,iters))
survival_control_estimates <- array(, dim = c(5,iters))

l <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

sandwich_mrd <- array(, dim = c(5,sampling_size,iters))
bootstrap_mrd <-array(, dim = c(5,bootstrap_iter, iters))
LEF_outcome_mrd <- array(, dim = c(5,bootstrap_iter, iters))
LEF_both_mrd <- array(, dim = c(5,bootstrap_iter, iters))
jackknife_wald_mrd <- array(, dim = c(5, as.numeric(scenarios[l,1]), iters))
jackknife_mvn_mrd <- array(, dim = c(5, sampling_size, iters))

coeff_dim <- array(,dim = c(iters))
bootstrap_nas <- array(0,dim = c(iters))
jackknife_nas <- array(0,dim = c(iters))
jackknife_SEs <- array(, dim = c(5,iters))

not_pos_def <- 0.0

data_direction <- paste("~/rds/hpc-work/Project1/models_scenario_high_",l,sep = "")
# Set number of cores. 67 is sufficient for 200 cores.
registerDoParallel(cores = 67)

for (i in 1:iters){
  tryCatch({
    print(i)
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
                       mrd = - survival_difference )
    
    estimates[,i] <- pull(predicted_probas_PP,mrd)
    survival_treatment_estimates[,i] <- pull(predicted_probas_PP,survival_treatment)
    survival_control_estimates[,i] <- pull(predicted_probas_PP,survival_control)
    
    ################################### Bootstrap C I ###################################
    
    # Move the call to bootstrap_construct outside the for loop to avoid problems with random numbers 
    boot_data <- list()
    for (k in 1:bootstrap_iter) {
      boot_data[[k]] <- sort(sample(unique(switch_data$id), length(unique(switch_data$id)), replace = TRUE))
    }
    print('Bootstrap')
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
      
      PP_boot <- TrialEmulation::trial_msm(data = boot_design_data,
                                           outcome_cov = my_covariates,
                                           model_var = c('assigned_treatment'),
                                           glm_function = 'glm',
                                           include_trial_period = ~1, include_followup_time = ~1,
                                           use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
      
      if(is.na(PP_boot$model$coefficients['t_4A']) == T){
        bootstrap_nas[i] < bootstrap_nas[i] + 1
        PP_boot$model <- update(PP_boot$model, . ~ . - t_4A, data = boot_design_data)
        if(is.na(PP_boot$model$coefficients['t_3A']) == T){
          PP_boot$model <- update(PP_boot$model, . ~ . - t_3A, data = boot_design_data)
        }
      }
      
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
    
    bootstrap_mrd[,,i] <- as.matrix(surv_PP_difference_boostrap_estimates)
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
    
    print('LEF outcome')
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
        dplyr::group_by(id, trial_period) %>% 
        dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                      cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
        dplyr::ungroup() %>% 
        dplyr::group_by(followup_time) %>% 
        dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                         survival_control = mean(cum_hazard_control))
      predicted_probas_PP_boot[,3] - predicted_probas_PP_boot[,2]
    }
    LEF_outcome_mrd[,,i] <- as.matrix(surv_PP_difference_LEF_outcome_estimates)
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
    print('LEF both')
    
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
        dplyr::group_by(id, trial_period) %>% 
        dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
                      cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
        dplyr::ungroup() %>% 
        dplyr::group_by(followup_time) %>% 
        dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
                         survival_control = mean(cum_hazard_control))
      
      predicted_probas_PP_boot[,3] - predicted_probas_PP_boot[,2]
    }
    LEF_both_mrd[,,i] <- as.matrix(surv_PP_difference_LEF_both_estimates)
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
    
    if(as.numeric(scenarios[l,1]) == 200){
      time <- proc.time()
      print('Jackknife')
      
      ############################# JACKKNIFE #############################
      # Move the call to bootstrap_construct outside the for loop to avoid problems with random numbers 
      boot_data <- list()
      lenid <- length(unique(switch_data$id))
      for (k in 1:lenid) {
        boot_data[[k]] <- unique(switch_data$id[switch_data$id != k])
      }
      
      jackknife_est_theta<- array(, dim = c(length(PP$model$coefficients),lenid))
      jackknife_est_mrd<- array(, dim = c(5,lenid))
      fail <- FALSE
      for (k in 1:lenid){
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
        
        PP_boot <- TrialEmulation::trial_msm(data = boot_design_data,
                                             outcome_cov = my_covariates,
                                             model_var = c('assigned_treatment'),
                                             glm_function = 'glm',
                                             include_trial_period = ~1, include_followup_time = ~1,
                                             use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
        if(is.na(PP_boot$model$coefficients['t_4A']) == T){
          fail <- TRUE
          jackknife_nas[i] <- jackknife_nas[i]+1
          PP_boot$model <- update(PP_boot$model, . ~ . - t_4A, data = boot_design_data)
          if(is.na(PP_boot$model$coefficients['t_3A']) == T){
            PP_boot$model <- update(PP_boot$model, . ~ . - t_3A, data = boot_design_data)
          }
        }
        
        if(length(PP_boot$model$coefficients) == length(PP$model$coefficients)){
          jackknife_est_theta[,k] <-PP_boot$model$coefficients
        }
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
      
      mrd <-pull(predicted_probas_PP,mrd)
      jacknife_mrd_se <- array(,dim = c(5))
      for (t in 1:5){
        jacknife_mrd_se[t] <- sqrt(((as.numeric(scenarios[l,1])-1)/as.numeric(scenarios[l,1]))*sum((jackknife_est_mrd[t,] - mrd[t])^2))
      }
      
      jackknife_SEs[,i]<- jacknife_mrd_se
      CI_jackknife_wald_PP_red[,,i] <- cbind(mrd - 1.96*jacknife_mrd_se, mrd + 1.96*jacknife_mrd_se)
      
      computation_time[4,i] <- (proc.time() - time)[[3]]
      jackknife_wald_mrd[,,i]<- jackknife_est_mrd
      
      if(fail == FALSE){
        jackknife_theta_var <- array(0,dim = c(length(PP$model$coefficients),length(PP$model$coefficients)))
        b_bar <- array(0,dim = c(length(PP$model$coefficients)))
        for (k in 1:lenid){
          b_bar <- b_bar + (as.numeric(scenarios[l,1])*PP$model$coefficients - (as.numeric(scenarios[l,1])-1)*jackknife_est_theta[,k])/as.numeric(scenarios[l,1])
        }
        for (k in 1:lenid){
          bi <- as.numeric(scenarios[l,1])*PP$model$coefficients - (as.numeric(scenarios[l,1])-1)*jackknife_est_theta[,k]
          jackknife_theta_var <- jackknife_theta_var + (1/((as.numeric(scenarios[l,1])-1)*as.numeric(scenarios[l,1])))*outer(bi- b_bar, bi-b_bar)
        }
        
        jackknife_theta_var[is.na(jackknife_theta_var)] <- 0
        
        coeffs_sample <- mvrnorm(sampling_size,PP$model$coefficients, jackknife_theta_var)
        
        surv_PP_difference_jackknife_estimates <- foreach(k = 1:sampling_size, .combine=cbind) %dopar% {
          
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
        
        CI_jackknife_mvn_PP_red[,,i] <- cbind(apply(surv_PP_difference_jackknife_estimates,
                                                    1,
                                                    quantile,
                                                    probs = c(0.025)), 
                                              apply(surv_PP_difference_jackknife_estimates,
                                                    1,
                                                    quantile,
                                                    probs = c(0.975)))
        
        computation_time[5,i] <- (proc.time() - time)[[3]]
        jackknife_mvn_mrd[,,i] <- as.matrix(surv_PP_difference_jackknife_estimates)
      }
    }
    print('Sandwich')
    
    time <- proc.time()
    ############################SANDWICH #######################################
    covariance_mat <-PP$robust$matrix
    if (all(eigen(covariance_mat)$values > 0) == F){
      not_pos_def <- not_pos_def + 1.0
      next
    }
    
    
    coeffs_sample <- mvrnorm(sampling_size,PP$model$coefficients[!is.na(PP$model$coefficients)], covariance_mat)
    
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
    sandwich_mrd[,,i]<- as.matrix(surv_PP_difference_sandwich_estimates)
    
    #Step 3 -- calculating lower and upper bounds by 2.5% and 97.5% quantiles
    surv_PP_difference_sandwich_estimates$lb <- apply(surv_PP_difference_sandwich_estimates,
                                                      1,
                                                      quantile,
                                                      probs = c(0.025))
    surv_PP_difference_sandwich_estimates$ub <- apply(surv_PP_difference_sandwich_estimates,
                                                      1,
                                                      quantile,
                                                      probs = c(0.975))
    
    
    computation_time[6,i] <- (proc.time() - time)[[3]]
    CI_sandwich_PP_red[,1,i] <- surv_PP_difference_sandwich_estimates$lb
    CI_sandwich_PP_red[,2,i] <- surv_PP_difference_sandwich_estimates$ub
    
    print('All done')
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
print(paste0("% not pos def: ", not_pos_def*100/iters))
save(CI_bootstrap_PP_red, file = paste("NewSimusJ/J_CI_bootstrap_PP_red_high_",as.character(l),".rda", sep = ""))
save(CI_sandwich_PP_red, file = paste("NewSimusJ/J_CI_sandwich_PP_red_high_",as.character(l),".rda", sep = ""))
save(CI_LEF_outcome_PP_red, file = paste("NewSimusJ/J_CI_LEF_outcome_PP_red_high_",as.character(l),".rda", sep = ""))
save(CI_LEF_both_PP_red, file = paste("NewSimusJ/J_CI_LEF_both_PP_red_high_",as.character(l),".rda", sep = ""))
save(CI_jackknife_mvn_PP_red, file = paste("NewSimusJ/J_CI_jackknife_mvn_PP_red_high_",as.character(l),".rda", sep = ""))
save(CI_jackknife_wald_PP_red, file = paste("NewSimusJ/J_CI_jackknife_wald_PP_red_high_",as.character(l),".rda", sep = ""))
save(computation_time, file = paste("NewSimusJ/J_computation_time_high_",as.character(l),".rda", sep = ""))
save(estimates, file = paste("NewSimusJ/J_estimates_red_high_", as.character(l), ".rda", sep = ""))
save(survival_treatment_estimates, file = paste("NewSimusJ/J_survival_treatment_estimates_high_", as.character(l), ".rda", sep = ""))
save(survival_control_estimates, file = paste("NewSimusJ/J_survival_control_estimates_high_", as.character(l), ".rda", sep = ""))
save(jackknife_SEs, file = paste("NewSimusJ/J_jackknife_SEs_high_", as.character(l), ".rda", sep = ""))
save(bootstrap_mrd, file = paste("NewSimusJ/J_bootstrap_mrd_high_", as.character(l), ".rda", sep = "")) 
save(LEF_outcome_mrd, file = paste("NewSimusJ/J_LEF_outcome_mrd_high_", as.character(l), ".rda", sep = "")) 
save(LEF_both_mrd, file = paste("NewSimusJ/J_LEF_both_mrd_high_", as.character(l), ".rda", sep = "")) 
save(jackknife_wald_mrd, file = paste("NewSimusJ/J_jackknife_wald_mrd_high_", as.character(l), ".rda", sep = "")) 
save(jackknife_mvn_mrd, file = paste("NewSimusJ/J_jackknife_mvn_mrd_high_", as.character(l), ".rda", sep = "")) 
save(sandwich_mrd, file = paste("NewSimusJ/J_sandwich_mrd_high_", as.character(l), ".rda", sep = "")) 
save(coeff_dim, file = paste("NewSimusJ/coeff_dim_high_", as.character(l), ".rda", sep = "")) 
save(bootstrap_nas, file = paste("NewSimusJ/bootstrap_nas_high_", as.character(l), ".rda", sep = "")) 
save(jackknife_nas, file = paste("NewSimusJ/jackknife_nas_high_", as.character(l), ".rda", sep = "")) 



