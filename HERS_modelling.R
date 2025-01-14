load('hers.Rdata')
library(modelr)
library(reshape2)
library(tidyverse)
library(tidyr)
source("simulate_MSM_simplified.R")
source("weight_func.R")
data_direction <- getwd()
library(TrialEmulation, lib.loc = '/home/jml219/R/x86_64-redhat-linux-gnu-library/4.3')
library(MASS)
library(sandwich)
library(doParallel)
library(doRNG)
library(pammtools)
library(lmtest)
library(xtable)

############# DATA PREPARATION ##################
HERS$Y <- as.factor(HERS$Y)
HERS$t <- HERS$visit - 8
HERS$SITE1 <- as.factor(HERS$SITE1)
HERS$SITE2 <- as.factor(HERS$SITE2)
HERS$SITE3 <- as.factor(HERS$SITE3)
HERS$WHITE <- as.factor(HERS$WHITE)
HERS$OTHER <- as.factor(HERS$OTHER)

HERS$CD4 <- (sqrt(as.numeric(HERS$CD4)) - mean(sqrt(HERS$CD4)))/sd(sqrt(HERS$CD4))
HERS$CD4_1 <- (sqrt(as.numeric(HERS$CD4_1)) - mean(sqrt(HERS$CD4_1)))/sd(sqrt(HERS$CD4_1))
HERS$CD4_2 <- (sqrt(as.numeric(HERS$CD4_2)) - mean(sqrt(HERS$CD4_2)))/sd(sqrt(HERS$CD4_2))

# HERS$HIVsym <- as.factor(HERS$HIVsym)
# HERS$HIVsym_1 <- as.factor(HERS$HIVsym_1)
# HERS$HIVsym_2 <- as.factor(HERS$HIVsym_2)

HERS$viral <- (log10(HERS$viral) - mean(log10(HERS$viral)))/sd(log10(HERS$viral))
HERS$viral_1 <- (log10(HERS$viral_1) - mean(log10(HERS$viral_1)))/sd(log10(HERS$viral_1))
HERS$viral_2 <- (log10(HERS$viral_2) - mean(log10(HERS$viral_2)))/sd(log10(HERS$viral_2))

HERS <- HERS %>% 
  dplyr::arrange(id,t) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(CAp = cumsum(as.numeric(haart_1 == 1 | haart_2 == 1))) %>% 
  dplyr::ungroup()

HERS$A <- HERS$haart
HERS$Ap <- HERS$haart_1
HERS$App <- HERS$haart_2
HERS[,'ID'] <- HERS$id
HERS <- HERS %>% 
  dplyr::mutate(SITE = as.factor(ifelse(SITE1 == 1, 1, ifelse(SITE2 == 1,2,3))),
                ETHNICITY = as.factor(ifelse(WHITE == 0 & OTHER == 0, 0, ifelse(WHITE == 1, 1, 2)))) %>% 
  dplyr::select(ID, t, A, Ap, App,CAp, CD4, CD4_1,CD4_2,
         viral,viral_1,viral_2,HIVsym,HIVsym_1,HIVsym_2,
         SITE, WHITE, OTHER,ETHNICITY, Y, C)
HERS$eligible <- as.numeric(HERS$CAp == 0)

summary_HERS <- HERS %>% 
  summarize(total_cense = sum(as.numeric(Y)-1))
#################SEQUENTIAL TRIALS IPW AND MSM #########################
PP_prep <- TrialEmulation::data_preparation(data = HERS, id='ID', period='t', treatment='A', outcome='Y', 
                                            eligible ='eligible', cense = 'C',
                                            switch_d_cov = ~ CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE + WHITE + OTHER,
                                            cense_d_cov = ~ CD4_1 + CD4  + viral + viral_1  + HIVsym + HIVsym_1+ SITE + WHITE + OTHER,
                                            outcome_cov = ~ CD4 + CD4_1 + CD4_2 + viral+ viral_1 + viral_2 + HIVsym + HIVsym_1 + HIVsym_2+ SITE + WHITE + OTHER,
                                            model_var = c('assigned_treatment'),
                                            use_weight=TRUE, use_censor=TRUE, quiet = F,
                                            save_weight_models = T, 
                                            data_dir = data_direction)
switch_data <- PP_prep$data %>% 
  dplyr::mutate(haartCD4_1 = assigned_treatment*CD4_1)

summary_trial0 <- switch_data %>% 
  dplyr::filter(trial_period == 4,followup_time == 0) %>% 
  dplyr::count(assigned_treatment)

summary_trial0 <- HERS %>% 
  dplyr::group_by(t) %>% 
  dplyr::count(A)

PP <- TrialEmulation::trial_msm(data = switch_data,
                                     outcome_cov = ~ CD4_1 + CD4_2 + viral_1 + viral_2 + SITE + WHITE + OTHER + assigned_treatment+
                                       haartCD4_1,
                                     model_var = c('assigned_treatment'),
                                     glm_function = 'glm',
                                    include_trial_period = ~1, include_followup_time = ~1,
                                     use_weight=T, use_censor=T, quiet = F, use_sample_weights =  F)

print(xtable(as.data.frame(summary(PP$model)$coefficients), type = "latex", digits = 3))

switch_d0 <- readRDS(paste(data_direction,'/weight_model_switch_d0.rds', sep = ""))
switch_n0 <- readRDS(paste(data_direction,'/weight_model_switch_n0.rds', sep = ""))
switch_d1 <- readRDS(paste(data_direction,'/weight_model_switch_d1.rds', sep = ""))
switch_n1 <- readRDS(paste(data_direction,'/weight_model_switch_n1.rds', sep = ""))

cense_d0 <- readRDS(paste(data_direction,'/cense_model_d0.rds', sep = ""))
cense_d1 <- readRDS(paste(data_direction,'/cense_model_d1.rds', sep = ""))
cense_n0 <- readRDS(paste(data_direction,'/cense_model_n0.rds', sep = ""))
cense_n1 <- readRDS(paste(data_direction,'/cense_model_n1.rds', sep = ""))

#################### TRADITIONAL IPW-MSM ('SINGLE TRIAL') ####################################
# weight_model_switch_d1 <- glm(A ~ CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE + WHITE + OTHER,
#                               family = binomial(link = "logit"), data = HERS[HERS$Ap == 1])
# weight_model_switch_d0 <- glm(A ~ CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE + WHITE + OTHER,
#                               family = binomial(link = "logit"), data = HERS[HERS$Ap == 0])
# weight_model_switch_n1 <- glm(A ~ 1,
#                               family = binomial(link = "logit"), data = HERS[HERS$Ap == 1])
# weight_model_switch_n0 <- glm(A ~ 1,
#                               family = binomial(link = "logit"), data = HERS[HERS$Ap == 0])
# weight_model_censor_d1 <- glm(C ~ CD4_1 + CD4 + viral + viral_1 + HIVsym + HIVsym_1+ SITE + WHITE + OTHER,
#                                family = binomial(link = "logit"), data = HERS[HERS$Ap == 1])
# weight_model_censor_d0 <- glm(C ~ CD4_1 + CD4 + viral + viral_1 + HIVsym + HIVsym_1+ SITE + WHITE + OTHER,
#                               family = binomial(link = "logit"), data = HERS[HERS$Ap == 0])
# weight_model_censor_n1 <- glm(C ~ 1,
#                              family = binomial(link = "logit"), data = HERS[HERS$Ap == 1])
# weight_model_censor_n0 <- glm(C ~ 1,
#                               family = binomial(link = "logit"), data = HERS[HERS$Ap == 0])
# 
# 
# hers_data_tradi <- HERS %>% 
#   dplyr::arrange(ID, t) %>% 
#   dplyr::mutate(weight_A = ifelse(Ap == 1 & A == 1, 
#                                   predict.glm(weight_model_switch_n1,HERS[HERS$Ap == 1],
#                                               type = 'response'),
#                                   ifelse(Ap == 1 & A == 0,
#                                          1-predict.glm(weight_model_switch_n1,HERS[HERS$Ap == 1], type = 'response'),
#                                          ifelse(Ap ==0 & A == 1,
#                                                 predict.glm(weight_model_switch_n0,HERS[HERS$Ap == 0],
#                                                             type = 'response'),
#                                                 1-predict.glm(weight_model_switch_n0,HERS[HERS$Ap == 0],
#                                                               type = 'response'))))/ifelse(Ap == 1 & A == 1, 
#                                     predict.glm(weight_model_switch_d1,HERS[HERS$Ap == 1],
#                                                 type = 'response'),
#                                     ifelse(Ap == 1 & A == 0,
#                                            1-predict.glm(weight_model_switch_d1,HERS[HERS$Ap == 1], type = 'response'),
#                                            ifelse(Ap ==0 & A == 1,
#                                                   predict.glm(weight_model_switch_d0,HERS[HERS$Ap == 0],
#                                                                               type = 'response'),
#                                                   1-predict.glm(weight_model_switch_d0,HERS[HERS$Ap == 0],
#                                                               type = 'response')))),
#                 weight_C = ifelse(Ap == 1 & C == 1, 
#                                   predict.glm(weight_model_censor_n1,HERS[HERS$Ap == 1],
#                                               type = 'response'),
#                                   ifelse(Ap == 1 & C == 0,1-predict.glm(weight_model_censor_n1,HERS[HERS$Ap == 1], type = 'response'),
#                                          ifelse(Ap ==0 & C == 1,predict.glm(weight_model_switch_n0,HERS[HERS$Ap == 0],
#                                                                             type = 'response'),
#                                                 1-predict.glm(weight_model_censor_n0,HERS[HERS$Ap == 0],
#                                                               type = 'response'))))/ifelse(Ap == 1 & C == 1, 
#                                   predict.glm(weight_model_censor_d1,HERS[HERS$Ap == 1],
#                                               type = 'response'),
#                                   ifelse(Ap == 1 & C == 0,1-predict.glm(weight_model_censor_d1,HERS[HERS$Ap == 1], type = 'response'),
#                                          ifelse(Ap ==0 & C == 1,predict.glm(weight_model_switch_d0,HERS[HERS$Ap == 0],
#                                                                             type = 'response'),
#                                                 1-predict.glm(weight_model_censor_d0,HERS[HERS$Ap == 0],
#                                                               type = 'response')))),
#                 weight = weight_A * weight_C) %>% 
#   dplyr::group_by(ID) %>% 
#   dplyr::mutate(weight = replace(weight, 1,1)) %>% 
#   dplyr::mutate(weight = cumprod(weight),CA = cumsum(A),
#                 CD4_1 = first(CD4_1),
#                 CD4_2 = first(CD4_2),
#                 viral_1 = first(viral_1),
#                 viral_2 = first(viral_2),
#                 SITE = first(SITE),
#                 WHITE = first(WHITE),
#                 OTHER = first(OTHER),
#                 haartCD4_1 = A*CD4_1) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::mutate(tA = as.factor(t * A),
#                 t = as.factor(t),
#                 followup_time = t,
#                 assigned_treatment = A,
#                 outcome = Y,
#                 CA = cumsum(assigned_treatment)) %>% 
#   dplyr::group_by(ID) %>% 
#   dplyr::filter(
#                 first(CAp) == 0)
# 
# summary(hers_data_tradi$weight)
# sd(hers_data_tradi$weight)
# quantile(hers_data_tradi$weight, c(0.01,0.99))

summary(switch_data$weight)
sd(switch_data$weight)
quantile(switch_data$weight, c(0.01,0.99))


# fit_tradi <- glm(formula = outcome ~  CA +CD4_1 + CD4_2 + 
#                        viral_1 + viral_2 + SITE + WHITE + OTHER + 
#                        CA*CD4_1 , family = binomial(link = "logit"), data = hers_data_tradi, 
#                      weights = hers_data_tradi[["weight"]])
# 
# print(xtable(as.data.frame(summary(fit_tradi)$coefficients), type = "latex", digits = 3))

design_mat <- expand.grid(id = 1:609,
                          trial_period = 0:4,
                          followup_time = 0:4) 
design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]

fitting_data_treatment <-  switch_data %>% 
  dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
  dplyr::select(id,trial_period, followup_time,CD4_1 , CD4_2 , viral_1 , 
                viral_2 , SITE , WHITE , OTHER , assigned_treatment,
                  haartCD4_1) %>% 
  merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>% 
  dplyr::group_by(id) %>% 
  tidyr::fill(CD4_1 , CD4_2 , viral_1 , viral_2 , SITE , WHITE , OTHER , assigned_treatment,
              haartCD4_1,.direction = "down") %>% 
  dplyr::ungroup() %>% 
  dplyr::select(id, trial_period, followup_time,CD4_1 , CD4_2 , viral_1 , 
                viral_2 , SITE , WHITE , OTHER , assigned_treatment,
                haartCD4_1) %>% 
  merge(data.frame(id = switch_data$id, trial_period = switch_data$trial_period), 
        by = c("id", "trial_period"), all.y = TRUE) %>% 
  dplyr::arrange(id, trial_period, followup_time) %>% 
  dplyr::mutate(haartCD4_1 = assigned_treatment*CD4_1) %>% 
  distinct() %>% 
  dplyr::filter(trial_period == 0)

fitting_data_treatment <- fitting_data_treatment[!duplicated(fitting_data_treatment),]  

fitting_data_control <- fitting_data_treatment %>% 
  dplyr::mutate(assigned_treatment = assigned_treatment*0,
                haartCD4_1 = haartCD4_1*0)%>% 
  dplyr::distinct()
####################### ESTIMATES ####################
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
                   risk_difference = survival_control - survival_treatment)

# Y_pred_PP_treatment_tradi <- predict.glm(fit_tradi, 
#                                    fitting_data_treatment, 
#                                    type = "response")
# Y_pred_PP_control_tradi <- predict.glm(fit_tradi, 
#                                  fitting_data_control,
#                                  type = "response")
# predicted_probas_PP_tradi <- fitting_data_treatment %>% 
#   dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment_tradi,
#                 predicted_proba_control = Y_pred_PP_control_tradi) %>% 
#   dplyr::group_by(id, trial_period) %>% 
#   dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
#                 cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::group_by(followup_time) %>% 
#   dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
#                    survival_control = mean(cum_hazard_control),
#                    risk_difference = survival_control - survival_treatment)
# 
# ggplot() +geom_line(aes(x = 0:4,y = pull(predicted_probas_PP_tradi,risk_difference), 
#                         color = 'Traditional IPW-MSM analysis'))+
#   geom_point(aes(x = 0:4,y = pull(predicted_probas_PP_tradi,risk_difference)))+ 
#   geom_line(aes(x = 0:4,y = pull(predicted_probas_PP,risk_difference), 
#                 color = 'Sequential trial emulation analysis'))+
#   geom_point(aes(x = 0:4,y = pull(predicted_probas_PP,risk_difference))) + 
#   scale_color_manual(name = "Per-protocol analysis method",
#                      values = c("Traditional IPW-MSM analysis"= "red",
#                                 "Sequential trial emulation analysis" = "blue")) +
#   labs(x = 'Follow-up time', 
#        y = "Marginal risk difference", title = "HERS data analysis: comparison of\nsequential trial emulation and IPW-MSM")
# 

boot_data_conf <- list()
for (k in 1:500) {
  boot_data_conf[[k]] <- sort(sample(unique(switch_data$id), length(unique(switch_data$id)), replace = TRUE))
}
registerDoParallel(cores = 2)

############################# DIRECT BOOTSTRAP #############################
surv_PP_difference_boostrap_estimates_conf <-as.data.frame(matrix(,5,500))
surv_PP_difference_boostrap_estimates_conf_tradi <-as.data.frame(matrix(,5,500))

for (k in 1:500){
  print(k)
  weights_table_boot <- data.frame(id = 1:609) %>% 
    rowwise() %>% 
    dplyr::mutate(weight_boot = length(boot_data_conf[[k]][boot_data_conf[[k]] == id])) #bootstrap weight is number of times they were sampled
  IP_model <- weight_func_bootstrap(data = HERS, expanded_data = switch_data, 
                                    treatment = 'A',
                                    switch_d_cov = ~CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE + WHITE + OTHER,
                                    cense_d_cov = ~CD4_1 + CD4 + viral + viral_1 + HIVsym + HIVsym_1+ SITE + WHITE + OTHER,
                                    cense = 'C',
                                    weight_model_d0 = switch_d0,
                                    weight_model_n0 = switch_n0,
                                    weight_model_d1 = switch_d1,
                                    weight_model_n1 = switch_n1,
                                    cense_model_d0 = cense_d0,
                                    cense_model_n0 = cense_n0,
                                    cense_model_d1 = cense_d1,
                                    cense_model_n1 = cense_n1, 
                                    boot_idx = boot_data_conf[[k]], remodel = TRUE, quiet = TRUE)
  #calculate IP weights from bootstrap sample
  boot_design_data <- IP_model$data %>%
    merge(weights_table_boot, by = 'id', all.y = TRUE) %>% 
    dplyr::mutate(weight = ifelse(weight_boot !=0,weight*weight_boot,0))
  
  # ############
  # weight_model_switch_d1_boot <- glm(A ~ CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE + WHITE + OTHER,
  #                               family = binomial(link = "logit"), data = HERS[HERS$Ap == 1 & HERS$ID %in% boot_data_conf[[k]]])
  # weight_model_switch_d0_boot <- glm(A ~ CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE + WHITE + OTHER,
  #                               family = binomial(link = "logit"), data = HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]])
  # weight_model_switch_n1_boot <- glm(A ~ 1,
  #                                    family = binomial(link = "logit"), data = HERS[HERS$Ap == 1 & HERS$ID %in% boot_data_conf[[k]]])
  # weight_model_switch_n0_boot <- glm(A ~ 1,
  #                                    family = binomial(link = "logit"), data = HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]])
  # weight_model_censor_d1_boot <- glm(C ~ CD4_1 + CD4 + viral + viral_1 + HIVsym + HIVsym_1+ SITE + WHITE + OTHER,
  #                               family = binomial(link = "logit"), data = HERS[HERS$Ap == 1& HERS$ID %in% boot_data_conf[[k]]])
  # weight_model_censor_d0_boot <- glm(C ~ CD4_1 + CD4 + viral + viral_1 + HIVsym + HIVsym_1+ SITE + WHITE + OTHER,
  #                               family = binomial(link = "logit"), data = HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]])
  # weight_model_censor_n1_boot <- glm(C ~1,
  #                                    family = binomial(link = "logit"), data = HERS[HERS$Ap == 1& HERS$ID %in% boot_data_conf[[k]]])
  # weight_model_censor_n0_boot <- glm(C ~ 1,
  #                                    family = binomial(link = "logit"), data = HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]])
  # boot_design_data_tradi <-  HERS[HERS$ID %in% boot_data_conf[[k]]] %>% 
  #   dplyr::arrange(ID, t) %>% 
  #   dplyr::mutate(weight_A = ifelse(Ap == 1 & A == 1, 
  #                                   predict.glm(weight_model_switch_n1_boot,HERS[HERS$Ap == 1& HERS$ID %in% boot_data_conf[[k]]],
  #                                               type = 'response'),
  #                                   ifelse(Ap == 1 & A == 0,
  #                                          1-predict.glm(weight_model_switch_n1_boot,HERS[HERS$Ap == 1& HERS$ID %in% boot_data_conf[[k]]], type = 'response'),
  #                                          ifelse(Ap ==0 & A == 1,
  #                                                 predict.glm(weight_model_switch_n0_boot,HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]],
  #                                                             type = 'response'),
  #                                                 1-predict.glm(weight_model_switch_n0_boot,HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]],
  #                                                               type = 'response'))))/ifelse(Ap == 1 & A == 1, 
  #                                     predict.glm(weight_model_switch_d1_boot,HERS[HERS$Ap == 1& HERS$ID %in% boot_data_conf[[k]]],
  #                                                 type = 'response'),
  #                                     ifelse(Ap == 1 & A == 0,
  #                                            1-predict.glm(weight_model_switch_d1_boot,HERS[HERS$Ap == 1& HERS$ID %in% boot_data_conf[[k]]], type = 'response'),
  #                                            ifelse(Ap ==0 & A == 1,
  #                                                   predict.glm(weight_model_switch_d0_boot,HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]],
  #                                                               type = 'response'),
  #                                                   1-predict.glm(weight_model_switch_d0_boot,HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]],
  #                                                                 type = 'response')))),
  #                 weight_C = ifelse(Ap == 1 & C == 1, 
  #                                   predict.glm(weight_model_censor_n1_boot,HERS[HERS$Ap == 1& HERS$ID %in% boot_data_conf[[k]]],
  #                                               type = 'response'),
  #                                   ifelse(Ap == 1 & C == 0,1-predict.glm(weight_model_censor_n1_boot,HERS[HERS$Ap == 1& HERS$ID %in% boot_data_conf[[k]]], type = 'response'),
  #                                          ifelse(Ap ==0 & C == 1,predict.glm(weight_model_switch_n0_boot,HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]],
  #                                                                             type = 'response'),
  #                                                 1-predict.glm(weight_model_censor_n0_boot,HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]],
  #                                                               type = 'response'))))/ifelse(Ap == 1 & C == 1, 
  #                                     predict.glm(weight_model_censor_d1_boot,HERS[HERS$Ap == 1& HERS$ID %in% boot_data_conf[[k]]],
  #                                                 type = 'response'),
  #                                     ifelse(Ap == 1 & C == 0,1-predict.glm(weight_model_censor_d1_boot,HERS[HERS$Ap == 1& HERS$ID %in% boot_data_conf[[k]]], type = 'response'),
  #                                            ifelse(Ap ==0 & C == 1,predict.glm(weight_model_switch_d0_boot,HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]],
  #                                                                               type = 'response'),
  #                                                   1-predict.glm(weight_model_censor_d0_boot,HERS[HERS$Ap == 0& HERS$ID %in% boot_data_conf[[k]]],
  #                                                                 type = 'response')))),
  #                 weight = weight_A * weight_C) %>% 
  #   dplyr::group_by(ID) %>% 
  #   dplyr::mutate(weight = replace(weight, 1,1)) %>% 
  #   dplyr::mutate(weight = cumprod(weight),CA = cumsum(A),
  #                 CD4_1 = first(CD4_1),
  #                 CD4_2 = first(CD4_2),
  #                 viral_1 = first(viral_1),
  #                 viral_2 = first(viral_2),
  #                 SITE = first(SITE),
  #                 WHITE = first(WHITE),
  #                 OTHER = first(OTHER),
  #                 haartCD4_1 = A*CD4_1) %>% 
  #   dplyr::ungroup() %>% 
  #   dplyr::mutate(tA = as.factor(t * A),
  #                 t = as.factor(t),
  #                 followup_time = t,
  #                 assigned_treatment = A,
  #                 outcome = Y,
  #                 CA = cumsum(assigned_treatment)) %>% 
  #   dplyr::group_by(ID) %>% 
  #   dplyr::filter(
  #     first(CAp) == 0) %>%
  #   dplyr::mutate(id = ID) %>% 
  #   merge(weights_table_boot, by = 'id', all.y = TRUE) %>% 
  #   dplyr::mutate(weight = ifelse(weight_boot !=0,weight*weight_boot,0))
  # 
  # fit_tradi_boot <- glm(formula = outcome ~ assigned_treatment + CD4_1 + CD4_2 + 
  #                        viral_1 + viral_2 + SITE + WHITE + OTHER + 
  #                        haartCD4_1 + CA , family = binomial(link = "logit"), data = boot_design_data_tradi, 
  #                      weights = boot_design_data_tradi[["weight"]])
  ###########
  PP_boot <- TrialEmulation::trial_msm(data = boot_design_data,
                                            outcome_cov = ~ CD4_1 + CD4_2 + viral_1 + viral_2 + SITE + WHITE + OTHER + assigned_treatment+
                                              haartCD4_1,
                                            model_var = c('assigned_treatment'),
                                            glm_function = 'glm',
                                            include_trial_period = ~1, include_followup_time = ~1,
                                            use_weight=T, use_censor=T, quiet = T, use_sample_weights =  F)
  

  design_mat <- expand.grid(id = 1:tail(boot_design_data$id, n = 1), 
                            trial_period = 0:4,
                            followup_time = 0:4)
  design_mat <- design_mat[which(5 -design_mat$trial_period > design_mat$followup_time),]
  
  # fitting_data_treatment_boot <-  boot_design_data %>% 
  #   dplyr::mutate(assigned_treatment = followup_time*0 + 1) %>% 
  #   dplyr::select(id,trial_period, followup_time,CD4_1 , CD4_2 , viral_1 , 
  #                 viral_2 , SITE , WHITE , OTHER , assigned_treatment,
  #                 haartCD4_1) %>% 
  #   merge(design_mat, by = c("id", "trial_period", "followup_time"), all.y = TRUE) %>% 
  #   dplyr::group_by(id) %>% 
  #   tidyr::fill(CD4_1 , CD4_2 , viral_1 , viral_2 , SITE , WHITE , OTHER , assigned_treatment,
  #               haartCD4_1,.direction = "down") %>% 
  #   dplyr::ungroup() %>% 
  #   dplyr::select(id, trial_period, followup_time,CD4_1 , CD4_2 , viral_1 , 
  #                 viral_2 , SITE , WHITE , OTHER , assigned_treatment,
  #                 haartCD4_1) %>% 
  #   merge(data.frame(id = switch_data$id, trial_period = switch_data$trial_period), 
  #         by = c("id", "trial_period"), all.y = TRUE) %>% 
  #   dplyr::arrange(id, trial_period, followup_time) %>% 
  #   dplyr::mutate(haartCD4_1 = assigned_treatment*CD4_1,
  #                 tA = as.factor(followup_time * assigned_treatment),
  #                 followup_time = as.factor(followup_time)) %>% 
  #   distinct() %>% 
  #   dplyr::filter(trial_period == 0) %>% 
  #   dplyr::group_by(id) %>% 
  #   dplyr::mutate(CA = cumsum(assigned_treatment)) %>% 
  #   dplyr::ungroup()
  # 
  # fitting_data_control_boot <- fitting_data_treatment_boot %>% 
  #   dplyr::mutate(assigned_treatment = assigned_treatment*0,
  #                 tA = as.factor(0),
  #                 CA = cumsum(assigned_treatment)) %>% 
  #   dplyr::distinct()
  # 
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
                     survival_control = mean(cum_hazard_control),
                     risk_difference = survival_control - survival_treatment)
  
  # Y_pred_PP_treatment_tradi_boot <- predict.glm(fit_tradi_boot, 
  #                                              fitting_data_treatment_boot, 
  #                                              type = "response")
  # Y_pred_PP_control_tradi_boot <- predict.glm(fit_tradi_boot, 
  #                                            fitting_data_control_boot,
  #                                            type = "response")
  # predicted_probas_PP_tradi_boot <- fitting_data_treatment_boot %>% 
  #   dplyr::mutate(predicted_proba_treatment = Y_pred_PP_treatment_tradi_boot,
  #                 predicted_proba_control = Y_pred_PP_control_tradi_boot) %>% 
  #   dplyr::group_by(id, trial_period) %>% 
  #   dplyr::mutate(cum_hazard_treatment = cumprod(1-predicted_proba_treatment),
  #                 cum_hazard_control = cumprod(1-predicted_proba_control)) %>% 
  #   dplyr::ungroup() %>% 
  #   dplyr::group_by(followup_time) %>% 
  #   dplyr::summarise(survival_treatment = mean(cum_hazard_treatment),
  #                    survival_control = mean(cum_hazard_control),
  #                    risk_difference = survival_control - survival_treatment)
  
  surv_PP_difference_boostrap_estimates_conf[,k] <-pull(predicted_probas_PP_boot,risk_difference)
  # surv_PP_difference_boostrap_estimates_conf_tradi[,k] <-pull(predicted_probas_PP_tradi_boot,
                                                                  # risk_difference)
}

surv_PP_difference_boostrap_estimates_conf$lb <- apply(surv_PP_difference_boostrap_estimates_conf,
                                                       1,
                                                       quantile,
                                                       probs = c(0.025))
surv_PP_difference_boostrap_estimates_conf$ub <- apply(surv_PP_difference_boostrap_estimates_conf,
                                                       1,
                                                       quantile,
                                                       probs = c(0.975))
CI_bootstrap_coefs_PP_red <- array(, dim = c(5,2))
CI_bootstrap_coefs_PP_red[,1] <-  2*pull(predicted_probas_PP,risk_difference) - surv_PP_difference_boostrap_estimates_conf$ub
CI_bootstrap_coefs_PP_red[,2] <- 2*pull(predicted_probas_PP,risk_difference) - surv_PP_difference_boostrap_estimates_conf$lb

# surv_PP_difference_boostrap_estimates_conf_tradi$lb <- apply(surv_PP_difference_boostrap_estimates_conf_tradi,
#                                                        1,
#                                                        quantile,
#                                                        probs = c(0.025))
# surv_PP_difference_boostrap_estimates_conf_tradi$ub <- apply(surv_PP_difference_boostrap_estimates_conf_tradi,
#                                                        1,
#                                                        quantile,
#                                                        probs = c(0.975))
# CI_bootstrap_coefs_PP_red_tradi <- array(, dim = c(5,2))
# CI_bootstrap_coefs_PP_red_tradi[,1] <-  2*pull(predicted_probas_PP_tradi,risk_difference) - surv_PP_difference_boostrap_estimates_conf_tradi$ub
# CI_bootstrap_coefs_PP_red_tradi[,2] <- 2*pull(predicted_probas_PP_tradi,risk_difference) - surv_PP_difference_boostrap_estimates_conf_tradi$lb
# 
# ggplot() +geom_line(aes(x = 0:4,y = pull(predicted_probas_PP_tradi,risk_difference), 
#                         color = 'Traditional IPW-MSM analysis'))+
#   geom_point(aes(x = 0:4,y = pull(predicted_probas_PP_tradi,risk_difference)))+ 
#   geom_line(aes(x = 0:4,y = pull(predicted_probas_PP,risk_difference), 
#                 color = 'Sequential trial emulation analysis'))+
#   geom_point(aes(x = 0:4,y = pull(predicted_probas_PP,risk_difference))) + 
#   geom_stepribbon(aes(x = 0:4,ymin = CI_bootstrap_coefs_PP_red[,1], 
#                       ymax = CI_bootstrap_coefs_PP_red[,2], color = "Sequential trial emulation analysis"), alpha = 0.1) +
#   geom_stepribbon(aes(x = 0:4,ymin = CI_bootstrap_coefs_PP_red_tradi[,1], 
#                       ymax = CI_bootstrap_coefs_PP_red_tradi[,2], color = "Traditional IPW-MSM analysis"), alpha = 0.1) +
#   scale_color_manual(name = "Per-protocol analysis method",
#                      values = c("Traditional IPW-MSM analysis"= "red",
#                                 "Sequential trial emulation analysis" = "blue")) +
#   labs(x = 'Follow-up time', 
#        y = "Marginal risk difference")


################### LEF OUTCOME ONLY ##################################
X <- model.matrix(PP$model)
e <- PP$model$y - PP$model$fitted.values

surv_PP_difference_LEF_outcome_estimates_conf <- as.data.frame(matrix(,5,500))
for (k in 1:500){
  print(k)
  weights_table_boot <- data.frame(id = 1:609) %>% 
    rowwise() %>% 
    dplyr::mutate(weight_boot = length(boot_data_conf[[k]][boot_data_conf[[k]] == id])) #bootstrap weight is number of times they were sampled
  
  IP_model <- weight_func_bootstrap(data = HERS, expanded_data = switch_data, 
                                    switch_d_cov = ~CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE + WHITE + OTHER,
                                    cense_d_cov = ~CD4_1 + CD4 + viral + viral_1 + HIVsym + HIVsym_1+ SITE + WHITE + OTHER,
                                    cense = 'C',
                                    weight_model_d0 = switch_d0,
                                    weight_model_n0 = switch_n0,
                                    weight_model_d1 = switch_d1,
                                    weight_model_n1 = switch_n1,
                                    cense_model_d0 = cense_d0,
                                    cense_model_n0 = cense_n0,
                                    cense_model_d1 = cense_d1,
                                    cense_model_n1 = cense_n1, 
                                    boot_idx = boot_data_conf[[k]], remodel = TRUE, quiet = TRUE)
  #calculate IP weights from bootstrap sample
  
  boot_design_data <- IP_model$data %>%
    merge(weights_table_boot, by = 'id', all.y = TRUE) %>% 
    dplyr::mutate(weight = ifelse(weight_boot !=0,weight*weight_boot,0)) %>% 
    dplyr::filter(!is.na(trial_period))
  
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
                     survival_control = mean(cum_hazard_control),
                     survival_difference = survival_treatment - survival_control)
  
  surv_PP_difference_LEF_outcome_estimates_conf[,k] <- pull(predicted_probas_PP_boot,survival_difference)
}
surv_PP_difference_LEF_outcome_estimates_conf$lb <-apply(surv_PP_difference_LEF_outcome_estimates_conf,
                                                         1,
                                                         quantile,
                                                         probs = c(0.025))
surv_PP_difference_LEF_outcome_estimates_conf$ub <- apply(surv_PP_difference_LEF_outcome_estimates_conf,
                                                          1,
                                                          quantile,
                                                          probs = c(0.975))

CI_LEF_outcome_coefs_PP_red <- array(,dim = c(5,2))
CI_LEF_outcome_coefs_PP_red[,1] <- 2*pull(predicted_probas_PP,risk_difference) - surv_PP_difference_LEF_outcome_estimates_conf$lb
CI_LEF_outcome_coefs_PP_red[,2] <- 2*pull(predicted_probas_PP,risk_difference) - surv_PP_difference_LEF_outcome_estimates_conf$ub

#################### LEF WEIGHT AND OUTCOME  ##############################
X_sw_d0 <- model.matrix(switch_d0)
e_sw_d0 <- switch_d0$y - switch_d0$fitted.values

X_sw_n0 <- model.matrix(switch_n0)
e_sw_n0 <- switch_n0$y - switch_n0$fitted.values

X_sw_d1 <- model.matrix(switch_d1)
e_sw_d1 <- switch_d1$y - switch_d1$fitted.values

X_sw_n1 <- model.matrix(switch_n1)
e_sw_n1 <- switch_n1$y - switch_n1$fitted.values

X_c_d0 <- model.matrix(cense_d0)
e_c_d0 <- cense_d0$y - cense_d0$fitted.values

X_c_n0 <- model.matrix(cense_n0)
e_c_n0 <- cense_n0$y - cense_n0$fitted.values

X_c_d1 <- model.matrix(cense_d1)
e_c_d1 <- cense_d1$y - cense_d1$fitted.values

X_c_n1 <- model.matrix(cense_n1)
e_c_n1 <- cense_n1$y - cense_n1$fitted.values

surv_PP_difference_LEF_both_estimates_conf <- as.data.frame(matrix(,5,500))
for (k in 1:500){
  print(k)
  weights_table_boot <- data.frame(id = 1:609) %>% 
    rowwise() %>% 
    dplyr::mutate(weight_boot = length(boot_data_conf[[k]][boot_data_conf[[k]] == id])) #bootstrap weight is number of times they were sampled
  
  data_0 <- merge(weights_table_boot, switch_d0$data, on = id, all.y = TRUE)
  data_1 <- merge(weights_table_boot, switch_d1$data, on = id, all.y = TRUE)
  
  LEF_sw_d0_boot <- t(X_sw_d0)%*%(data_0$weight_boot*e_sw_d0)
  LEF_sw_n0_boot <- t(X_sw_n0)%*%(data_0$weight_boot*e_sw_n0)
  LEF_sw_d1_boot <- t(X_sw_d1)%*%(data_1$weight_boot*e_sw_d1)
  LEF_sw_n1_boot <- t(X_sw_n1)%*%(data_1$weight_boot*e_sw_n1)
  LEF_c_d0_boot <- t(X_c_d0)%*%(data_0$weight_boot*e_c_d0)
  LEF_c_n0_boot <- t(X_c_n0)%*%(data_0$weight_boot*e_c_n0)
  LEF_c_d1_boot <- t(X_c_d1)%*%(data_1$weight_boot*e_c_d1)
  LEF_c_n1_boot <- t(X_c_n1)%*%(data_1$weight_boot*e_c_n1)
  
  #Calculate \hat \beta(b)
  beta_sw_d0 <- switch_d0$coefficients + vcov(switch_d0)%*%LEF_sw_d0_boot
  beta_sw_n0 <- switch_n0$coefficients + vcov(switch_n0)%*%LEF_sw_n0_boot
  beta_sw_d1 <- switch_d1$coefficients + vcov(switch_d1)%*%LEF_sw_d1_boot
  beta_sw_n1 <- switch_n1$coefficients + vcov(switch_n1)%*%LEF_sw_n1_boot
  beta_c_d0 <- cense_d0$coefficients + vcov(cense_d0)%*%LEF_c_d0_boot
  beta_c_n0 <- cense_n0$coefficients + vcov(cense_n0)%*%LEF_c_n0_boot
  beta_c_d1 <- cense_d1$coefficients + vcov(cense_d1)%*%LEF_c_d1_boot
  beta_c_n1 <- cense_n1$coefficients + vcov(cense_n1)%*%LEF_c_n1_boot
  
  IP_model <- weight_func_bootstrap(data = HERS, expanded_data = switch_data, 
                                    switch_d_cov = ~CD4_1 + CD4_2 +  viral_1 + viral_2 + HIVsym_1+ SITE + WHITE + OTHER,
                                    cense_d_cov = ~CD4_1 + CD4 + viral + viral_1 + HIVsym + HIVsym_1+ SITE + WHITE + OTHER,
                                    cense = 'C',
                                    weight_model_d0 = switch_d0,
                                    weight_model_n0 = switch_n0,
                                    weight_model_d1 = switch_d1,
                                    weight_model_n1 = switch_n1,
                                    cense_model_d0 = cense_d0,
                                    cense_model_n0 = cense_n0,
                                    cense_model_d1 = cense_d1,
                                    cense_model_n1 = cense_n1, 
                                    new_coef_sw_d0 = beta_sw_d0,
                                    new_coef_sw_n0 = beta_sw_n0,
                                    new_coef_sw_d1 = beta_sw_d1,
                                    new_coef_sw_n1 = beta_sw_n1,
                                    new_coef_c_d0 = beta_c_d0,
                                    new_coef_c_n0 = beta_c_n0,
                                    new_coef_c_d1 = beta_c_d1,
                                    new_coef_c_n1 = beta_c_n1,
                                    boot_idx = boot_data_conf[[k]], remodel = FALSE, quiet = TRUE)
  
  #calculate IP weights from bootstrap sample
  
  boot_design_data <- IP_model$data %>%
    merge(weights_table_boot, by = 'id', all.y = TRUE) %>% 
    dplyr::mutate(weight = ifelse(weight_boot !=0,weight*weight_boot,0))%>% 
    dplyr::filter(!is.na(trial_period))
  
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
                     survival_control = mean(cum_hazard_control),
                     survival_difference = survival_treatment - survival_control)
  
  surv_PP_difference_LEF_both_estimates_conf[,k] <- pull(predicted_probas_PP_boot, survival_difference)
}
surv_PP_difference_LEF_both_estimates_conf$lb <-apply(surv_PP_difference_LEF_both_estimates_conf,
                                                      1,
                                                      quantile,
                                                      probs = c(0.025))
surv_PP_difference_LEF_both_estimates_conf$ub <- apply(surv_PP_difference_LEF_both_estimates_conf,
                                                       1,
                                                       quantile,
                                                       probs = c(0.975))


CI_LEF_both_coefs_PP_red <- array(,dim = c(5,2))
CI_LEF_both_coefs_PP_red[,1] <- 2*pull(predicted_probas_PP, risk_difference) - surv_PP_difference_LEF_both_estimates_conf$lb
CI_LEF_both_coefs_PP_red[,2] <- 2*pull(predicted_probas_PP, risk_difference) - surv_PP_difference_LEF_both_estimates_conf$ub

############################SANDWICH #######################################
covariance_mat <-PP$robust$matrix
if (all(eigen(covariance_mat)$values > 0) == F){
  not_pos_def <- not_pos_def + 1.0
  next
}
#Step 1 of algorithm  -- sampling Y_n1, ..., Y_nB ~ MN(coeffs,sandwich covariance)
sampling_size <- 500
coeffs_sample <- mvrnorm(sampling_size,PP$model$coefficients, covariance_mat)
surv_PP_difference_sandwich_estimates <- as.data.frame(matrix(,5,500))
for (k in 1:500){
  print(k)
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
                     survival_control = mean(cum_hazard_control),
                     risk_difference = survival_control - survival_treatment)
  
  surv_PP_difference_sandwich_estimates[,k] <- pull(predicted_probas_PP_sample, risk_difference)
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

CI_sandwich_coefs_PP_red <- array(, dim = c(5,2))
CI_sandwich_coefs_PP_red[,1] <- surv_PP_difference_sandwich_estimates$lb
CI_sandwich_coefs_PP_red[,2] <- surv_PP_difference_sandwich_estimates$ub
#
#################### CURVE PLOT AND CIS ##################################

ggplot() +geom_line(aes(x = 0:4,y = pull(predicted_probas_PP_tradi,risk_difference), 
                        color = 'Traditional IPW-MSM analysis'))+
  geom_point(aes(x = 0:4,y = pull(predicted_probas_PP_tradi,risk_difference)))+ 
  geom_line(aes(x = 0:4,y = pull(predicted_probas_PP,risk_difference), 
                color = 'Sequential trial emulation analysis'))+
  geom_point(aes(x = 0:4,y = pull(predicted_probas_PP,risk_difference))) + 
  geom_stepribbon(aes(x = 0:4,ymin = CI_bootstrap_coefs_PP_red[,1], 
                      ymax = CI_bootstrap_coefs_PP_red[,2], color = "Sequential trial emulation analysis"), alpha = 0.1) +
  geom_stepribbon(aes(x = 0:4,ymin = CI_bootstrap_coefs_PP_red_tradi[,1], 
                      ymax = CI_bootstrap_coefs_PP_red_tradi[,2], color = "Traditional IPW-MSM analysis"), alpha = 0.1) +
  scale_color_manual(name = "Per-protocol analysis method",
                     values = c("Traditional IPW-MSM analysis"= "red",
                                "Sequential trial emulation analysis" = "blue")) +
  labs(x = 'Follow-up time', 
       y = "Marginal risk difference", title = "HERS data analysis")

ggplot(data = predicted_probas_PP,aes(x = 0:4)) +
  geom_step(aes(x = 0:4,y = risk_difference)) +
  geom_stepribbon(aes(ymin = CI_LEF_outcome_coefs_PP_red[,1],
                      ymax = CI_LEF_outcome_coefs_PP_red[,2], color = "LEF outcome"), alpha = 0.1) +
  geom_stepribbon(aes(ymin = CI_LEF_both_coefs_PP_red[,1],
                      ymax = CI_LEF_both_coefs_PP_red[,2], color = "LEF both"), alpha = 0.1) +
  geom_stepribbon(aes(ymin = CI_sandwich_coefs_PP_red[,1],
                      ymax = CI_sandwich_coefs_PP_red[,2], color = "Sandwich"), alpha = 0.1) +
  geom_stepribbon(aes(ymin = CI_bootstrap_coefs_PP_red[,1], 
                      ymax = CI_bootstrap_coefs_PP_red[,2], color = "Nonparametric Boot."), alpha = 0.1) +
  scale_color_manual(name = "CI method", values = c("LEF outcome"= "green", "Nonparametric Boot." = "red", "LEF both" = 'purple', "Sandwich" = 'blue')) +
  labs(x = 'Follow-up time', 
       y = "Marginal risk difference")

 