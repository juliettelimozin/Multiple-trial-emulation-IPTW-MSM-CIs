library(modelr)
library(tidyverse)
library(tidyr)
source("simulate_MSM.R")
set.seed(20222022)
library(MASS)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(pammtools)

conf <- 1:9/10
true_value_conf <- array(,dim = c(10,2,9))
true_value_treat <- array(,dim = c(10,2,9))

for (i in conf){
  simdata_censored_conf_treat<-DATA_GEN_censored(1000000, 10, conf = i, all_treat = T, censor = F)
  simdata_censored_conf_control<-DATA_GEN_censored(1000000, 10, conf = i, all_control = T, censor = F)
  
  surv_data_treat <- simdata_censored_conf_treat[ !duplicated(simdata_censored_conf_treat[, c("ID")], fromLast=T),] %>% 
    dplyr::mutate(status = Y) %>% 
    dplyr::select(ID, t, status)
  
  f1 <- survfit(Surv(t, status) ~ 1, data = surv_data_treat)

  surv_data_control <- simdata_censored_conf_control[ !duplicated(simdata_censored_conf_control[, c("ID")], fromLast=T),] %>% 
    dplyr::mutate(status = Y) %>% 
    dplyr::select(ID, t, status)
  
  f2 <- survfit(Surv(t, status) ~ 1, data = surv_data_control)
  
  true_value_conf[,1,10*i] <- f1$surv
  true_value_conf[,2,10*i] <- f2$surv
  
  simdata_censored_treat_treat<-DATA_GEN_censored(1000000, 10, treat_prev = i, all_treat = T, censor = F)
  simdata_censored_treat_control<-DATA_GEN_censored(1000000, 10, treat_prev = i, all_control = T, censor = F)
  
  surv_data_treat <- simdata_censored_treat_treat[ !duplicated(simdata_censored_treat_treat[, c("ID")], fromLast=T),] %>% 
    dplyr::mutate(status = Y) %>% 
    dplyr::select(ID, t, status)
  
  f1 <- survfit(Surv(t, status) ~ 1, data = surv_data_treat)
  
  surv_data_control <- simdata_censored_treat_control[ !duplicated(simdata_censored_treat_control[, c("ID")], fromLast=T),] %>% 
    dplyr::mutate(status = Y) %>% 
    dplyr::select(ID, t, status)
  
  f2 <- survfit(Surv(t, status) ~ 1, data = surv_data_control)
  
  true_value_treat[,1,10*i] <- f1$surv
  true_value_treat[,2,10*i] <- f2$surv
}

save(true_value_conf, file = "true_value_conf.rda")
save(true_value_treat, file = "true_value_treat.rda")
