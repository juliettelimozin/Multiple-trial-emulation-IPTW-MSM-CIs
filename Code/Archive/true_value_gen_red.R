library(modelr)
library(tidyverse)
library(tidyr)
source("simulate_MSM_simplified.R")
setwd('~/rds/hpc-work')
set.seed(20222022)
library(MASS)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(pammtools)

conf <- 1:9/10
treat <- c(-1,-0.8,-0.5,-0.2,0,0.2,0.5,0.8,1)

true_value_conf_red <- array(,dim = c(5,2,9))
true_value_treat_red <- array(,dim = c(5,2,9))

for (i in conf){
  simdata_censored_conf_treat<-DATA_GEN_censored_reduced(1000000, 5, conf = i, all_treat = T, censor = F)
  simdata_censored_conf_control<-DATA_GEN_censored_reduced(1000000, 5, conf = i, all_control = T, censor = F)
  
  surv_data_treat <- simdata_censored_conf_treat[ !duplicated(simdata_censored_conf_treat[, c("ID")], fromLast=T),] %>% 
    dplyr::mutate(status = Y) %>% 
    dplyr::select(ID, t, status)
  
  f1 <- survfit(Surv(t, status) ~ 1, data = surv_data_treat)
  
  surv_data_control <- simdata_censored_conf_control[ !duplicated(simdata_censored_conf_control[, c("ID")], fromLast=T),] %>% 
    dplyr::mutate(status = Y) %>% 
    dplyr::select(ID, t, status)
  
  f2 <- survfit(Surv(t, status) ~ 1, data = surv_data_control)
  
  true_value_conf_red[,1,10*i] <- f1$surv
  true_value_conf_red[,2,10*i] <- f2$surv
}
save(true_value_conf_red, file = "true_value_conf_red_low.rda")

for (i in 1:9){
  simdata_censored_treat_treat<-DATA_GEN_censored_reduced(1000000, 5, treat_prev = treat[i], all_treat = T, censor = F)
  simdata_censored_treat_control<-DATA_GEN_censored_reduced(1000000, 5, treat_prev = treat[i], all_control = T, censor = F)
  
  surv_data_treat <- simdata_censored_treat_treat[ !duplicated(simdata_censored_treat_treat[, c("ID")], fromLast=T),] %>% 
    dplyr::mutate(status = Y) %>% 
    dplyr::select(ID, t, status)
  
  f1 <- survfit(Surv(t, status) ~ 1, data = surv_data_treat)
  
  surv_data_control <- simdata_censored_treat_control[ !duplicated(simdata_censored_treat_control[, c("ID")], fromLast=T),] %>% 
    dplyr::mutate(status = Y) %>% 
    dplyr::select(ID, t, status)
  
  f2 <- survfit(Surv(t, status) ~ 1, data = surv_data_control)
  
  true_value_treat_red[,1,i] <- f1$surv
  true_value_treat_red[,2,i] <- f2$surv
  
}

save(true_value_treat_red, file = "true_value_treat_red_low.rda")