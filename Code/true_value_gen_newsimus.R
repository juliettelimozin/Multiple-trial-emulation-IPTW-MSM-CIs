library(modelr)
library(tidyverse)
library(tidyr)
setwd('~/rds/hpc-work/Project1')
source("simulate_MSM_simplified.R")

set.seed(20222022)
library(MASS)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(pammtools)

treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)
outcome_prev <- c(-4.7,-3.8,-3)

scenarios <- tidyr::crossing(conf, treat)

true_value_red <- array(,dim = c(5,9,3))
surv0 <- array(,dim = c(5,9,3))
surv1 <- array(,dim = c(5,9,3))

for (l in 1:9){
  for (j in 1:3){
  simdata_censored_treat<-DATA_GEN_censored_reduced(1000000, 5, 
                                                         conf = as.numeric(scenarios[l,1]), 
                                                         treat_prev = as.numeric(scenarios[l,2]),
                                                         outcome_prev = outcome_prev[j],
                                                         all_treat = T,
                                                         censor = F)
  simdata_censored_control<-DATA_GEN_censored_reduced(1000000,5, 
                                                           conf = as.numeric(scenarios[l,1]), 
                                                           treat_prev = as.numeric(scenarios[l,2]),
                                                           outcome_prev = outcome_prev[j],
                                                           all_control = T,
                                                           censor = F)
  
  surv_data_treat <- simdata_censored_treat[ !duplicated(simdata_censored_treat[, c("ID")], fromLast=T),] %>% 
    dplyr::mutate(status = Y) %>% 
    dplyr::select(ID, t, status)
  
  f1 <- survfit(Surv(t, status) ~ 1, data = surv_data_treat)
  
  surv_data_control <- simdata_censored_control[ !duplicated(simdata_censored_control[, c("ID")], fromLast=T),] %>% 
    dplyr::mutate(status = Y) %>% 
    dplyr::select(ID, t, status)
  
  f2 <- survfit(Surv(t, status) ~ 1, data = surv_data_control)
  
  
  true_value_red[,l,j] <- f1$surv - f2$surv
  
  surv0[,l,j]<- f2$surv
  surv1[,l,j] <- f1$surv
  }
}
save(true_value_red, file = "true_value_red_newsimus.rda")
save(surv0, file = "true_value_surv0.rda")
save(surv1, file = "true_value_surv1.rda")

