library(modelr)
library(tidyverse)
library(tidyr)
source("simulate_MSM.R")
set.seed(20222022)
library(RandomisedTrialsEmulation)
library(MASS)
library(sandwich)
library(foreach)
library(doParallel)
library(parallel)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(pammtools)

simdata_censored_conf_treat<-DATA_GEN_censored(1000000, 10, conf = j, all_treat = T, censor = F)
simdata_censored_conf_control<-DATA_GEN_censored(1000000, 10, conf = j, all_control = T, censor = F)

surv_data_treat <- simdata_censored_conf_treat[ !duplicated(simdata_censored_conf_treat[, c("ID")], fromLast=T),] %>% 
  dplyr::mutate(status = Y) %>% 
  dplyr::select(ID, t, status)

f1 <- survfit(Surv(t, status) ~ 1, data = surv_data_treat)

plot(f1, 
     xlab = "Days", 
     ylab = "Overall survival probability ")


surv_data_control <- simdata_censored_conf_control[ !duplicated(simdata_censored_conf_control[, c("ID")], fromLast=T),] %>% 
  dplyr::mutate(status = Y) %>% 
  dplyr::select(ID, t, status)

f2 <- survfit(Surv(t, status) ~ 1, data = surv_data_control)

plot(f2, 
     xlab = "Days", 
     ylab = "Overall survival probability")


ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = f1$surv, color = "Treatment")) +
  geom_step(aes(y = f2$surv, color = "Control")) + 
  scale_color_manual(name = "Treatment assignment", values = c("Treatment"= "red", "Control" = "blue")) +
  labs(x = 'Follow-up time', y = "Survival function", color = "Assigned treatment", title = "PP Analysis")
