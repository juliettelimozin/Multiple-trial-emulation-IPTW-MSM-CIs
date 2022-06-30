library(modelr)
library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/MPhil-dissertation/Code")
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

load("true_value_conf.rda")
load("true_value_treat.rda")
load("HPC output/CI_bootstrap_coefs_PP_2.rda")
load("HPC output/CI_sandwich_coefs_PP_2.rda")
load("HPC output/CI_bootstrap_treat_PP_1.rda")
load("HPC output/CI_sandwich_treat_PP_1.rda")

for (iter in 1:20){
  print(ggplot(,aes(x = 1:10)) +
    geom_step(aes(y = true_value_treat[,1,1] - true_value_treat[,2,1])) +
    geom_stepribbon(aes(ymin = CI_bootstrap_treat_PP[,1,iter], ymax = CI_bootstrap_treat_PP[,2,iter]), alpha = 0.3) +
    geom_stepribbon(aes(ymin = CI_sandwich_treat_PP[,1,iter], ymax = CI_sandwich_treat_PP[,2,iter]), alpha = 0.1) +
    labs(x = 'Follow-up time', y = "Survival function", color = "Assigned treatment", title = paste0("PP Analysis ",iter)))
}

ggplot(data = predicted_probas_PP, aes(followup_time)) +
  geom_step(aes(y = survival_difference)) +
  geom_stepribbon(aes(ymin = survival_difference_lb_sandwich, 
                      ymax = survival_difference_ub_sandwich), alpha = 0.1, colour = 'red') +
  geom_stepribbon(aes(ymin = CI_bootstrap_coefs_PP[,1,1], ymax = CI_bootstrap_coefs_PP[,2,1]), alpha = 0.3) +
  geom_stepribbon(aes(ymin = CI_sandwich_coefs_PP[,1,1], ymax = CI_sandwich_coefs_PP[,2,1]), alpha = 0.1, colour = 'blue') +
  labs(x = 'Follow-up time', 
       y = "Difference between treatment and control survival functions", title = "PP analysis- sandwich CIs")
