#!/usr/bin R
library(modelr)
library(tidyverse)
library(tidyr)
setwd("Code")
source("simulate_MSM_simplified.R")
source("weight_func.R")
set.seed(NULL)
library(TrialEmulation)
library(MASS)
library(sandwich)
library(doParallel)
library(doRNG)

simdata_censored<-DATA_GEN_censored_reduced(200, 5, 
                                            conf = 0.1, 
                                            treat_prev = -1,
                                            outcome_prev = -4.7,
                                            censor = F)

event_summary <- simdata_censored %>% 
  dplyr::group_by(t, A) %>% 
  dplyr::summarise(sum_y = sum(Y))

PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                            eligible ='eligible',
                                            switch_d_cov = ~X2 + X4,
                                            outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                            use_weight=T, use_censor=T, quiet = T,
                                            save_weight_models = T,
                                            data_dir = getwd())
event_summary <- PP_prep$data %>% 
  dplyr::group_by(followup_time, assigned_treatment) %>% 
  dplyr::summarise(sum_y = sum(outcome))
  
  