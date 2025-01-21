#!/usr/bin R
.libPaths()
library(modelr)
library(tidyverse)
library(tidyr)
setwd("~/rds/hpc-work/Multiple-trial-emulation-IPTW-MSM-CIs/Code")
source("simulate_MSM_simplified.R")
source("weight_func.R")
set.seed(NULL)
library(TrialEmulation, lib.loc = '/home/jml219/R/x86_64-redhat-linux-gnu-library/4.3')
library(MASS)
library(sandwich)
library(doParallel)
library(doRNG)
library(rlist)

size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)
scenarios <- tidyr::crossing(size,conf, treat)

assigned_treatment <- 0:1
IPW_summary_treated <- scenarios %>% 
  dplyr::mutate(assigned_treatment = 1,Minimum = 0, Q1 = 0, Mean = 0, Median = 0, Q3 = 0, Maximum = 0)
IPW_summary_control <- scenarios %>% 
  dplyr::mutate(assigned_treatment = 0,Minimum = 0, Q1 = 0, Mean = 0, Median = 0, Q3 = 0, Maximum = 0)

for(l in 1:27){
simdata_censored<-DATA_GEN_censored_reduced(as.numeric(scenarios[l,1]), 5, 
                                            conf = as.numeric(scenarios[l,2]), 
                                            treat_prev = as.numeric(scenarios[l,3]),
                                            outcome_prev = -4.7,
                                            censor = F)
PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                            eligible ='eligible',
                                            switch_d_cov = ~X2 + X4,
                                            outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                            use_weight=T, use_censor=T, quiet = T,
                                            save_weight_models = F,
                                            data_dir = getwd())
IPW_summary_treated[l,5:10] <- t(rbind(min(PP_prep$data[PP_prep$data$assigned_treatment == 1,]$weight), 
                            quantile(PP_prep$data[PP_prep$data$assigned_treatment == 1,]$weight, probs = 0.25), 
                            mean(PP_prep$data[PP_prep$data$assigned_treatment == 1,]$weight), 
                            median(PP_prep$data[PP_prep$data$assigned_treatment == 1,]$weight),
                            quantile(PP_prep$data[PP_prep$data$assigned_treatment == 1,]$weight, probs = 0.75),
                            max(PP_prep$data[PP_prep$data$assigned_treatment == 1,]$weight)))
IPW_summary_control[l,5:10] <- t(rbind(min(PP_prep$data[PP_prep$data$assigned_treatment == 0,]$weight), 
                                      quantile(PP_prep$data[PP_prep$data$assigned_treatment ==0,]$weight, probs = 0.25), 
                                      mean(PP_prep$data[PP_prep$data$assigned_treatment == 0,]$weight), 
                                      median(PP_prep$data[PP_prep$data$assigned_treatment == 0,]$weight),
                                      quantile(PP_prep$data[PP_prep$data$assigned_treatment == 0,]$weight, probs = 0.75),
                                      max(PP_prep$data[PP_prep$data$assigned_treatment == 0,]$weight)))
}

IPW_summary <- rbind(IPW_summary_control, IPW_summary_treated) %>% 
  dplyr::arrange(size, conf, treat, assigned_treatment)

xftbl <- xtableFtable(IPW_summary, method = "compact")
print(xtable(IPW_summary, type = 'latex', digits = c(0,0,1,0,0,2,2,2,2,2,2)), include.rownames = F) 
