library(dplyr)
library(tidyverse)
library(tidyr)
setwd("~/Multiple-trial-emulation-IPTW-MSM-CIs/Code")
source("simulate_MSM_simplified.R")
library(TrialEmulation)
library(ggplot2)
library(ggpubr)
library(modelr)
library(MASS)
library(sandwich)
library(foreach)
library(doParallel)
library(parallel)
library(survival)
library(survminer)
library(lubridate)
library(pammtools)
library(doRNG)
library(matrixStats)
library(latex2exp)
library(grDevices)
library(xtable)
library(grDevices)
treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)
outcome_prev <- c(-4.7,-3.8,-3)

scenarios <- tidyr::crossing(conf, treat)

true_value_red <- array(,dim = c(5,9,3))
surv0 <- array(,dim = c(5,9,3))
surv1 <- array(,dim = c(5,9,3))

l = 9
j = 3

simdata_censored <-DATA_GEN_censored_reduced(1000000, 5, 
                                             conf = as.numeric(scenarios[l,1]), 
                                             treat_prev = as.numeric(scenarios[l,2]),
                                             outcome_prev = outcome_prev[j],
                                             censor = F)
PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                            eligible ='eligible',
                                            estimand_type = 'PP',
                                            switch_d_cov = ~X2 + X4,
                                            outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                            quiet = T,
                                            save_weight_models = F)

con4<-xtabs(~assigned_treatment + outcome + followup_time, data=PP_prep$data)
ftable(con4)

xftbl <- xtableFtable(ftable(con4), method = "compact")
print.xtableFtable(xftbl, booktabs = T) 