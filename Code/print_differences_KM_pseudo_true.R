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

load('true_value_red_newsimus_pseudo_true_1-5M.rda')
true_value_MSM_1_5M <- true_value_red

load('true_value_red_newsimus.rda')
true_value_red <- -true_value_red

differences_1_5M_KM <- rbind(cbind(0:4,round(true_value_MSM_1_5M[,,1] - true_value_red[,,1], digits = 3)),
                             cbind(0:4,round(true_value_MSM_1_5M[,,2] - true_value_red[,,2], digits = 3)),
                                   cbind(0:4,round(true_value_MSM_1_5M[,,3] - true_value_red[,,3], digits = 3)))

print(xtable(differences_1_5M_KM, type = 'latex', digits = 3), include.rownames = F)

load('true_value_red_pseudo_true_boot_200it.rda')
true_value_boot_200 <- true_value_boot

differences_200boot_KM <- rbind(cbind(0:4,round(true_value_boot_200[,,1] - true_value_red[,,1], digits = 4)),
                                cbind(0:4,round(true_value_boot_200[,,2] - true_value_red[,,2], digits = 4)),
                                cbind(0:4,round(true_value_boot_200[,,3] - true_value_red[,,3], digits = 4)))
print(xtable(differences_200boot_KM, type = 'latex', digits = c(0,0,4,4,4,4,4,4,4,4,4)), include.rownames = F)

load('true_value_red_pseudo_true_boot_500it.rda')
true_value_boot_500 <- true_value_boot

differences_500boot_KM <- rbind(cbind(0:4,round(true_value_boot_500[,,1] - true_value_red[,,1], digits = 4)),
                                cbind(0:4,round(true_value_boot_500[,,2] - true_value_red[,,2], digits = 4)),
                                cbind(0:4,round(true_value_boot_500[,,3] - true_value_red[,,3], digits = 4)))
print(xtable(differences_500boot_KM, type = 'latex', digits = c(0,0,4,4,4,4,4,4,4,4,4)), include.rownames = F)

load("~/Multiple-trial-emulation-IPTW-MSM-CIs/Code/true_value_red_pseudo_true_boot_200it_200000p.rda")

differences_200boot_200k_KM <- rbind(cbind(0:4,round(true_value_boot[,,1] - true_value_red[,,1], digits = 4)),
                                cbind(0:4,round(true_value_boot[,,2] - true_value_red[,,2], digits = 4)),
                                cbind(0:4,round(true_value_boot[,,3] - true_value_red[,,3], digits = 4)))
print(xtable(differences_200boot_200k_KM, type = 'latex', digits = c(0,0,4,4,4,4,4,4,4,4,4)), include.rownames = F)
