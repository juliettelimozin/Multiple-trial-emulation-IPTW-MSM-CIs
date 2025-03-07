library(dplyr)
library(tidyverse)
library(tidyr)
setwd("/home/jml219/rds/hpc-work/Project1")
library(ggplot2)
library(ggpubr)
#load("Rdata.RData")
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
iters <- 1000
bootstrap_iter <- 500
sampling_size <- 500

treat_pos <- c(-1,-0.8,-0.5,-0.2,0,0.2,0.5,0.8,1)
outcomes <- c("low", 'med', 'high')

sandwich_SE <- array(, dim = c(5,iters,27,3))
bootstrap_SE <-array(, dim = c(5, iters,27,3))
LEF_outcome_SE <- array(, dim = c(5, iters,27,3))
LEF_both_SE <- array(, dim = c(5, iters,27,3))
jackknife_mvn_SE <- array(, dim = c(5, iters,27,3))

size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)

scenarios <- as.data.frame(tidyr::crossing(size,conf, treat))

for (i in 1:9){
  for (j in 1:3){
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ_fixed/J_bootstrap_mrd_',outcomes[j],'_',i, '.rda'))
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ_fixed/J_LEF_outcome_mrd_',outcomes[j],'_',i, '.rda'))
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ_fixed/J_LEF_both_mrd_',outcomes[j],'_',i, '.rda'))
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ_fixed/J_sandwich_mrd_',outcomes[j],'_',i, '.rda'))
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ_fixed/J_jackknife_mvn_mrd_',outcomes[j],'_',i, '.rda'))

    for (k in 1:5){
      sandwich_SE[k,,i,j] <-colSds(sandwich_mrd[k,,], na.rm = TRUE)
      bootstrap_SE[k,,i,j] <- colSds(bootstrap_mrd[k,,], na.rm = TRUE)
      LEF_outcome_SE[k,,i,j] <- colSds(LEF_outcome_mrd[k,,], na.rm = TRUE)
      LEF_both_SE[k,,i,j] <- colSds(LEF_both_mrd[k,,], na.rm = TRUE)
      if (i %in% 1:9){
        jackknife_mvn_SE[k,,i,j] <- colSds(jackknife_mvn_mrd[k,,], na.rm = TRUE)}}
  }}
save(sandwich_SE, file = '~/rds/hpc-work/Project1/NewSimusJ_fixed/J_sandwich_SE.rda')
save(bootstrap_SE, file = '~/rds/hpc-work/Project1/NewSimusJ_fixed/J_bootstrap_SE.rda')
save(LEF_outcome_SE, file = '~/rds/hpc-work/Project1/NewSimusJ_fixed/J_LEF_outcome_SE.rda')
save(LEF_both_SE, file = '~/rds/hpc-work/Project1/NewSimusJ_fixed/J_LEF_both_SE.rda')
save(jackknife_mvn_SE, file = '~/rds/hpc-work/Project1/NewSimusJ_fixed/J_jackknife_mvn_SE.rda')