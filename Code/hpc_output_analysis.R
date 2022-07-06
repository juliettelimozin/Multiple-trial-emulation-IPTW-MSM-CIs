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
library(doRNG)

load("true_value_conf.rda")
load("true_value_treat.rda")
load("HPC output/est_true_value_conf.rda")
load("HPC output/est_true_value_treat.rda")

bootstrap_coefs <- array(,dim = c(10,2,1000,9))
sandwich_coefs <- array(,dim = c(10,2,1000,9))
bootstrap_treat <- array(,dim = c(10,2,1000,9))
sandwich_treat <- array(,dim = c(10,2,1000,9))

for (i in 1:9){
  load(paste0("HPC output/CI_bootstrap_coefs_PP_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_coefs_PP_", i, ".rda"))
  load(paste0("HPC output/CI_bootstrap_treat_PP_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_treat_PP_", i, ".rda"))
  bootstrap_coefs[,,,i] <- CI_bootstrap_coefs_PP
  sandwich_coefs[,,,i] <- CI_sandwich_coefs_PP
  bootstrap_treat[,,,i] <- CI_bootstrap_treat_PP
  sandwich_treat[,,,i] <- CI_sandwich_treat_PP
}

mean_lengths <- matrix(0, nrow = 4, ncol = 9)

for (i in 1:9){
  mean_lengths[1,i] <- mean(rowMeans(bootstrap_coefs[,2,,i]- bootstrap_coefs[,1,,i], na.rm = TRUE))
  mean_lengths[2,i] <- mean(rowMeans(sandwich_coefs[,2,,i]- sandwich_coefs[,1,,i], na.rm = TRUE))
  mean_lengths[3,i] <- mean(rowMeans(bootstrap_treat[,2,,i]- bootstrap_treat[,1,,i], na.rm = TRUE))
  mean_lengths[4,i] <- mean(rowMeans(sandwich_treat[,2,,i]- sandwich_treat[,1,,i], na.rm = TRUE))
}

coverage_ind <- matrix(0,nrow = 4, ncol = 9)
coverage_ind_est <- matrix(0,nrow = 4, ncol = 9)
success <- matrix(0,nrow = 4, ncol = 9)
for (i in 1:1000){
  for (j in 1:9){
    
    if (all(is.na(bootstrap_coefs[10,1,i,j])) == F){
      success[1,j] <- success[1,j] + 1
      if (all(bootstrap_coefs[10,1,i,j] < true_value_conf[10,1,j] - true_value_conf[10,2,j]) 
          & all(bootstrap_coefs[10,2,i,j] > true_value_conf[10,1,j] - true_value_conf[10,2,j])){
        coverage_ind[1,j] <- coverage_ind[1,j] + 1
      }
      if (all(bootstrap_coefs[10,1,i,j] < est_true_value_conf[10,1,j] - est_true_value_conf[10,2,j]) 
          & all(bootstrap_coefs[10,2,i,j] > est_true_value_conf[10,1,j] - est_true_value_conf[10,2,j])){
        coverage_ind_est[1,j] <- coverage_ind_est[1,j] + 1
      }
    }
    
    if (all(is.na(sandwich_coefs[10,1,i,j])) == F){
      success[2,j] <- success[2,j] + 1
      if (all(sandwich_coefs[10,1,i,j] < true_value_conf[10,1,j] - true_value_conf[10,2,j]) 
        & all(sandwich_coefs[10,2,i,j] > true_value_conf[10,1,j] - true_value_conf[10,2,j])){
        coverage_ind[2,j] <- coverage_ind[2,j] + 1
      }
      if (all(sandwich_coefs[10,1,i,j] < est_true_value_conf[10,1,j] - est_true_value_conf[10,2,j]) 
          & all(sandwich_coefs[10,2,i,j] > est_true_value_conf[10,1,j] - est_true_value_conf[10,2,j])){
        coverage_ind_est[2,j] <- coverage_ind_est[2,j] + 1
      }
    }
    if (all(is.na(bootstrap_treat[10,1,i,j])) == F){
      success[3,j] <- success[3,j] + 1
      if (all(bootstrap_treat[10,1,i,j] < true_value_treat[10,1,j] - true_value_treat[10,2,j]) 
          & all(bootstrap_treat[10,2,i,j] > true_value_treat[10,1,j] - true_value_treat[10,2,j])){
        coverage_ind[3,j] <- coverage_ind[3,j] + 1
      }
      if (all(bootstrap_treat[10,1,i,j] < est_true_value_treat[10,1,j] - est_true_value_treat[10,2,j]) 
          & all(bootstrap_treat[10,2,i,j] > est_true_value_treat[10,1,j] - est_true_value_treat[10,2,j])){
        coverage_ind_est[3,j] <- coverage_ind_est[3,j] + 1
      }
    }
    if (all(is.na(sandwich_treat[10,1,i,j])) == F){
      success[4,j] <- success[4,j] + 1
      if (all(sandwich_treat[10,1,i,j] < true_value_treat[10,1,j] - true_value_treat[10,2,j]) 
          & all(sandwich_treat[10,2,i,j] > true_value_treat[10,1,j] - true_value_treat[10,2,j])){
        coverage_ind[4,j] <- coverage_ind[4,j] + 1
      }
      if (all(sandwich_treat[10,1,i,j] < est_true_value_treat[10,1,j] - est_true_value_treat[10,2,j]) 
          & all(sandwich_treat[10,2,i,j] > est_true_value_treat[10,1,j] - est_true_value_treat[10,2,j])){
        coverage_ind_est[4,j] <- coverage_ind_est[4,j] + 1
      }
    }
  }
}

coverage_ind <- coverage_ind/success
coverage_ind_est <- coverage_ind_est/success

ggplot() +
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[1,], colour = "Bootstrap",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[2,], colour = "Sandwich",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[2,], colour = "Sandwich")) +
  
  geom_line(aes(x = 1:9/10, y = coverage_ind[1,], colour = "Bootstrap",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind[2,], colour = "Sandwich",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[2,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  scale_linetype_manual(name = "Comparison type", values = c("Kaplan-Meier"= "solid", "Extreme sample" = "dashed")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Confounding strength', 
       y = "Empirical coverage rate", color = "CI type", 
       linetype = "Comparison type",
       title = "Outcome model with interaction terms")

ggplot() +
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[3,], colour = "Bootstrap",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[3,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[4,], colour = "Sandwich",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[4,], colour = "Sandwich")) +
  
  geom_line(aes(x = 1:9/10, y = coverage_ind[3,], colour = "Bootstrap",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[3,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind[4,], colour = "Sandwich",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[4,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  scale_linetype_manual(name = "Comparison type", values = c("Kaplan-Meier"= "solid", "Extreme sample" = "dashed")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Prevalence of treatment', 
       y = "Empirical coverage rate", color = "CI type", 
       linetype = "Comparison type",
       title = "Outcome model with interaction terms")
