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

est_true_value_conf_full <- est_true_value_conf
est_true_value_treat_full <- est_true_value_treat

load("HPC output/est_true_value_conf_noint.rda")
load("HPC output/est_true_value_treat_noint.rda")

est_true_value_conf_noint <- est_true_value_conf
est_true_value_treat_noint <- est_true_value_treat

bootstrap_coefs <- array(,dim = c(10,2,1000,9))
sandwich_coefs <- array(,dim = c(10,2,1000,9))
bootstrap_treat <- array(,dim = c(10,2,1000,9))
sandwich_treat <- array(,dim = c(10,2,1000,9))
bootstrap_size <- array(,dim = c(10,2,1000,6))
sandwich_size <- array(,dim = c(10,2,1000,6))

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

sizes <- c(2,5,8,10,25,50)
for (i in 1:6){
  load(paste0("HPC output/CI_bootstrap_size_PP_", sizes[i], ".rda"))
  load(paste0("HPC output/CI_sandwich_size_PP_", sizes[i], ".rda"))
  bootstrap_size[,,,i] <- CI_bootstrap_size_PP
  sandwich_size[,,,i] <- CI_sandwich_size_PP
}

bootstrap_coefs_noint <- array(,dim = c(10,2,1000,9))
sandwich_coefs_noint <- array(,dim = c(10,2,1000,9))
bootstrap_treat_noint <- array(,dim = c(10,2,1000,9))
sandwich_treat_noint <- array(,dim = c(10,2,1000,9))
bootstrap_size_noint <- array(,dim = c(10,2,1000,6))
sandwich_size_noint <- array(,dim = c(10,2,1000,6))

for (i in 1:9){
  load(paste0("HPC output/CI_bootstrap_coefs_PP_noint_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_coefs_PP_noint_", i, ".rda"))
  load(paste0("HPC output/CI_bootstrap_treat_PP_noint_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_treat_PP_noint_", i, ".rda"))
  bootstrap_coefs_noint[,,,i] <- CI_bootstrap_coefs_PP
  sandwich_coefs_noint[,,,i] <- CI_sandwich_coefs_PP
  bootstrap_treat_noint[,,,i] <- CI_bootstrap_treat_PP
  sandwich_treat_noint[,,,i] <- CI_sandwich_treat_PP
}

for (i in 1:6){
  load(paste0("HPC output/CI_bootstrap_size_PP_noint_", sizes[i], ".rda"))
  load(paste0("HPC output/CI_sandwich_size_PP_noint_", sizes[i], ".rda"))
  bootstrap_size_noint[,,,i] <- CI_bootstrap_size_PP
  sandwich_size_noint[,,,i] <- CI_sandwich_size_PP
}


mean_lengths <- matrix(0, nrow = 8, ncol = 9)

for (i in 1:9){
  mean_lengths[1,i] <- mean(rowMeans(bootstrap_coefs[,2,,i]- bootstrap_coefs[,1,,i], na.rm = TRUE))
  mean_lengths[2,i] <- mean(rowMeans(sandwich_coefs[,2,,i]- sandwich_coefs[,1,,i], na.rm = TRUE))
  mean_lengths[3,i] <- mean(rowMeans(bootstrap_treat[,2,,i]- bootstrap_treat[,1,,i], na.rm = TRUE))
  mean_lengths[4,i] <- mean(rowMeans(sandwich_treat[,2,,i]- sandwich_treat[,1,,i], na.rm = TRUE))
  mean_lengths[5,i] <- mean(rowMeans(bootstrap_coefs_noint[,2,,i]- bootstrap_coefs_noint[,1,,i], na.rm = TRUE))
  mean_lengths[6,i] <- mean(rowMeans(sandwich_coefs_noint[,2,,i]- sandwich_coefs_noint[,1,,i], na.rm = TRUE))
  mean_lengths[7,i] <- mean(rowMeans(bootstrap_treat_noint[,2,,i]- bootstrap_treat_noint[,1,,i], na.rm = TRUE))
  mean_lengths[8,i] <- mean(rowMeans(sandwich_treat_noint[,2,,i]- sandwich_treat_noint[,1,,i], na.rm = TRUE))
}

mean_lengths_size <- matrix(0,nrow = 4, ncol = 6)
for (i in 1:6){
  mean_lengths_size[1,i] <- mean(rowMeans(bootstrap_size[,2,,i]- bootstrap_size[,1,,i], na.rm = TRUE))
  mean_lengths_size[2,i] <- mean(rowMeans(sandwich_size[,2,,i]- sandwich_size[,1,,i], na.rm = TRUE))
  mean_lengths_size[3,i] <- mean(rowMeans(bootstrap_size_noint[,2,,i]- bootstrap_size_noint[,1,,i], na.rm = TRUE))
  mean_lengths_size[4,i] <- mean(rowMeans(sandwich_size_noint[,2,,i]- sandwich_size_noint[,1,,i], na.rm = TRUE))
}

p1 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_lengths[1,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths[2,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[2,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Confounding strength', 
       y = "Mean CI length",
       title = "Outcome model with interaction terms")

p2 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_lengths[3,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[3,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths[4,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[4,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Prevalence of treatment', 
       y = "Mean CI length",
       title = "Outcome model with interaction terms")

p3 <- ggplot() +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[1,], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = mean_lengths_size[1,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[2,], colour = "Sandwich")) +
  geom_point(aes(x =sizes*100, y = mean_lengths_size[2,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Sample size', 
       y = "Mean CI length",
       title = "Outcome model with interaction terms")

p4 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_lengths[5,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[5,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths[6,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[6,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Confounding strength', 
       y = "Mean CI length",
       title = "Outcome model with no interaction term")

p5 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_lengths[7,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[7,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths[8,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[8,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Prevalence of treatment', 
       y = "Mean CI length",
       title = "Outcome model with no interaction term")

p6 <- ggplot() +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[3,], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = mean_lengths_size[3,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[4,], colour = "Sandwich")) +
  geom_point(aes(x =sizes*100, y = mean_lengths_size[4,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Sample size', 
       y = "Mean CI length",
       title = "Outcome model with no interaction term")
ggarrange(p1, p2,p3,p4,p5,p6, ncol=3, nrow=2, common.legend = TRUE, legend="right")

coverage_ind <- matrix(0,nrow = 8, ncol = 9)
coverage_ind_est <- matrix(0,nrow = 8, ncol = 9)

coverage_ind_size <- matrix(0,nrow = 4, ncol = 6)
coverage_ind_est_size <- matrix(0,nrow = 4, ncol = 6)

success <- matrix(0,nrow = 8, ncol = 9)
success_size <- matrix(0,nrow = 4, ncol = 6)

for (i in 1:1000){
  for (j in 1:9){
    
    if (all(is.na(bootstrap_coefs[10,1,i,j])) == F){
      success[1,j] <- success[1,j] + 1
      if (all(bootstrap_coefs[10,1,i,j] < true_value_conf[10,1,j] - true_value_conf[10,2,j]) 
          & all(bootstrap_coefs[10,2,i,j] > true_value_conf[10,1,j] - true_value_conf[10,2,j])){
        coverage_ind[1,j] <- coverage_ind[1,j] + 1
      }
      if (all(bootstrap_coefs[10,1,i,j] < est_true_value_conf_full[10,1,j] - est_true_value_conf_full[10,2,j]) 
          & all(bootstrap_coefs[10,2,i,j] > est_true_value_conf_full[10,1,j] - est_true_value_conf_full[10,2,j])){
        coverage_ind_est[1,j] <- coverage_ind_est[1,j] + 1
      }
    }
    
    if (all(is.na(sandwich_coefs[10,1,i,j])) == F){
      success[2,j] <- success[2,j] + 1
      if (all(sandwich_coefs[10,1,i,j] < true_value_conf[10,1,j] - true_value_conf[10,2,j]) 
        & all(sandwich_coefs[10,2,i,j] > true_value_conf[10,1,j] - true_value_conf[10,2,j])){
        coverage_ind[2,j] <- coverage_ind[2,j] + 1
      }
      if (all(sandwich_coefs[10,1,i,j] < est_true_value_conf_full[10,1,j] - est_true_value_conf_full[10,2,j]) 
          & all(sandwich_coefs[10,2,i,j] > est_true_value_conf_full[10,1,j] - est_true_value_conf_full[10,2,j])){
        coverage_ind_est[2,j] <- coverage_ind_est[2,j] + 1
      }
    }
    if (all(is.na(bootstrap_treat[10,1,i,j])) == F){
      success[3,j] <- success[3,j] + 1
      if (all(bootstrap_treat[10,1,i,j] < true_value_treat[10,1,j] - true_value_treat[10,2,j]) 
          & all(bootstrap_treat[10,2,i,j] > true_value_treat[10,1,j] - true_value_treat[10,2,j])){
        coverage_ind[3,j] <- coverage_ind[3,j] + 1
      }
      if (all(bootstrap_treat[10,1,i,j] < est_true_value_treat_full[10,1,j] - est_true_value_treat_full[10,2,j]) 
          & all(bootstrap_treat[10,2,i,j] > est_true_value_treat_full[10,1,j] - est_true_value_treat_full[10,2,j])){
        coverage_ind_est[3,j] <- coverage_ind_est[3,j] + 1
      }
    }
    if (all(is.na(sandwich_treat[10,1,i,j])) == F){
      success[4,j] <- success[4,j] + 1
      if (all(sandwich_treat[10,1,i,j] < true_value_treat[10,1,j] - true_value_treat[10,2,j]) 
          & all(sandwich_treat[10,2,i,j] > true_value_treat[10,1,j] - true_value_treat[10,2,j])){
        coverage_ind[4,j] <- coverage_ind[4,j] + 1
      }
      if (all(sandwich_treat[10,1,i,j] < est_true_value_treat_full[10,1,j] - est_true_value_treat_full[10,2,j]) 
          & all(sandwich_treat[10,2,i,j] > est_true_value_treat_full[10,1,j] - est_true_value_treat_full[10,2,j])){
        coverage_ind_est[4,j] <- coverage_ind_est[4,j] + 1
      }
    }
    
    if (all(is.na(bootstrap_coefs_noint[10,1,i,j])) == F){
      success[5,j] <- success[5,j] + 1
      if (all(bootstrap_coefs_noint[10,1,i,j] < true_value_conf[10,1,j] - true_value_conf[10,2,j]) 
          & all(bootstrap_coefs_noint[10,2,i,j] > true_value_conf[10,1,j] - true_value_conf[10,2,j])){
        coverage_ind[5,j] <- coverage_ind[5,j] + 1
      }
      if (all(bootstrap_coefs_noint[10,1,i,j] < est_true_value_conf_noint[10,1,j] - est_true_value_conf_noint[10,2,j]) 
          & all(bootstrap_coefs_noint[10,2,i,j] > est_true_value_conf_noint[10,1,j] - est_true_value_conf_noint[10,2,j])){
        coverage_ind_est[5,j] <- coverage_ind_est[5,j] + 1
      }
    }
    
    if (all(is.na(sandwich_coefs_noint[10,1,i,j])) == F){
      success[6,j] <- success[6,j] + 1
      if (all(sandwich_coefs_noint[10,1,i,j] < true_value_conf[10,1,j] - true_value_conf[10,2,j]) 
          & all(sandwich_coefs_noint[10,2,i,j] > true_value_conf[10,1,j] - true_value_conf[10,2,j])){
        coverage_ind[6,j] <- coverage_ind[6,j] + 1
      }
      if (all(sandwich_coefs_noint[10,1,i,j] < est_true_value_conf_noint[10,1,j] - est_true_value_conf_noint[10,2,j]) 
          & all(sandwich_coefs_noint[10,2,i,j] > est_true_value_conf_noint[10,1,j] - est_true_value_conf_noint[10,2,j])){
        coverage_ind_est[6,j] <- coverage_ind_est[6,j] + 1
      }
    }
    if (all(is.na(bootstrap_treat_noint[10,1,i,j])) == F){
      success[7,j] <- success[7,j] + 1
      if (all(bootstrap_treat_noint[10,1,i,j] < true_value_treat[10,1,j] - true_value_treat[10,2,j]) 
          & all(bootstrap_treat_noint[10,2,i,j] > true_value_treat[10,1,j] - true_value_treat[10,2,j])){
        coverage_ind[7,j] <- coverage_ind[7,j] + 1
      }
      if (all(bootstrap_treat_noint[10,1,i,j] < est_true_value_treat_noint[10,1,j] - est_true_value_treat_noint[10,2,j]) 
          & all(bootstrap_treat_noint[10,2,i,j] > est_true_value_treat_noint[10,1,j] - est_true_value_treat_noint[10,2,j])){
        coverage_ind_est[7,j] <- coverage_ind_est[7,j] + 1
      }
    }
    if (all(is.na(sandwich_treat_noint[10,1,i,j])) == F){
      success[8,j] <- success[8,j] + 1
      if (all(sandwich_treat_noint[10,1,i,j] < true_value_treat[10,1,j] - true_value_treat[10,2,j]) 
          & all(sandwich_treat_noint[10,2,i,j] > true_value_treat[10,1,j] - true_value_treat[10,2,j])){
        coverage_ind[8,j] <- coverage_ind[8,j] + 1
      }
      if (all(sandwich_treat_noint[10,1,i,j] < est_true_value_treat_noint[10,1,j] - est_true_value_treat_noint[10,2,j]) 
          & all(sandwich_treat_noint[10,2,i,j] > est_true_value_treat_noint[10,1,j] - est_true_value_treat_noint[10,2,j])){
        coverage_ind_est[8,j] <- coverage_ind_est[8,j] + 1
      }
    }
  }
  for (j in 1:6){
    
    if (all(is.na(bootstrap_size[10,1,i,j])) == F){
      success_size[1,j] <- success_size[1,j] + 1
      if (all(bootstrap_size[10,1,i,j] < true_value_conf[10,1,5] - true_value_conf[10,2,5]) 
          & all(bootstrap_size[10,2,i,j] > true_value_conf[10,1,j] - true_value_conf[10,2,j])){
        coverage_ind_size[1,j] <- coverage_ind_size[1,j] + 1
      }
      if (all(bootstrap_size[10,1,i,j] < est_true_value_conf_full[10,1,5] - est_true_value_conf_full[10,2,5]) 
          & all(bootstrap_size[10,2,i,j] > est_true_value_conf_full[10,1,5] - est_true_value_conf_full[10,2,5])){
        coverage_ind_est_size[1,j] <- coverage_ind_est_size[1,j] + 1
      }
    }
    
    if (all(is.na(sandwich_size[10,1,i,j])) == F){
      success_size[2,j] <- success_size[2,j] + 1
      if (all(sandwich_size[10,1,i,j] < true_value_conf[10,1,5] - true_value_conf[10,2,5]) 
          & all(sandwich_size[10,2,i,j] > true_value_conf[10,1,5] - true_value_conf[10,2,5])){
        coverage_ind_size[2,j] <- coverage_ind_size[2,j] + 1
      }
      if (all(sandwich_size[10,1,i,j] < est_true_value_conf_full[10,1,5] - est_true_value_conf_full[10,2,5]) 
          & all(sandwich_size[10,2,i,j] > est_true_value_conf_full[10,1,5] - est_true_value_conf_full[10,2,5])){
        coverage_ind_est_size[2,j] <- coverage_ind_est_size[2,j] + 1
      }
    }
    
    if (all(is.na(bootstrap_size_noint[10,1,i,j])) == F){
      success_size[3,j] <- success_size[3,j] + 1
      if (all(bootstrap_size_noint[10,1,i,j] < true_value_conf[10,1,5] - true_value_conf[10,2,5]) 
          & all(bootstrap_size_noint[10,2,i,j] > true_value_conf[10,1,j] - true_value_conf[10,2,j])){
        coverage_ind_size[3,j] <- coverage_ind_size[3,j] + 1
      }
      if (all(bootstrap_size_noint[10,1,i,j] < est_true_value_conf_noint[10,1,5] - est_true_value_conf_noint[10,2,5]) 
          & all(bootstrap_size_noint[10,2,i,j] > est_true_value_conf_noint[10,1,5] - est_true_value_conf_noint[10,2,5])){
        coverage_ind_est_size[3,j] <- coverage_ind_est_size[3,j] + 1
      }
    }
    
    if (all(is.na(sandwich_size_noint[10,1,i,j])) == F){
      success_size[4,j] <- success_size[4,j] + 1
      if (all(sandwich_size_noint[10,1,i,j] < true_value_conf[10,1,5] - true_value_conf[10,2,5]) 
          & all(sandwich_size_noint[10,2,i,j] > true_value_conf[10,1,5] - true_value_conf[10,2,5])){
        coverage_ind_size[4,j] <- coverage_ind_size[4,j] + 1
      }
      if (all(sandwich_size_noint[10,1,i,j] < est_true_value_conf_noint[10,1,5] - est_true_value_conf_noint[10,2,5]) 
          & all(sandwich_size_noint[10,2,i,j] > est_true_value_conf_noint[10,1,5] - est_true_value_conf_noint[10,2,5])){
        coverage_ind_est_size[4,j] <- coverage_ind_est_size[4,j] + 1
      }
    }
  }
}

coverage_ind <- coverage_ind/success
coverage_ind_est <- coverage_ind_est/success
coverage_ind_size <- coverage_ind_size/success_size
coverage_ind_est_size <- coverage_ind_est_size/success_size

p7 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = coverage_ind[1,], colour = "Bootstrap",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind[2,], colour = "Sandwich",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[2,], colour = "Sandwich")) +
  
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[1,], colour = "Bootstrap",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[2,], colour = "Sandwich",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[2,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  scale_linetype_manual(name = "True value estimate", values = c("Kaplan-Meier"= "solid", "Extreme sample" = "dashed")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Confounding strength', 
       y = "Empirical coverage rate",
       title = "Outcome model with interaction terms")

p8 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = coverage_ind[3,], colour = "Bootstrap",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[3,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind[4,], colour = "Sandwich",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[4,], colour = "Sandwich")) +
  
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[3,], colour = "Bootstrap",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[3,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[4,], colour = "Sandwich",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[4,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  scale_linetype_manual(name = "True value estimate", values = c("Kaplan-Meier"= "solid", "Extreme sample" = "dashed")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Prevalence of treatment', 
       y = "Empirical coverage rate",
       title = "Outcome model with interaction terms")

p9 <- ggplot() +
  geom_line(aes(x = sizes*100, y = coverage_ind_size[1,], colour = "Bootstrap",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = sizes*100, y = coverage_ind_size[1,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = coverage_ind_size[2,], colour = "Sandwich",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = sizes*100, y = coverage_ind_size[2,], colour = "Sandwich")) +
  
  geom_line(aes(x = sizes*100, y = coverage_ind_est_size[1,], colour = "Bootstrap",linetype = "Extreme sample")) +
  geom_point(aes(x = sizes*100, y = coverage_ind_est_size[1,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = coverage_ind_est_size[2,], colour = "Sandwich",linetype = "Extreme sample")) +
  geom_point(aes(x = sizes*100, y = coverage_ind_est_size[2,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  scale_linetype_manual(name = "True value estimate", values = c("Kaplan-Meier"= "solid", "Extreme sample" = "dashed")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Sample size', 
       y = "Empirical coverage rate",
       title = "Outcome model with interaction terms")

p10 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = coverage_ind[5,], colour = "Bootstrap",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[5,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind[6,], colour = "Sandwich",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[6,], colour = "Sandwich")) +
  
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[5,], colour = "Bootstrap",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[5,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[6,], colour = "Sandwich",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[6,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  scale_linetype_manual(name = "True value estimate", values = c("Kaplan-Meier"= "solid", "Extreme sample" = "dashed")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Confounding strength', 
       y = "Empirical coverage rate",
       title = "Outcome model with no interaction term")

p11 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = coverage_ind[7,], colour = "Bootstrap",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[7,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind[8,], colour = "Sandwich",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind[8,], colour = "Sandwich")) +
  
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[7,], colour = "Bootstrap",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[7,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = coverage_ind_est[8,], colour = "Sandwich",linetype = "Extreme sample")) +
  geom_point(aes(x = 1:9/10, y = coverage_ind_est[8,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  scale_linetype_manual(name = "True value estimate", values = c("Kaplan-Meier"= "solid", "Extreme sample" = "dashed")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Prevalence of treatment', 
       y = "Empirical coverage rate",
       title = "Outcome model with no interaction term")

p12 <- ggplot() +
  geom_line(aes(x = sizes*100, y = coverage_ind_size[3,], colour = "Bootstrap",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = sizes*100, y = coverage_ind_size[3,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = coverage_ind_size[4,], colour = "Sandwich",linetype = "Kaplan-Meier")) +
  geom_point(aes(x = sizes*100, y = coverage_ind_size[4,], colour = "Sandwich")) +
  
  geom_line(aes(x = sizes*100, y = coverage_ind_est_size[3,], colour = "Bootstrap",linetype = "Extreme sample")) +
  geom_point(aes(x = sizes*100, y = coverage_ind_est_size[3,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = coverage_ind_est_size[4,], colour = "Sandwich",linetype = "Extreme sample")) +
  geom_point(aes(x = sizes*100, y = coverage_ind_est_size[4,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  scale_linetype_manual(name = "True value estimate", values = c("Kaplan-Meier"= "solid", "Extreme sample" = "dashed")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Sample size', 
       y = "Empirical coverage rate",
       title = "Outcome model with no interaction term")
ggarrange(p7, p8,p9,p10,p11,p12, ncol=3, nrow=2, common.legend = TRUE, legend="right")
