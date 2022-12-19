library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/MPhil-dissertation/Code")
library(ggplot2)

load("HPC output/true_value_conf_red.rda")
load("HPC output/true_value_treat_red.rda")

bootstrap_coefs <- array(,dim = c(5,2,1000,9))
LEF_outcome_coefs <- array(,dim = c(5,2,1000,9))
LEF_both_coefs <- array(,dim = c(5,2,1000,9))
sandwich_coefs <- array(,dim = c(5,2,1000,9))
time_coefs <- array(,dim = c(4,1000,9))

bootstrap_treat <- array(,dim = c(5,2,1000,9))
LEF_outcome_treat <- array(,dim = c(5,2,1000,9))
LEF_both_treat <- array(,dim = c(5,2,1000,9))
sandwich_treat <- array(,dim = c(5,2,1000,9))
time_treat <- array(,dim = c(4,1000,9))

bootstrap_size <- array(,dim = c(5,2,1000,7))
LEF_outcome_size <- array(,dim = c(5,2,1000,7))
LEF_both_size <- array(,dim = c(5,2,1000,7))
sandwich_size <- array(,dim = c(5,2,1000,7))
time_size <- array(,dim = c(4,1000,7))

bootstrap_coefs_big <- array(,dim = c(5,2,1000,9))
LEF_outcome_coefs_big <- array(,dim = c(5,2,1000,9))
LEF_both_coefs_big <- array(,dim = c(5,2,1000,9))
sandwich_coefs_big <- array(,dim = c(5,2,1000,9))
time_coefs_big <- array(,dim = c(4,1000,9))

bootstrap_treat_big <- array(,dim = c(5,2,1000,9))
LEF_outcome_treat_big <- array(,dim = c(5,2,1000,9))
LEF_both_treat_big <- array(,dim = c(5,2,1000,9))
sandwich_treat_big <- array(,dim = c(5,2,1000,9))
time_treat_big <-  array(,dim = c(4,1000,9))

for (i in 1:9){
  load(paste0("HPC output/CI_bootstrap_coefs_PP_red_", i, ".rda"))
  load(paste0("HPC output/CI_LEF_outcome_coefs_PP_red_", i, ".rda"))
  load(paste0("HPC output/CI_LEF_both_coefs_PP_red_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_coefs_PP_red_", i, ".rda"))
  load(paste0("HPC output/computation_time_coefs_", i, ".rda"))
  
  load(paste0("HPC output/CI_bootstrap_treat_PP_red_", i, ".rda"))
  load(paste0("HPC output/CI_LEF_outcome_treat_PP_red_", i, ".rda"))
  load(paste0("HPC output/CI_LEF_both_treat_PP_red_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_treat_PP_red_", i, ".rda"))
  load(paste0("HPC output/computation_time_treat_", i, ".rda"))
  
  bootstrap_coefs[,,,i] <- CI_bootstrap_coefs_PP_red
  LEF_outcome_coefs[,,,i] <- CI_LEF_outcome_coefs_PP_red
  LEF_both_coefs[,,,i] <- CI_LEF_both_coefs_PP_red
  sandwich_coefs[,,,i] <- CI_sandwich_coefs_PP_red
  time_coefs[,,i] <- computation_time_coefs
  
  bootstrap_treat[,,,i] <- CI_bootstrap_treat_PP_red
  LEF_outcome_treat[,,,i] <- CI_LEF_outcome_treat_PP_red
  LEF_both_treat[,,,i] <- CI_LEF_both_treat_PP_red
  sandwich_treat[,,,i] <- CI_sandwich_treat_PP_red
  time_treat[,,i] <- computation_time_treat
  
  load(paste0("HPC output/CI_bootstrap_coefs_PP_red_big_", i, ".rda"))
  load(paste0("HPC output/CI_LEF_outcome_coefs_PP_red_big_", i, ".rda"))
  load(paste0("HPC output/CI_LEF_both_coefs_PP_red_big_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_coefs_PP_red_big_", i, ".rda"))
  load(paste0("HPC output/computation_time_coefs_big_", i, ".rda"))
  
  load(paste0("HPC output/CI_bootstrap_treat_PP_red_big_", i, ".rda"))
  load(paste0("HPC output/CI_LEF_outcome_treat_PP_red_big_", i, ".rda"))
  load(paste0("HPC output/CI_LEF_both_treat_PP_red_big_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_treat_PP_red_big_", i, ".rda"))
  load(paste0("HPC output/computation_time_treat_big_", i, ".rda"))
  
  bootstrap_coefs_big[,,,i] <- CI_bootstrap_coefs_PP_red_big
  LEF_outcome_coefs_big[,,,i] <- CI_LEF_outcome_coefs_PP_red_big
  LEF_both_coefs_big[,,,i] <- CI_LEF_both_coefs_PP_red_big
  sandwich_coefs_big[,,,i] <- CI_sandwich_coefs_PP_red_big
  time_coefs_big[,,i] <- computation_time_coefs_big
  
  bootstrap_treat_big[,,,i] <- CI_bootstrap_treat_PP_red_big
  LEF_outcome_treat_big[,,,i] <- CI_LEF_outcome_treat_PP_red_big
  LEF_both_treat_big[,,,i] <- CI_LEF_both_treat_PP_red_big
  sandwich_treat_big[,,,i] <- CI_sandwich_treat_PP_red_big
  time_treat_big[,,i] <- computation_time_treat_big
  
}

sizes <- c(2,5,8,10,25,50,100)
for (i in 1:7){
  load(paste0("HPC output/CI_bootstrap_size_PP_red_", sizes[i], ".rda"))
  load(paste0("HPC output/CI_LEF_outcome_size_PP_red_", sizes[i], ".rda"))
  load(paste0("HPC output/CI_LEF_both_size_PP_red_", sizes[i], ".rda"))
  load(paste0("HPC output/CI_sandwich_size_PP_red_", sizes[i], ".rda"))
  load(paste0("HPC output/computation_time_size_", sizes[i], ".rda"))
  
  bootstrap_size[,,,i] <- CI_bootstrap_size_PP_red
  LEF_outcome_size[,,,i] <- CI_LEF_outcome_size_PP_red
  LEF_both_size[,,,i] <- CI_LEF_both_size_PP_red
  sandwich_size[,,,i] <- CI_sandwich_size_PP_red
  time_size[,,i] <- computation_time_size
}

mean_lengths_coefs <- matrix(0, nrow = 4, ncol = 9)
mean_lengths_treat <- matrix(0, nrow = 8, ncol = 9)
mean_lengths_size <- matrix(0, nrow = 8, ncol = 7)

sd_lengths_coefs <- matrix(0, nrow = 4, ncol = 9)
sd_lengths_treat <- matrix(0, nrow = 8, ncol = 9)
sd_lengths_size <- matrix(0, nrow = 8, ncol = 7)

mean_lengths_coefs_big <- matrix(0, nrow = 4, ncol = 9)
mean_lengths_treat_big <- matrix(0, nrow = 8, ncol = 9)

sd_lengths_coefs_big <- matrix(0, nrow = 4, ncol = 9)
sd_lengths_treat_big <- matrix(0, nrow = 8, ncol = 9)

for (i in 1:9){
  mean_lengths_coefs[1,i] <- mean(rowMeans(bootstrap_coefs[,2,,i]- bootstrap_coefs[,1,,i], na.rm = TRUE))
  mean_lengths_coefs[2,i] <- mean(rowMeans(LEF_outcome_coefs[,2,,i]- LEF_outcome_coefs[,1,,i], na.rm = TRUE))
  mean_lengths_coefs[3,i] <- mean(rowMeans(LEF_both_coefs[,2,,i]- LEF_both_coefs[,1,,i], na.rm = TRUE))
  mean_lengths_coefs[4,i] <- mean(rowMeans(sandwich_coefs[,2,,i]- sandwich_coefs[,1,,i], na.rm = TRUE))
  
  mean_lengths_treat[1,i] <- mean(rowMeans(bootstrap_treat[,2,,i]- bootstrap_treat[,1,,i], na.rm = TRUE))
  mean_lengths_treat[2,i] <- mean(rowMeans(LEF_outcome_treat[,2,,i]- LEF_outcome_treat[,1,,i], na.rm = TRUE))
  mean_lengths_treat[3,i] <- mean(rowMeans(LEF_both_treat[,2,,i]- LEF_both_treat[,1,,i], na.rm = TRUE))
  mean_lengths_treat[4,i] <- mean(rowMeans(sandwich_treat[,2,,i]- sandwich_treat[,1,,i], na.rm = TRUE))
  
  sd_lengths_coefs[1,i] <- sd(rowMeans(bootstrap_coefs[,2,,i]- bootstrap_coefs[,1,,i], na.rm = TRUE))
  sd_lengths_coefs[2,i] <- sd(rowMeans(LEF_outcome_coefs[,2,,i]- LEF_outcome_coefs[,1,,i], na.rm = TRUE))
  sd_lengths_coefs[3,i] <- sd(rowMeans(LEF_both_coefs[,2,,i]- LEF_both_coefs[,1,,i], na.rm = TRUE))
  sd_lengths_coefs[4,i] <- sd(rowMeans(sandwich_coefs[,2,,i]- sandwich_coefs[,1,,i], na.rm = TRUE))
  
  sd_lengths_treat[1,i] <- sd(rowMeans(bootstrap_treat[,2,,i]- bootstrap_treat[,1,,i], na.rm = TRUE))
  sd_lengths_treat[2,i] <- sd(rowMeans(LEF_outcome_treat[,2,,i]- LEF_outcome_treat[,1,,i], na.rm = TRUE))
  sd_lengths_treat[3,i] <- sd(rowMeans(LEF_both_treat[,2,,i]- LEF_both_treat[,1,,i], na.rm = TRUE))
  sd_lengths_treat[4,i] <- sd(rowMeans(sandwich_treat[,2,,i]- sandwich_treat[,1,,i], na.rm = TRUE))
  
  mean_lengths_coefs_big[1,i] <- mean(rowMeans(bootstrap_coefs_big[,2,,i]- bootstrap_coefs_big[,1,,i], na.rm = TRUE))
  mean_lengths_coefs_big[2,i] <- mean(rowMeans(LEF_outcome_coefs_big[,2,,i]- LEF_outcome_coefs_big[,1,,i], na.rm = TRUE))
  mean_lengths_coefs_big[3,i] <- mean(rowMeans(LEF_both_coefs_big[,2,,i]- LEF_both_coefs_big[,1,,i], na.rm = TRUE))
  mean_lengths_coefs_big[4,i] <- mean(rowMeans(sandwich_coefs_big[,2,,i]- sandwich_coefs_big[,1,,i], na.rm = TRUE))
  
  mean_lengths_treat_big[1,i] <- mean(rowMeans(bootstrap_treat_big[,2,,i]- bootstrap_treat_big[,1,,i], na.rm = TRUE))
  mean_lengths_treat_big[2,i] <- mean(rowMeans(LEF_outcome_treat_big[,2,,i]- LEF_outcome_treat_big[,1,,i], na.rm = TRUE))
  mean_lengths_treat_big[3,i] <- mean(rowMeans(LEF_both_treat_big[,2,,i]- LEF_both_treat_big[,1,,i], na.rm = TRUE))
  mean_lengths_treat_big[4,i] <- mean(rowMeans(sandwich_treat_big[,2,,i]- sandwich_treat_big[,1,,i], na.rm = TRUE))
  
  sd_lengths_coefs_big[1,i] <- sd(rowMeans(bootstrap_coefs_big[,2,,i]- bootstrap_coefs_big[,1,,i], na.rm = TRUE))
  sd_lengths_coefs_big[2,i] <- sd(rowMeans(LEF_outcome_coefs_big[,2,,i]- LEF_outcome_coefs_big[,1,,i], na.rm = TRUE))
  sd_lengths_coefs_big[3,i] <- sd(rowMeans(LEF_both_coefs_big[,2,,i]- LEF_both_coefs_big[,1,,i], na.rm = TRUE))
  sd_lengths_coefs_big[4,i] <- sd(rowMeans(sandwich_coefs_big[,2,,i]- sandwich_coefs_big[,1,,i], na.rm = TRUE))
  
  sd_lengths_treat_big[1,i] <- sd(rowMeans(bootstrap_treat_big[,2,,i]- bootstrap_treat_big[,1,,i], na.rm = TRUE))
  sd_lengths_treat_big[2,i] <- sd(rowMeans(LEF_outcome_treat_big[,2,,i]- LEF_outcome_treat_big[,1,,i], na.rm = TRUE))
  sd_lengths_treat_big[3,i] <- sd(rowMeans(LEF_both_treat_big[,2,,i]- LEF_both_treat_big[,1,,i], na.rm = TRUE))
  sd_lengths_treat_big[4,i] <- sd(rowMeans(sandwich_treat_big[,2,,i]- sandwich_treat_big[,1,,i], na.rm = TRUE))
  
}

for (i in 1:7){
  mean_lengths_size[1,i] <- mean(rowMeans(bootstrap_size[,2,,i]- bootstrap_size[,1,,i], na.rm = TRUE))
  mean_lengths_size[2,i] <- mean(rowMeans(LEF_outcome_size[,2,,i]- LEF_outcome_size[,1,,i], na.rm = TRUE))
  mean_lengths_size[3,i] <- mean(rowMeans(LEF_both_size[,2,,i]- LEF_both_size[,1,,i], na.rm = TRUE))
  mean_lengths_size[4,i] <- mean(rowMeans(sandwich_size[,2,,i]- sandwich_size[,1,,i], na.rm = TRUE))
  
  sd_lengths_size[1,i] <- sd(rowMeans(bootstrap_size[,2,,i]- bootstrap_size[,1,,i], na.rm = TRUE))
  sd_lengths_size[2,i] <- sd(rowMeans(LEF_outcome_size[,2,,i]- LEF_outcome_size[,1,,i], na.rm = TRUE))
  sd_lengths_size[3,i] <- sd(rowMeans(LEF_both_size[,2,,i]- LEF_both_size[,1,,i], na.rm = TRUE))
  sd_lengths_size[4,i] <- sd(rowMeans(sandwich_size[,2,,i]- sandwich_size[,1,,i], na.rm = TRUE))
}

p1 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_lengths[1,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths[2,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[2,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Confounding strength', 
       y = "Mean CI length")

p2 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_lengths[3,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[3,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths[4,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[4,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Prevalence of treatment', 
       y = "Mean CI length")

p3 <- ggplot() +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[1,], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = mean_lengths_size[1,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[2,], colour = "Sandwich")) +
  geom_point(aes(x =sizes*100, y = mean_lengths_size[2,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Sample size', 
       y = "Mean CI length")

figureA <- ggarrange(p1 + rremove("ylab"),p2 + rremove("ylab"),p3 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figureA <- annotate_figure(figureA,top = text_grob("Full model",size = 14), left = text_grob("Mean CI length",size = 10, rot =  90))

p4 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_lengths[5,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[5,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths[6,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[6,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Confounding strength', 
       y = "Mean CI length")

p5 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_lengths[7,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[7,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths[8,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths[8,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Prevalence of treatment', 
       y = "Mean CI length")

p6 <- ggplot() +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[3,], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = mean_lengths_size[3,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[4,], colour = "Sandwich")) +
  geom_point(aes(x =sizes*100, y = mean_lengths_size[4,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Sample size', 
       y = "Mean CI length")
figureB <- ggarrange(p4 + rremove("ylab"),p5 + rremove("ylab"),p6 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figureB <- annotate_figure(figureB,top = text_grob("Simple model",size = 14), left = text_grob("Mean CI length",size = 10, rot =  90))
ggarrange(heights = c(4, 4.9),figureA, figureB, nrow = 2, ncol = 1)

p1_sd <- ggplot() +
  geom_line(aes(x = 1:9/10, y = sd_lengths[1,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = sd_lengths[2,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths[2,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Confounding strength', 
       y = "CI length SE") +ylim(0,0.2)

p2_sd <- ggplot() +
  geom_line(aes(x = 1:9/10, y = sd_lengths[3,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths[3,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = sd_lengths[4,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths[4,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Prevalence of treatment', 
       y = "Mean CI length") +ylim(0,0.2)

p3_sd <- ggplot() +
  geom_line(aes(x = sizes*100, y = sd_lengths_size[1,], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = sd_lengths_size[1,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = sd_lengths_size[2,], colour = "Sandwich")) +
  geom_point(aes(x =sizes*100, y = sd_lengths_size[2,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Sample size', 
       y = "Mean CI length")+ylim(0,0.2)

figureA_sd <- ggarrange(p1_sd + rremove("ylab"),p2_sd + rremove("ylab"),p3_sd + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figureA_sd <- annotate_figure(figureA_sd,top = text_grob("Full model",size = 14), left = text_grob("CI length SE",size = 10, rot =  90))

p4_sd <- ggplot() +
  geom_line(aes(x = 1:9/10, y = sd_lengths[5,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths[5,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = sd_lengths[6,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths[6,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Confounding strength', 
       y = "Mean CI length")+ylim(0,0.2)

p5_sd <- ggplot() +
  geom_line(aes(x = 1:9/10, y = sd_lengths[7,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths[7,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = sd_lengths[8,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths[8,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Prevalence of treatment', 
       y = "Mean CI length")+ylim(0,0.2)

p6_sd <- ggplot() +
  geom_line(aes(x = sizes*100, y = sd_lengths_size[3,], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = sd_lengths_size[3,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = sd_lengths_size[4,], colour = "Sandwich")) +
  geom_point(aes(x =sizes*100, y = sd_lengths_size[4,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Sample size', 
       y = "Mean CI length")+ylim(0,0.2)
figureB_sd <- ggarrange(p4_sd + rremove("ylab"),p5_sd + rremove("ylab"),p6_sd + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figureB_sd <- annotate_figure(figureB_sd,top = text_grob("Simple model",size = 14), left = text_grob("CI length SE",size = 10, rot =  90))
ggarrange(heights = c(4, 4.9),figureA_sd, figureB_sd, nrow = 2, ncol = 1)


coverage_ind <- array(0,dim = c(8,9,10))
coverage_ind_est <- array(0,dim = c(8,9,10))

coverage_ind_size <- array(0,dim = c(4,6,10))
coverage_ind_est_size <- array(0,dim = c(4,6,10))

success <- array(0,dim = c(8,9,10))
success_size <- array(0,dim = c(4,6,10))

for (i in 1:1000){
  for (k in 1:10){
    for (j in 1:9){
      if (is.na(bootstrap_coefs[k,1,i,j]) == F){
        success[1,j,k] <- success[1,j,k] + 1
        if (all(bootstrap_coefs[k,1,i,j] < true_value_conf[k,1,j] - true_value_conf[k,2,j]) 
            & all(bootstrap_coefs[k,2,i,j] > true_value_conf[k,1,j] - true_value_conf[k,2,j])){
          coverage_ind[1,j,k] <- coverage_ind[1,j,k] + 1
        }
        if (all(bootstrap_coefs[k,1,i,j] < est_true_value_conf_full[k,1,j] - est_true_value_conf_full[k,2,j]) 
            & all(bootstrap_coefs[k,2,i,j] > est_true_value_conf_full[k,1,j] - est_true_value_conf_full[k,2,j])){
          coverage_ind_est[1,j,k] <- coverage_ind_est[1,j,k] + 1
        }
      }
      
      if (all(is.na(sandwich_coefs[k,1,i,j])) == F){
        success[2,j,k] <- success[2,j,k] + 1
        if (all(sandwich_coefs[k,1,i,j] < true_value_conf[k,1,j] - true_value_conf[k,2,j]) 
            & all(sandwich_coefs[k,2,i,j] > true_value_conf[k,1,j] - true_value_conf[k,2,j])){
          coverage_ind[2,j,k]<- coverage_ind[2,j,k]+ 1
        }
        if (all(sandwich_coefs[k,1,i,j] < est_true_value_conf_full[k,1,j] - est_true_value_conf_full[k,2,j]) 
            & all(sandwich_coefs[k,2,i,j] > est_true_value_conf_full[k,1,j] - est_true_value_conf_full[k,2,j])){
          coverage_ind_est[2,j,k]<- coverage_ind_est[2,j,k]+ 1
        }
      }
      if (all(is.na(bootstrap_treat[k,1,i,j])) == F){
        success[3,j,k]<- success[3,j,k]+ 1
        if (all(bootstrap_treat[k,1,i,j] < true_value_treat[k,1,j] - true_value_treat[k,2,j]) 
            & all(bootstrap_treat[k,2,i,j] > true_value_treat[k,1,j] - true_value_treat[k,2,j])){
          coverage_ind[3,j,k]<- coverage_ind[3,j,k]+ 1
        }
        if (all(bootstrap_treat[k,1,i,j] < est_true_value_treat_full[k,1,j] - est_true_value_treat_full[k,2,j]) 
            & all(bootstrap_treat[k,2,i,j] > est_true_value_treat_full[k,1,j] - est_true_value_treat_full[k,2,j])){
          coverage_ind_est[3,j,k]<- coverage_ind_est[3,j,k]+ 1
        }
      }
      if (all(is.na(sandwich_treat[k,1,i,j])) == F){
        success[4,j,k]<- success[4,j,k]+ 1
        if (all(sandwich_treat[k,1,i,j] < true_value_treat[k,1,j] - true_value_treat[k,2,j]) 
            & all(sandwich_treat[k,2,i,j] > true_value_treat[k,1,j] - true_value_treat[k,2,j])){
          coverage_ind[4,j,k]<- coverage_ind[4,j,k]+ 1
        }
        if (all(sandwich_treat[k,1,i,j] < est_true_value_treat_full[k,1,j] - est_true_value_treat_full[k,2,j]) 
            & all(sandwich_treat[k,2,i,j] > est_true_value_treat_full[k,1,j] - est_true_value_treat_full[k,2,j])){
          coverage_ind_est[4,j,k]<- coverage_ind_est[4,j,k]+ 1
        }
      }
      
      if (all(is.na(bootstrap_coefs_noint[k,1,i,j])) == F){
        success[5,j,k]<- success[5,j,k]+ 1
        if (all(bootstrap_coefs_noint[k,1,i,j] < true_value_conf[k,1,j] - true_value_conf[k,2,j]) 
            & all(bootstrap_coefs_noint[k,2,i,j] > true_value_conf[k,1,j] - true_value_conf[k,2,j])){
          coverage_ind[5,j,k]<- coverage_ind[5,j,k]+ 1
        }
        if (all(bootstrap_coefs_noint[k,1,i,j] < est_true_value_conf_noint[k,1,j] - est_true_value_conf_noint[k,2,j]) 
            & all(bootstrap_coefs_noint[k,2,i,j] > est_true_value_conf_noint[k,1,j] - est_true_value_conf_noint[k,2,j])){
          coverage_ind_est[5,j,k]<- coverage_ind_est[5,j,k]+ 1
        }
      }
      
      if (all(is.na(sandwich_coefs_noint[k,1,i,j])) == F){
        success[6,j,k]<- success[6,j,k]+ 1
        if (all(sandwich_coefs_noint[k,1,i,j] < true_value_conf[k,1,j] - true_value_conf[k,2,j]) 
            & all(sandwich_coefs_noint[k,2,i,j] > true_value_conf[k,1,j] - true_value_conf[k,2,j])){
          coverage_ind[6,j,k]<- coverage_ind[6,j,k]+ 1
        }
        if (all(sandwich_coefs_noint[k,1,i,j] < est_true_value_conf_noint[k,1,j] - est_true_value_conf_noint[k,2,j]) 
            & all(sandwich_coefs_noint[k,2,i,j] > est_true_value_conf_noint[k,1,j] - est_true_value_conf_noint[k,2,j])){
          coverage_ind_est[6,j,k]<- coverage_ind_est[6,j,k]+ 1
        }
      }
      if (all(is.na(bootstrap_treat_noint[k,1,i,j])) == F){
        success[7,j,k]<- success[7,j,k]+ 1
        if (all(bootstrap_treat_noint[k,1,i,j] < true_value_treat[k,1,j] - true_value_treat[k,2,j]) 
            & all(bootstrap_treat_noint[k,2,i,j] > true_value_treat[k,1,j] - true_value_treat[k,2,j])){
          coverage_ind[7,j,k]<- coverage_ind[7,j,k]+ 1
        }
        if (all(bootstrap_treat_noint[k,1,i,j] < est_true_value_treat_noint[k,1,j] - est_true_value_treat_noint[k,2,j]) 
            & all(bootstrap_treat_noint[k,2,i,j] > est_true_value_treat_noint[k,1,j] - est_true_value_treat_noint[k,2,j])){
          coverage_ind_est[7,j,k]<- coverage_ind_est[7,j,k]+ 1
        }
      }
      if (all(is.na(sandwich_treat_noint[k,1,i,j])) == F){
        success[8,j,k]<- success[8,j,k]+ 1
        if (all(sandwich_treat_noint[k,1,i,j] < true_value_treat[k,1,j] - true_value_treat[k,2,j]) 
            & all(sandwich_treat_noint[k,2,i,j] > true_value_treat[k,1,j] - true_value_treat[k,2,j])){
          coverage_ind[8,j,k]<- coverage_ind[8,j,k]+ 1
        }
        if (all(sandwich_treat_noint[k,1,i,j] < est_true_value_treat_noint[k,1,j] - est_true_value_treat_noint[k,2,j]) 
            & all(sandwich_treat_noint[k,2,i,j] > est_true_value_treat_noint[k,1,j] - est_true_value_treat_noint[k,2,j])){
          coverage_ind_est[8,j,k]<- coverage_ind_est[8,j,k]+ 1
        }
      }
    }
    for (j in 1:6){
      
      if (all(is.na(bootstrap_size[k,1,i,j])) == F){
        success_size[1,j,k]<- success_size[1,j,k]+ 1
        if (all(bootstrap_size[k,1,i,j] < true_value_conf[k,1,5] - true_value_conf[k,2,5]) 
            & all(bootstrap_size[k,2,i,j] > true_value_conf[k,1,5] - true_value_conf[k,2,5])){
          coverage_ind_size[1,j,k]<- coverage_ind_size[1,j,k]+ 1
        }
        if (all(bootstrap_size[k,1,i,j] < est_true_value_conf_full[k,1,5] - est_true_value_conf_full[k,2,5]) 
            & all(bootstrap_size[k,2,i,j] > est_true_value_conf_full[k,1,5] - est_true_value_conf_full[k,2,5])){
          coverage_ind_est_size[1,j,k]<- coverage_ind_est_size[1,j,k]+ 1
        }
      }
      
      if (all(is.na(sandwich_size[k,1,i,j])) == F){
        success_size[2,j,k]<- success_size[2,j,k]+ 1
        if (all(sandwich_size[k,1,i,j] < true_value_conf[k,1,5] - true_value_conf[k,2,5]) 
            & all(sandwich_size[k,2,i,j] > true_value_conf[k,1,5] - true_value_conf[k,2,5])){
          coverage_ind_size[2,j,k]<- coverage_ind_size[2,j,k]+ 1
        }
        if (all(sandwich_size[k,1,i,j] < est_true_value_conf_full[k,1,5] - est_true_value_conf_full[k,2,5]) 
            & all(sandwich_size[k,2,i,j] > est_true_value_conf_full[k,1,5] - est_true_value_conf_full[k,2,5])){
          coverage_ind_est_size[2,j,k]<- coverage_ind_est_size[2,j,k]+ 1
        }
      }
      
      if (all(is.na(bootstrap_size_noint[k,1,i,j])) == F){
        success_size[3,j,k]<- success_size[3,j,k]+ 1
        if (all(bootstrap_size_noint[k,1,i,j] < true_value_conf[k,1,5] - true_value_conf[k,2,5]) 
            & all(bootstrap_size_noint[k,2,i,j] > true_value_conf[k,1,5] - true_value_conf[k,2,5])){
          coverage_ind_size[3,j,k]<- coverage_ind_size[3,j,k]+ 1
        }
        if (all(bootstrap_size_noint[k,1,i,j] < est_true_value_conf_noint[k,1,5] - est_true_value_conf_noint[k,2,5]) 
            & all(bootstrap_size_noint[k,2,i,j] > est_true_value_conf_noint[k,1,5] - est_true_value_conf_noint[k,2,5])){
          coverage_ind_est_size[3,j,k]<- coverage_ind_est_size[3,j,k]+ 1
        }
      }
      
      if (all(is.na(sandwich_size_noint[k,1,i,j])) == F){
        success_size[4,j,k]<- success_size[4,j,k]+ 1
        if (all(sandwich_size_noint[k,1,i,j] < true_value_conf[k,1,5] - true_value_conf[k,2,5]) 
            & all(sandwich_size_noint[k,2,i,j] > true_value_conf[k,1,5] - true_value_conf[k,2,5])){
          coverage_ind_size[4,j,k]<- coverage_ind_size[4,j,k]+ 1
        }
        if (all(sandwich_size_noint[k,1,i,j] < est_true_value_conf_noint[k,1,5] - est_true_value_conf_noint[k,2,5]) 
            & all(sandwich_size_noint[k,2,i,j] > est_true_value_conf_noint[k,1,5] - est_true_value_conf_noint[k,2,5])){
          coverage_ind_est_size[4,j,k]<- coverage_ind_est_size[4,j,k]+ 1
        }
      }
    }
  }
}

coverage_ind <- coverage_ind/success
coverage_ind_est <- coverage_ind_est/success
coverage_ind_size <- coverage_ind_size/success_size
coverage_ind_est_size <- coverage_ind_est_size/success_size

p7 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[1,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[1,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[2,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[2,1,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.1")+ ylim(0.5,1) + theme(plot.title = element_text(size=10))
p8 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[1,5,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[1,5,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[2,5,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[2,5,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.5")+ ylim(0.5,1) + theme(plot.title = element_text(size=10))

p9 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[1,9,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[1,9,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[2,9,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[2,9,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.9")+ ylim(0.5,1) + theme(plot.title = element_text(size=10))

figure1 <- ggarrange(p7+ rremove("ylab"),p8+ rremove("ylab"),p9+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure1 <- annotate_figure(figure1,top = text_grob("Full model",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))

p10 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[5,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[5,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[6,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[6,1,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.1")+ ylim(0.5,1) + theme(plot.title = element_text(size=10))
p11 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[5,5,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[5,5,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[6,5,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[6,5,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.5") + ylim(0.5,1)+ theme(plot.title = element_text(size=10))

p12 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[5,9,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[5,9,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[6,9,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[6,9,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time',
       title = "Confounding = 0.9")+ ylim(0.5,1) + theme(plot.title = element_text(size=10)) 

figure2 <- ggarrange(p10 + rremove("ylab"),p11 + rremove("ylab"),p12 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T,legend = 'bottom')
figure2 <- annotate_figure(figure2,top = text_grob("Simple model",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))
ggarrange(heights = c(4, 4.9),figure1, figure2, nrow = 2, ncol = 1)
p13 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[3,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[3,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[4,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[4,1,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical CI coverage",
       title = "Prevalence of treatment = 0.1")+ ylim(0.5,1) + theme(plot.title = element_text(size=10))

p14 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[3,5,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[3,5,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[4,5,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[4,5,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical CI coverage",
       title = "Prevalence of treatment = 0.5")+ ylim(0.5,1) + theme(plot.title = element_text(size=10))
p15 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[3,9,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[3,9,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[4,9,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[4,9,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical CI coverage",
       title = "Prevalence of treatment = 0.9")+ ylim(0.5,1) + theme(plot.title = element_text(size=10))

figure3 <- ggarrange(p13 + rremove("ylab"),p14 + rremove("ylab"),p15 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure3 <- annotate_figure(figure3,top = text_grob("Full model",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))

p16 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[7,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[7,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[8,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[8,1,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical CI coverage",
       title = "Prevalence of treatment = 0.1")+ ylim(0.5,1) + theme(plot.title = element_text(size=10))

p17 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[7,5,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[7,5,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[8,5,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[8,5,], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical CI coverage",
       title = "Prevalence of treatment = 0.5")+ ylim(0.5,1) + theme(plot.title = element_text(size=10))
p18 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind[7,9,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind[7,9,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind[8,9,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind[8,9,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical CI coverage",
       title = "Prevalence of treatment = 0.9")+ ylim(0.5,1) + theme(plot.title = element_text(size=10))

figure4 <- ggarrange(p16 + rremove("ylab"),p17 + rremove("ylab"),p18 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figure4 <- annotate_figure(figure4,top = text_grob("Simple model",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))
ggarrange(heights = c(4, 4.9),figure3, figure4, nrow = 2, ncol = 1)

p19 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind_size[1,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[1,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind_size[2,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[2,1,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Sample size = 200")+ ylim(0,1) + theme(plot.title = element_text(size=10))
p20 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind_size[1,4,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[1,4,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind_size[2,4,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[2,4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Sample size = 1000")+ ylim(0,1) + theme(plot.title = element_text(size=10))
p21 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind_size[1,6,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[1,6,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind_size[2,6,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[2,6,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Sample size = 50000")+ ylim(0,1) + theme(plot.title = element_text(size=10))

figure5 <- ggarrange(p19 + rremove("ylab"),p20 + rremove("ylab"),p21 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure5 <- annotate_figure(figure5,top = text_grob("Full model",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))

p22 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind_size[3,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[3,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind_size[4,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[4,1,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Sample size = 200")+ ylim(0,1) + theme(plot.title = element_text(size=10))
p23 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind_size[3,4,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[3,4,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind_size[4,4,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[4,4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Sample size = 1000")+ ylim(0,1) + theme(plot.title = element_text(size=10))
p24 <- ggplot() +
  geom_line(aes(x = 0:9, y = coverage_ind_size[3,6,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[3,6,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = coverage_ind_size[4,6,], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = coverage_ind_size[4,6,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Sample size = 50000")+ ylim(0,1) + theme(plot.title = element_text(size=10))

figure6 <- ggarrange(p22 + rremove("ylab"),p23 + rremove("ylab"),p24 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figure6 <- annotate_figure(figure6,top = text_grob("Simple model",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))

ggarrange(heights = c(4, 4.9),figure5, figure6, nrow = 2, ncol = 1)

SE <- sqrt((coverage_ind*(1-coverage_ind))/1000)
SE_size <- sqrt((coverage_ind_size*(1-coverage_ind_size))/1000)

p25 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = SE[1,,1], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = SE[1,,1], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = SE[2,,1], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = SE[2,,1], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Confounding strength', 
       y = "Monte Carlo SE") +ylim(0,0.02)

p26 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = SE[3,,1], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = SE[3,,1], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = SE[4,,1], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = SE[4,,1], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Prevalence of treatment', 
       y = "Monte Carlo SE")+ylim(0,0.02)

p27 <- ggplot() +
  geom_line(aes(x = sizes*100, y = SE_size[1,,1], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = SE_size[1,,1], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = SE_size[2,,1], colour = "Sandwich")) +
  geom_point(aes(x = sizes*100, y = SE_size[2,,1], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Sample size', 
       y = "Monte Carlo SE")+ylim(0,0.02)

figure7 <- ggarrange(p25 + rremove("ylab"),p26 + rremove("ylab"),p27 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure7 <- annotate_figure(figure7,top = text_grob("Full model",size = 14), left = text_grob("Monte Carlo SE",size = 10, rot =  90))

p28 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = SE[5,,1], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = SE[5,,1], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = SE[6,,1], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = SE[6,,1], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Confounding strength', 
       y = "Monte Carlo SE")+ylim(0,0.02)

p29 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = SE[7,,1], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = SE[7,,1], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = SE[8,,1], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = SE[8,,1], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Prevalence of treatment', 
       y = "Monte Carlo SE")+ylim(0,0.02)

p30 <- ggplot() +
  geom_line(aes(x = sizes*100, y = SE_size[3,,1], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = SE_size[3,,1], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = SE_size[4,,1], colour = "Sandwich")) +
  geom_point(aes(x = sizes*100, y = SE_size[4,,1], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  labs(x = 'Sizes', 
       y = "Monte Carlo SE")+ylim(0,0.02)

figure8 <- ggarrange(p28 + rremove("ylab"),p29 + rremove("ylab"),p30 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figure8 <- annotate_figure(figure8,top = text_grob("Simple model",size = 14), left = text_grob("Monte Carlo SE",size = 10, rot =  90))

ggarrange(heights = c(4,4.9), figure7, figure8, nrow = 2, ncol = 1)


