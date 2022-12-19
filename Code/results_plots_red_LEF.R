library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/MPhil-dissertation/Code")
library(ggplot2)
library(ggpubr)
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

treat_pos <- c(-1,-0.8,-0.5,-0.2,0,0.2,0.5,0.8,1)

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

mean_time_coefs <- array(,dim = c(4,9))
mean_time_treat <- array(,dim = c(4,9))
mean_time_size <- array(,dim = c(4,7))

mean_time_coefs_big <- array(,dim = c(4,9))
mean_time_treat_big <- array(,dim = c(4,9))

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
  
  mean_time_coefs[,i] <- rowMeans(time_coefs[,,i], na.rm = TRUE)
  mean_time_treat[,i] <- rowMeans(time_treat[,,i], na.rm = TRUE)
  
  mean_time_coefs_big[,i] <- rowMeans(time_coefs_big[,,i], na.rm = TRUE)
  mean_time_treat_big[,i] <- rowMeans(time_treat_big[,,i], na.rm = TRUE)
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
  
  mean_time_size[,i] <- rowMeans(time_size[,,i], na.rm = TRUE)
}



######### MEAN SD LENGTH PLOTS ###################
p1 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_lengths_coefs[1,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths_coefs[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths_coefs[2,], colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths_coefs[2,], colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths_coefs[3,], colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths_coefs[3,], colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths_coefs[4,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths_coefs[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Confounding strength\n(N = 1000, Treat. prev = 0)', 
       y = "Mean CI length")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p2 <- ggplot() +
  geom_line(aes(x = treat_pos, y = mean_lengths_treat[1,], colour = "Bootstrap")) +
  geom_point(aes(x = treat_pos, y = mean_lengths_treat[1,], colour = "Bootstrap")) +
  geom_line(aes(x = treat_pos, y = mean_lengths_treat[2,], colour = "LEF outcome")) +
  geom_point(aes(x = treat_pos, y = mean_lengths_treat[2,], colour = "LEF outcome")) +
  geom_line(aes(x = treat_pos, y = mean_lengths_treat[3,], colour = "LEF both")) +
  geom_point(aes(x = treat_pos, y = mean_lengths_treat[3,], colour = "LEF both")) +
  geom_line(aes(x = treat_pos, y = mean_lengths_treat[4,], colour = "Sandwich")) +
  geom_point(aes(x = treat_pos, y = mean_lengths_treat[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Treatment prevalence\n(N = 1000, Confounding = 0.5)', 
       y = "Mean CI length")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p3 <- ggplot() +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[1,], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = mean_lengths_size[1,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[2,], colour = "LEF outcome")) +
  geom_point(aes(x = sizes*100, y = mean_lengths_size[2,], colour = "LEF outcome")) +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[3,], colour = "LEF both")) +
  geom_point(aes(x = sizes*100, y = mean_lengths_size[3,], colour = "LEF both")) +
  geom_line(aes(x = sizes*100, y = mean_lengths_size[4,], colour = "Sandwich")) +
  geom_point(aes(x = sizes*100, y = mean_lengths_size[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Patient sample size N\n(Confounding  = 0.5, Treat. prev = 0)', 
       y = "Mean CI length")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figureA <- ggarrange(p1 + rremove("ylab"),p2 + rremove("ylab"),p3 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figureA <- annotate_figure(figureA,top = text_grob("Mean CI length",size = 14), left = text_grob("Mean CI length",size = 10, rot =  90))

p1_sd <- ggplot() +
  geom_line(aes(x = 1:9/10, y = sd_lengths_coefs[1,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths_coefs[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = sd_lengths_coefs[2,], colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths_coefs[2,], colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = sd_lengths_coefs[3,], colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths_coefs[3,], colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = sd_lengths_coefs[4,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths_coefs[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Confounding strength\n(N = 1000, Treat. prev = 0)', 
       y = "Standard deviation of CI length")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p2_sd <- ggplot() +
  geom_line(aes(x = treat_pos, y = sd_lengths_treat[1,], colour = "Bootstrap")) +
  geom_point(aes(x = treat_pos, y = sd_lengths_treat[1,], colour = "Bootstrap")) +
  geom_line(aes(x = treat_pos, y = sd_lengths_treat[2,], colour = "LEF outcome")) +
  geom_point(aes(x = treat_pos, y = sd_lengths_treat[2,], colour = "LEF outcome")) +
  geom_line(aes(x = treat_pos, y = sd_lengths_treat[3,], colour = "LEF both")) +
  geom_point(aes(x = treat_pos, y = sd_lengths_treat[3,], colour = "LEF both")) +
  geom_line(aes(x = treat_pos, y = sd_lengths_treat[4,], colour = "Sandwich")) +
  geom_point(aes(x = treat_pos, y = sd_lengths_treat[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Treatment prevalence\n(N = 1000, Confounding = 0.5)', 
       y = "Standard deviation of CI length")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p3_sd <- ggplot() +
  geom_line(aes(x = sizes*100, y = sd_lengths_size[1,], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = sd_lengths_size[1,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = sd_lengths_size[2,], colour = "LEF outcome")) +
  geom_point(aes(x = sizes*100, y = sd_lengths_size[2,], colour = "LEF outcome")) +
  geom_line(aes(x = sizes*100, y = sd_lengths_size[3,], colour = "LEF both")) +
  geom_point(aes(x = sizes*100, y = sd_lengths_size[3,], colour = "LEF both")) +
  geom_line(aes(x = sizes*100, y = sd_lengths_size[4,], colour = "Sandwich")) +
  geom_point(aes(x = sizes*100, y = sd_lengths_size[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Patient sample size N\n(Confounding  = 0.5, Treat. prev = 0)', 
       y = "Standard deviation of CI length")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figureA_sd <- ggarrange(p1_sd + rremove("ylab"),p2_sd + rremove("ylab"),p3_sd + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figureA_sd <- annotate_figure(figureA_sd,top = text_grob("Standard deviation of CI length",size = 14), left = text_grob("CI length SD",size = 10, rot =  90))
ggarrange(heights = c(4, 4.9),figureA, figureA_sd, nrow = 2, ncol = 1)

######### MEAN SD LENGTH PLOTS BIG #################
p1 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_lengths_coefs_big[1,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths_coefs_big[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths_coefs_big[2,], colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths_coefs_big[2,], colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths_coefs_big[3,], colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths_coefs_big[3,], colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = mean_lengths_coefs_big[4,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_lengths_coefs_big[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Confounding strength\n(N = 5000, Treat. prev = 0)', 
       y = "Mean CI length")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p2 <- ggplot() +
  geom_line(aes(x = treat_pos, y = mean_lengths_treat_big[1,], colour = "Bootstrap")) +
  geom_point(aes(x = treat_pos, y = mean_lengths_treat_big[1,], colour = "Bootstrap")) +
  geom_line(aes(x = treat_pos, y = mean_lengths_treat_big[2,], colour = "LEF outcome")) +
  geom_point(aes(x = treat_pos, y = mean_lengths_treat_big[2,], colour = "LEF outcome")) +
  geom_line(aes(x = treat_pos, y = mean_lengths_treat_big[3,], colour = "LEF both")) +
  geom_point(aes(x = treat_pos, y = mean_lengths_treat_big[3,], colour = "LEF both")) +
  geom_line(aes(x = treat_pos, y = mean_lengths_treat_big[4,], colour = "Sandwich")) +
  geom_point(aes(x = treat_pos, y = mean_lengths_treat_big[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Treatment prevalence\n(N = 5000, Confounding = 0.5)', 
       y = "Mean CI length")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figureA <- ggarrange(p1 + rremove("ylab"),p2 + rremove("ylab"), nrow = 1, ncol = 2, common.legend = T, legend = 'none')
figureA <- annotate_figure(figureA,top = text_grob("Mean CI length",size = 14), left = text_grob("Mean CI length",size = 10, rot =  90))

p1_sd <- ggplot() +
  geom_line(aes(x = 1:9/10, y = sd_lengths_coefs_big[1,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths_coefs_big[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = sd_lengths_coefs_big[2,], colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths_coefs_big[2,], colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = sd_lengths_coefs_big[3,], colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths_coefs_big[3,], colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = sd_lengths_coefs_big[4,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = sd_lengths_coefs_big[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Confounding strength\n(N = 5000, Treat. prev = 0)', 
       y = "Standard deviation of CI length")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p2_sd <- ggplot() +
  geom_line(aes(x = treat_pos, y = sd_lengths_treat_big[1,], colour = "Bootstrap")) +
  geom_point(aes(x = treat_pos, y = sd_lengths_treat_big[1,], colour = "Bootstrap")) +
  geom_line(aes(x = treat_pos, y = sd_lengths_treat_big[2,], colour = "LEF outcome")) +
  geom_point(aes(x = treat_pos, y = sd_lengths_treat_big[2,], colour = "LEF outcome")) +
  geom_line(aes(x = treat_pos, y = sd_lengths_treat_big[3,], colour = "LEF both")) +
  geom_point(aes(x = treat_pos, y = sd_lengths_treat_big[3,], colour = "LEF both")) +
  geom_line(aes(x = treat_pos, y = sd_lengths_treat_big[4,], colour = "Sandwich")) +
  geom_point(aes(x = treat_pos, y = sd_lengths_treat_big[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Treatment prevalence\n(N = 5000, Confounding = 0.5)', 
       y = "Standard deviation of CI length")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figureA_sd <- ggarrange(p1_sd + rremove("ylab"),p2_sd + rremove("ylab"), nrow = 1, ncol = 2, common.legend = T, legend = 'bottom')
figureA_sd <- annotate_figure(figureA_sd,top = text_grob("Standard deviation of CI length",size = 14), left = text_grob("CI length SD",size = 10, rot =  90))
ggarrange(heights = c(4, 4.9),figureA, figureA_sd, nrow = 2, ncol = 1)

################# COMPUTATION TIME PLOTS ##################
p1 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_time_coefs[1,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_time_coefs[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_time_coefs[2,], colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = mean_time_coefs[2,], colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = mean_time_coefs[3,], colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = mean_time_coefs[3,], colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = mean_time_coefs[4,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_time_coefs[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Confounding strength\n(N = 1000, Treat. prev. = 0)', 
       y = "Mean computation time") + theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p2 <- ggplot() +
  geom_line(aes(x = treat_pos, y = mean_time_treat[1,], colour = "Bootstrap")) +
  geom_point(aes(x = treat_pos, y = mean_time_treat[1,], colour = "Bootstrap")) +
  geom_line(aes(x = treat_pos, y = mean_time_treat[2,], colour = "LEF outcome")) +
  geom_point(aes(x = treat_pos, y = mean_time_treat[2,], colour = "LEF outcome")) +
  geom_line(aes(x = treat_pos, y = mean_time_treat[3,], colour = "LEF both")) +
  geom_point(aes(x = treat_pos, y = mean_time_treat[3,], colour = "LEF both")) +
  geom_line(aes(x = treat_pos, y = mean_time_treat[4,], colour = "Sandwich")) +
  geom_point(aes(x = treat_pos, y = mean_time_treat[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Treatment prevalence\n(N = 1000, Confounding = 0.5)', 
       y = "Mean computation time") + theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p3 <- ggplot() +
  geom_line(aes(x = sizes*100, y = mean_time_size[1,], colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = mean_time_size[1,], colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = mean_time_size[2,], colour = "LEF outcome")) +
  geom_point(aes(x = sizes*100, y = mean_time_size[2,], colour = "LEF outcome")) +
  geom_line(aes(x = sizes*100, y = mean_time_size[3,], colour = "LEF both")) +
  geom_point(aes(x = sizes*100, y = mean_time_size[3,], colour = "LEF both")) +
  geom_line(aes(x = sizes*100, y = mean_time_size[4,], colour = "Sandwich")) +
  geom_point(aes(x = sizes*100, y = mean_time_size[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Patient sample size N\n(Confounding  = 0.5, Treat. prev = 0)', 
       y = "Mean computation time") + theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figureT <- ggarrange(p1 + rremove("ylab"),p2 + rremove("ylab"),p3 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figureT <- annotate_figure(figureT,top = text_grob("Mean CI computation time",size = 14), left = text_grob("Mean computation time (sec.)",size = 10, rot =  90))
figureT

############ COMPUTATION TIME BIG #################
p1 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = mean_time_coefs_big[1,], colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = mean_time_coefs_big[1,], colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = mean_time_coefs_big[2,], colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = mean_time_coefs_big[2,], colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = mean_time_coefs_big[3,], colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = mean_time_coefs_big[3,], colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = mean_time_coefs_big[4,], colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = mean_time_coefs_big[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Confounding strength\n(N = 5000, Treat. prev. = 0)', 
       y = "Mean computation time")+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p2 <- ggplot() +
  geom_line(aes(x = treat_pos, y = mean_time_treat_big[1,], colour = "Bootstrap")) +
  geom_point(aes(x = treat_pos, y = mean_time_treat_big[1,], colour = "Bootstrap")) +
  geom_line(aes(x = treat_pos, y = mean_time_treat_big[2,], colour = "LEF outcome")) +
  geom_point(aes(x = treat_pos, y = mean_time_treat_big[2,], colour = "LEF outcome")) +
  geom_line(aes(x = treat_pos, y = mean_time_treat_big[3,], colour = "LEF both")) +
  geom_point(aes(x = treat_pos, y = mean_time_treat_big[3,], colour = "LEF both")) +
  geom_line(aes(x = treat_pos, y = mean_time_treat_big[4,], colour = "Sandwich")) +
  geom_point(aes(x = treat_pos, y = mean_time_treat_big[4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  labs(x = 'Treatment prevalence\n(N = 5000, Confounding = 0.5)', 
       y = "Mean computation time",)+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
figureT <- ggarrange(p1 + rremove("ylab"),p2 + rremove("ylab"), nrow = 1, ncol = 2, common.legend = T, legend = 'bottom')
figureT <- annotate_figure(figureT,top = text_grob("Mean CI computation time",size = 14), left = text_grob("Mean computation time (sec.)",size = 10, rot =  90))
figureT
############# COVERAGE #####################

coverage_ind_coefs <- array(0,dim = c(4,9,5))
coverage_ind_treat <- array(0,dim = c(4,9,5))
coverage_ind_size <- array(0,dim = c(4,7,5))

coverage_ind_coefs_big <- array(0,dim = c(4,9,5))
coverage_ind_treat_big <- array(0,dim = c(4,9,5))

success_coefs <- array(0,dim = c(4,9,5))
success_treat <- array(0,dim = c(4,9,5))
success_size <- array(0,dim = c(4,7,5))

success_coefs_big <- array(0,dim = c(4,9,5))
success_treat_big <- array(0,dim = c(4,9,5))

for (i in 1:1000){
  for (k in 1:5){
    for (j in 1:9){
      if (is.na(bootstrap_coefs[k,1,i,j]) == F){
        success_coefs[1,j,k] <- success_coefs[1,j,k] + 1
        if (all(bootstrap_coefs[k,1,i,j] <= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j]) 
            & all(bootstrap_coefs[k,2,i,j] >= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j])){
          coverage_ind_coefs[1,j,k] <- coverage_ind_coefs[1,j,k] + 1
        }
      }
      
      if (is.na(LEF_outcome_coefs[k,1,i,j]) == F){
        success_coefs[2,j,k] <- success_coefs[2,j,k] + 1
        if (all(LEF_outcome_coefs[k,1,i,j] <= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j]) 
            & all(LEF_outcome_coefs[k,2,i,j] >= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j])){
          coverage_ind_coefs[2,j,k] <- coverage_ind_coefs[2,j,k] + 1
        }
      }
      
      if (is.na(LEF_both_coefs[k,1,i,j]) == F){
        success_coefs[3,j,k] <- success_coefs[3,j,k] + 1
        if (all(LEF_both_coefs[k,1,i,j] <= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j]) 
            & all(LEF_both_coefs[k,2,i,j] >= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j])){
          coverage_ind_coefs[3,j,k] <- coverage_ind_coefs[3,j,k] + 1
        }
      }
      
      if (all(is.na(sandwich_coefs[k,1,i,j])) == F){
        success_coefs[4,j,k] <- success_coefs[4,j,k] + 1
        if (all(sandwich_coefs[k,1,i,j] <= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j]) 
            & all(sandwich_coefs[k,2,i,j] >= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j])){
          coverage_ind_coefs[4,j,k]<- coverage_ind_coefs[4,j,k]+ 1
        }
      }
      
      if (all(is.na(bootstrap_treat[k,1,i,j])) == F){
        success_treat[1,j,k]<- success_treat[1,j,k]+ 1
        if (all(bootstrap_treat[k,1,i,j] <= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j]) 
            & all(bootstrap_treat[k,2,i,j] >= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j])){
          coverage_ind_treat[1,j,k]<- coverage_ind_treat[1,j,k]+ 1
        }
      }
      if (is.na(LEF_outcome_treat[k,1,i,j]) == F){
        success_treat[2,j,k] <- success_treat[2,j,k] + 1
        if (all(LEF_outcome_treat[k,1,i,j] <= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j]) 
            & all(LEF_outcome_treat[k,2,i,j] >= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j])){
          coverage_ind_treat[2,j,k] <- coverage_ind_treat[2,j,k] + 1
        }
      }
      
      if (is.na(LEF_both_treat[k,1,i,j]) == F){
        success_treat[3,j,k] <- success_treat[3,j,k] + 1
        if (all(LEF_both_treat[k,1,i,j] <= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j]) 
            & all(LEF_both_treat[k,2,i,j] >= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j])){
          coverage_ind_treat[3,j,k] <- coverage_ind_treat[3,j,k] + 1
        }
      }
      if (all(is.na(sandwich_treat[k,1,i,j])) == F){
        success_treat[4,j,k]<- success_treat[4,j,k]+ 1
        if (all(sandwich_treat[k,1,i,j] <= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j]) 
            & all(sandwich_treat[k,2,i,j] >= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j])){
          coverage_ind_treat[4,j,k]<- coverage_ind_treat[4,j,k]+ 1
        }
      }
    ### big
      if (is.na(bootstrap_coefs_big[k,1,i,j]) == F){
        success_coefs_big[1,j,k] <- success_coefs_big[1,j,k] + 1
        if (all(bootstrap_coefs_big[k,1,i,j] <= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j]) 
            & all(bootstrap_coefs_big[k,2,i,j] >= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j])){
          coverage_ind_coefs_big[1,j,k] <- coverage_ind_coefs_big[1,j,k] + 1
        }
      }
      
      if (is.na(LEF_outcome_coefs_big[k,1,i,j]) == F){
        success_coefs_big[2,j,k] <- success_coefs_big[2,j,k] + 1
        if (all(LEF_outcome_coefs_big[k,1,i,j] <= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j]) 
            & all(LEF_outcome_coefs_big[k,2,i,j] >= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j])){
          coverage_ind_coefs_big[2,j,k] <- coverage_ind_coefs_big[2,j,k] + 1
        }
      }
      
      if (is.na(LEF_both_coefs_big[k,1,i,j]) == F){
        success_coefs_big[3,j,k] <- success_coefs_big[3,j,k] + 1
        if (all(LEF_both_coefs_big[k,1,i,j] <= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j]) 
            & all(LEF_both_coefs_big[k,2,i,j] >= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j])){
          coverage_ind_coefs_big[3,j,k] <- coverage_ind_coefs_big[3,j,k] + 1
        }
      }
      
      if (all(is.na(sandwich_coefs_big[k,1,i,j])) == F){
        success_coefs_big[4,j,k] <- success_coefs_big[4,j,k] + 1
        if (all(sandwich_coefs_big[k,1,i,j] <= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j]) 
            & all(sandwich_coefs_big[k,2,i,j] >= true_value_conf_red[k,1,j] - true_value_conf_red[k,2,j])){
          coverage_ind_coefs_big[4,j,k]<- coverage_ind_coefs_big[4,j,k]+ 1
        }
      }
      
      if (all(is.na(bootstrap_treat_big[k,1,i,j])) == F){
        success_treat_big[1,j,k]<- success_treat_big[1,j,k]+ 1
        if (all(bootstrap_treat_big[k,1,i,j] <= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j]) 
            & all(bootstrap_treat_big[k,2,i,j] >= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j])){
          coverage_ind_treat_big[1,j,k]<- coverage_ind_treat_big[1,j,k]+ 1
        }
      }
      if (is.na(LEF_outcome_treat_big[k,1,i,j]) == F){
        success_treat_big[2,j,k] <- success_treat_big[2,j,k] + 1
        if (all(LEF_outcome_treat_big[k,1,i,j] <= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j]) 
            & all(LEF_outcome_treat_big[k,2,i,j] >= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j])){
          coverage_ind_treat_big[2,j,k] <- coverage_ind_treat_big[2,j,k] + 1
        }
      }
      
      if (is.na(LEF_both_treat_big[k,1,i,j]) == F){
        success_treat_big[3,j,k] <- success_treat_big[3,j,k] + 1
        if (all(LEF_both_treat_big[k,1,i,j] <= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j]) 
            & all(LEF_both_treat_big[k,2,i,j] >= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j])){
          coverage_ind_treat_big[3,j,k] <- coverage_ind_treat_big[3,j,k] + 1
        }
      }
      if (all(is.na(sandwich_treat_big[k,1,i,j])) == F){
        success_treat_big[4,j,k]<- success_treat_big[4,j,k]+ 1
        if (all(sandwich_treat_big[k,1,i,j] <= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j]) 
            & all(sandwich_treat_big[k,2,i,j] >= true_value_treat_red[k,1,j] - true_value_treat_red[k,2,j])){
          coverage_ind_treat_big[4,j,k]<- coverage_ind_treat_big[4,j,k]+ 1
        }
      }
    } 
    for (j in 1:7){
      if (all(is.na(bootstrap_size[k,1,i,j])) == F){
        success_size[1,j,k]<- success_size[1,j,k]+ 1
        if (all(bootstrap_size[k,1,i,j] <= true_value_conf_red[k,1,5] - true_value_conf_red[k,2,5]) 
            & all(bootstrap_size[k,2,i,j] >= true_value_conf_red[k,1,5] - true_value_conf_red[k,2,5])){
          coverage_ind_size[1,j,k]<- coverage_ind_size[1,j,k]+ 1
        }
      }
      if (is.na(LEF_outcome_size[k,1,i,j]) == F){
        success_size[2,j,k] <- success_size[2,j,k] + 1
        if (all(LEF_outcome_size[k,1,i,j] <= true_value_conf_red[k,1,5] - true_value_conf_red[k,2,5]) 
            & all(LEF_outcome_size[k,2,i,j] >= true_value_conf_red[k,1,5] - true_value_conf_red[k,2,5])){
          coverage_ind_size[2,j,k] <- coverage_ind_size[2,j,k] + 1
        }
      }
      
      if (is.na(LEF_both_size[k,1,i,j]) == F){
        success_size[3,j,k] <- success_size[3,j,k] + 1
        if (all(LEF_both_size[k,1,i,j] <= true_value_conf_red[k,1,5] - true_value_conf_red[k,2,5]) 
            & all(LEF_both_size[k,2,i,j] >= true_value_conf_red[k,1,5] - true_value_conf_red[k,2,5])){
          coverage_ind_size[3,j,k] <- coverage_ind_size[3,j,k] + 1
        }
      }
      if (all(is.na(sandwich_size[k,1,i,j])) == F){
        success_size[4,j,k]<- success_size[4,j,k]+ 1
        if (all(sandwich_size[k,1,i,j] <= true_value_conf_red[k,1,5] - true_value_conf_red[k,2,5]) 
            & all(sandwich_size[k,2,i,j] >= true_value_conf_red[k,1,5] - true_value_conf_red[k,2,5])){
          coverage_ind_size[4,j,k]<- coverage_ind_size[4,j,k]+ 1
        }
      }
    }
  }
}

coverage_ind_coefs <- coverage_ind_coefs/success_coefs
coverage_ind_treat <- coverage_ind_treat/success_treat
coverage_ind_size <- coverage_ind_size/success_size

coverage_ind_coefs_big <- coverage_ind_coefs_big/success_coefs_big
coverage_ind_treat_big <- coverage_ind_treat_big/success_treat_big

############# COVERAGE PLOTS ##################
p7 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[1,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[1,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[2,1,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[2,1,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[3,1,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[3,1,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[4,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[4,1,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.1")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
p8 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[1,5,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[1,5,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[2,5,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[2,5,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[3,5,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[3,5,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[4,5,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[4,5,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.5")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p9 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[1,9,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[1,9,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[2,9,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[2,9,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[3,9,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[3,9,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs[4,9,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs[4,9,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.9")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figure1 <- ggarrange(p7+ rremove("ylab"),p8+ rremove("ylab"),p9+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure1 <- annotate_figure(figure1,top = text_grob("N = 1000, Treatment prev. = 0",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))

p10 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[1,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[1,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[2,1,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[2,1,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[3,1,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[3,1,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[4,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[4,1,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Treatment prevalence = -1")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ 
  theme(aspect.ratio = 1, axis.title = element_text(size = 10))
p11 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[1,5,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[1,5,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[2,5,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[2,5,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[3,5,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[3,5,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[4,5,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[4,5,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Treatment prevalence = 0")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p12 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[1,9,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[1,9,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[2,9,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[2,9,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[3,9,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[3,9,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat[4,9,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat[4,9,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Treatment prevalence = 1")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figure2 <- ggarrange(p10 + rremove("ylab"),p11 + rremove("ylab"),p12 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T,legend = 'bottom')
figure2 <- annotate_figure(figure2,top = text_grob("N = 1000, Confounding strength = 0.5",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))
ggarrange(heights = c(4, 4.9),figure1, figure2, nrow = 2, ncol = 1)

############# COVERAGE BIG ###################
p7 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[1,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[1,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[2,1,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[2,1,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[3,1,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[3,1,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[4,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[4,1,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.1")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
p8 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[1,5,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[1,5,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[2,5,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[2,5,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[3,5,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[3,5,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[4,5,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[4,5,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.5")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p9 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[1,9,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[1,9,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[2,9,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[2,9,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[3,9,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[3,9,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_coefs_big[4,9,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_coefs_big[4,9,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Confounding = 0.9")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figure1 <- ggarrange(p7+ rremove("ylab"),p8+ rremove("ylab"),p9+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure1 <- annotate_figure(figure1,top = text_grob("N = 5000, Treatment prev. = 0",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))

p10 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[1,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[1,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[2,1,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[2,1,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[3,1,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[3,1,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[4,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[4,1,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Treatment prevalence = -1")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
p11 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[1,5,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[1,5,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[2,5,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[2,5,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[3,5,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[3,5,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[4,5,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[4,5,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Treatment prevalence = 0")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p12 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[1,9,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[1,9,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[2,9,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[2,9,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[3,9,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[3,9,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_treat_big[4,9,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_treat_big[4,9,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "Treatment prevalence = 1")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figure2 <- ggarrange(p10 + rremove("ylab"),p11 + rremove("ylab"),p12 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T,legend = 'bottom')
figure2 <- annotate_figure(figure2,top = text_grob("N = 5000, Confounding strength = 0.5",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))
ggarrange(heights = c(4, 4.9),figure1, figure2, nrow = 2, ncol = 1)

###################### COVERAGE SIZE ###########
p19 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_size[1,1,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[1,1,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[2,1,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[2,1,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[3,1,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[3,1,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[4,1,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[4,1,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple"))+
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "N = 200")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
p20 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_size[1,4,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[1,4,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[2,4,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[2,4,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[3,4,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[3,4,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[4,4,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[4,4,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple"))+
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "N = 1000")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
p21 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_size[1,6,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[1,6,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[2,6,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[2,6,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[3,6,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[3,6,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[4,6,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[4,6,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple"))+
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "N = 5000")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p22 <- ggplot() +
  geom_line(aes(x = 0:4, y = coverage_ind_size[1,7,], colour = "Bootstrap")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[1,7,], colour = "Bootstrap")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[2,7,], colour = "LEF outcome")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[2,7,], colour = "LEF outcome")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[3,7,], colour = "LEF both")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[3,7,], colour = "LEF both")) +
  geom_line(aes(x = 0:4, y = coverage_ind_size[4,7,], colour = "Sandwich")) +
  geom_point(aes(x = 0:4, y = coverage_ind_size[4,7,], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple"))+
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical coverage rate",
       title = "N = 10000")+ ylim(0.6,1) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figure5 <- ggarrange(p19 + rremove("ylab"),p20 + rremove("ylab"),
                     p21 + rremove("ylab"),p22 + rremove("ylab"), nrow = 1, ncol = 4, common.legend = T, legend = 'bottom')
figure5 <- annotate_figure(figure5,top = text_grob("Confounding = 0.5, Treat. prev. = 0",size = 14),
                           left = text_grob("Empirical CI coverage",size = 10, rot =  90))
figure5

################ MC SE PLOTS #######################
SE_coefs <- sqrt((coverage_ind_coefs*(1-coverage_ind_coefs))/1000)
SE_treat <- sqrt((coverage_ind_treat*(1-coverage_ind_treat))/1000)

SE_coefs_big <- sqrt((coverage_ind_coefs_big*(1-coverage_ind_coefs_big))/1000)
SE_treat_big <- sqrt((coverage_ind_treat_big*(1-coverage_ind_treat_big))/1000)

SE_size <- sqrt((coverage_ind_size*(1-coverage_ind_size))/1000)

p25 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_coefs[1,,]), colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_coefs[1,,]), colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_coefs[2,,]), colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_coefs[2,,]), colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_coefs[3,,]), colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_coefs[3,,]), colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_coefs[4,,]), colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_coefs[4,,]), colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple"))+
  labs(x = 'Confounding strength\n(N = 1000, Treat. prev = 0)', 
       y = "Monte Carlo SE") +ylim(0,0.02)+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p26 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_treat[1,,]), colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_treat[1,,]), colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_treat[2,,]), colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_treat[2,,]), colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_treat[3,,]), colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_treat[3,,]), colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_treat[4,,]), colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_treat[4,,]), colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple"))+
  labs(x = 'Treatment prevalence\n(N = 1000, Confounding = 0.5)', 
       y = "Monte Carlo SE") +ylim(0,0.02)+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p27 <- ggplot() +
  geom_line(aes(x = sizes*100, y = rowMeans(SE_size[1,,]), colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = rowMeans(SE_size[1,,]), colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = rowMeans(SE_size[2,,]), colour = "LEF outcome")) +
  geom_point(aes(x = sizes*100, y = rowMeans(SE_size[2,,]), colour = "LEF outcome")) +
  geom_line(aes(x = sizes*100, y = rowMeans(SE_size[3,,]), colour = "LEF both")) +
  geom_point(aes(x = sizes*100, y = rowMeans(SE_size[3,,]), colour = "LEF both")) +
  geom_line(aes(x = sizes*100, y = rowMeans(SE_size[4,,]), colour = "Sandwich")) +
  geom_point(aes(x = sizes*100, y = rowMeans(SE_size[4,,]), colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple"))+
  labs(x = 'Patient sample size N\n(Confounding  = 0.5, Treat. prev = 0)', 
       y = "Monte Carlo SE") +ylim(0,0.02)+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figure7 <- ggarrange(p25 + rremove("ylab"),p26 + rremove("ylab"),p27 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figure7 <- annotate_figure(figure7,top = text_grob("Monte Carlo Simulation Standard Error",size = 14), left = text_grob("Monte Carlo SE",size = 10, rot =  90))
figure7
############### MC ERROR BIG ######################
p28 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_coefs_big[1,,]), colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_coefs_big[1,,]), colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_coefs_big[2,,]), colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_coefs_big[2,,]), colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_coefs_big[3,,]), colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_coefs_big[3,,]), colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_coefs_big[4,,]), colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_coefs_big[4,,]), colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple"))+
  labs(x = 'Confounding strength\n(N = 5000, Treat. prev = 0)', 
       y = "Monte Carlo SE") +ylim(0,0.02)+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p29 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_treat_big[1,,]), colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_treat_big[1,,]), colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_treat_big[2,,]), colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_treat_big[2,,]), colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_treat_big[3,,]), colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_treat_big[3,,]), colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = rowMeans(SE_treat_big[4,,]), colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = rowMeans(SE_treat_big[4,,]), colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple"))+
  labs(x = 'Treatment prevalence\n(N = 1000, Confounding = 0.5)', 
       y = "Monte Carlo SE") +ylim(0,0.02)+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))


figure8 <- ggarrange(p28 + rremove("ylab"),p29 + rremove("ylab"), nrow = 1, ncol = 2, common.legend = T, legend = 'bottom')
figure8 <- annotate_figure(figure8,top = text_grob("Monte Carlo Simulation Standard Error",size = 14), left = text_grob("Monte Carlo SE",size = 10, rot =  90))
figure8

########## METHOD FAILURE PLOTS ##############
p7 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = (1000-success_coefs[1,,1])/1000, colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = (1000-success_coefs[1,,1])/1000, colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = (1000-success_coefs[2,,1])/1000, colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = (1000-success_coefs[2,,1])/1000, colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = (1000-success_coefs[3,,1])/1000, colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = (1000-success_coefs[3,,1])/1000, colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = (1000-success_coefs[4,,1])/1000, colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = (1000-success_coefs[4,,1])/1000, colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Confounding strength\n(N = 1000, Treat. prev. = 0)', 
       y = "Method failure rate")+ ylim(0,0.4) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p10 <- ggplot() +
  geom_line(aes(x = treat_pos, y = (1000-success_treat[1,,1])/1000, colour = "Bootstrap")) +
  geom_point(aes(x = treat_pos, y = (1000-success_treat[1,,1])/1000, colour = "Bootstrap")) +
  geom_line(aes(x = treat_pos, y = (1000-success_treat[2,,1])/1000, colour = "LEF outcome")) +
  geom_point(aes(x = treat_pos, y = (1000-success_treat[2,,1])/1000, colour = "LEF outcome")) +
  geom_line(aes(x = treat_pos, y = (1000-success_treat[3,,1])/1000, colour = "LEF both")) +
  geom_point(aes(x = treat_pos, y = (1000-success_treat[3,,1])/1000, colour = "LEF both")) +
  geom_line(aes(x = treat_pos, y = (1000-success_treat[4,,1])/1000, colour = "Sandwich")) +
  geom_point(aes(x = treat_pos, y = (1000-success_treat[4,,1])/1000, colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Treatment prevalence\n(N = 1000, Confounding = 0.5)', 
       y = "Method failure rate")+ ylim(0,0.4) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p12 <- ggplot() +
  geom_line(aes(x = sizes*100, y = (1000-success_size[1,,1])/1000, colour = "Bootstrap")) +
  geom_point(aes(x = sizes*100, y = (1000-success_size[1,,1])/1000, colour = "Bootstrap")) +
  geom_line(aes(x = sizes*100, y = (1000-success_size[2,,1])/1000, colour = "LEF outcome")) +
  geom_point(aes(x = sizes*100, y = (1000-success_size[2,,1])/1000, colour = "LEF outcome")) +
  geom_line(aes(x = sizes*100, y = (1000-success_size[3,,1])/1000, colour = "LEF both")) +
  geom_point(aes(x = sizes*100, y = (1000-success_size[3,,1])/1000, colour = "LEF both")) +
  geom_line(aes(x = sizes*100, y = (1000-success_size[4,,1])/1000, colour = "Sandwich")) +
  geom_point(aes(x = sizes*100, y = (1000-success_size[4,,1])/1000, colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Patient sample size N\n(Confounding = 0.5, Treat. prev. = 0)', 
       y = "Method failure rate")+ ylim(0,0.4) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

figure2 <- ggarrange(p7 + rremove("ylab"),p10 + rremove("ylab"),p12 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T,legend = 'bottom')
figure2 <- annotate_figure(figure2,top = text_grob("CI construction failure rate",size = 14),
                           left = text_grob("Failure rate",size = 10, rot =  90))
figure2

############ METHOD FAILURE BIG #####################
p7 <- ggplot() +
  geom_line(aes(x = 1:9/10, y = (1000-success_coefs_big[1,,1])/1000, colour = "Bootstrap")) +
  geom_point(aes(x = 1:9/10, y = (1000-success_coefs_big[1,,1])/1000, colour = "Bootstrap")) +
  geom_line(aes(x = 1:9/10, y = (1000-success_coefs_big[2,,1])/1000, colour = "LEF outcome")) +
  geom_point(aes(x = 1:9/10, y = (1000-success_coefs_big[2,,1])/1000, colour = "LEF outcome")) +
  geom_line(aes(x = 1:9/10, y = (1000-success_coefs_big[3,,1])/1000, colour = "LEF both")) +
  geom_point(aes(x = 1:9/10, y = (1000-success_coefs_big[3,,1])/1000, colour = "LEF both")) +
  geom_line(aes(x = 1:9/10, y = (1000-success_coefs_big[4,,1])/1000, colour = "Sandwich")) +
  geom_point(aes(x = 1:9/10, y = (1000-success_coefs_big[4,,1])/1000, colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Confounding strength\n(N = 5000, Treat. prev. = 0)', 
       y = "Method failure rate")+ ylim(0,0.4) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))

p10 <- ggplot() +
  geom_line(aes(x = treat_pos, y = (1000-success_treat_big[1,,1])/1000, colour = "Bootstrap")) +
  geom_point(aes(x = treat_pos, y = (1000-success_treat_big[1,,1])/1000, colour = "Bootstrap")) +
  geom_line(aes(x = treat_pos, y = (1000-success_treat_big[2,,1])/1000, colour = "LEF outcome")) +
  geom_point(aes(x = treat_pos, y = (1000-success_treat_big[2,,1])/1000, colour = "LEF outcome")) +
  geom_line(aes(x = treat_pos, y = (1000-success_treat_big[3,,1])/1000, colour = "LEF both")) +
  geom_point(aes(x = treat_pos, y = (1000-success_treat_big[3,,1])/1000, colour = "LEF both")) +
  geom_line(aes(x = treat_pos, y = (1000-success_treat_big[4,,1])/1000, colour = "Sandwich")) +
  geom_point(aes(x = treat_pos, y = (1000-success_treat_big[4,,1])/1000, colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue", 
                                                  "LEF outcome" = "green", "LEF both" = "purple")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(x = 'Treatment prevalence\n(N = 5000, Confounding = 0.5)', 
       y = "Method failure rate")+ ylim(0,0.4) + theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))


figure2 <- ggarrange(p7 + rremove("ylab"),p10 + rremove("ylab"), nrow = 1, ncol = 2, common.legend = T,legend = 'bottom')
figure2 <- annotate_figure(figure2,top = text_grob("CI construction failure rate",size = 14),
                           left = text_grob("Failure rate",size = 10, rot =  90))
figure2

