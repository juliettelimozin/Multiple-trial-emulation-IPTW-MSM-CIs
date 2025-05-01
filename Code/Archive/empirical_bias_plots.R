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

coefs <- array(, dim = c(10,1000,9))
treat <- array(, dim = c(10,1000,9))
size <- array(, dim = c(10,1000,6))

for (i in 1:9){
  load(paste0("HPC output/CI_bootstrap_coefs_PP_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_coefs_PP_", i, ".rda"))
  load(paste0("HPC output/CI_bootstrap_treat_PP_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_treat_PP_", i, ".rda"))
  load(paste0("HPC output/coefs_PP_",i,".rda"))
  load(paste0("HPC output/treat_PP_",i,".rda"))
  bootstrap_coefs[,,,i] <- CI_bootstrap_coefs_PP
  sandwich_coefs[,,,i] <- CI_sandwich_coefs_PP
  bootstrap_treat[,,,i] <- CI_bootstrap_treat_PP
  sandwich_treat[,,,i] <- CI_sandwich_treat_PP
  coefs[,,i] <- coefs_PP
  treat[,,i] <- treat_PP
}

sizes <- c(2,5,8,10,25,50)
for (i in 1:6){
  load(paste0("HPC output/CI_bootstrap_size_PP_", sizes[i], ".rda"))
  load(paste0("HPC output/CI_sandwich_size_PP_", sizes[i], ".rda"))
  load(paste0("HPC output/size_PP_",sizes[i],".rda"))
  bootstrap_size[,,,i] <- CI_bootstrap_size_PP
  sandwich_size[,,,i] <- CI_sandwich_size_PP
  size[,,i] <- size_PP
}

bootstrap_coefs_noint <- array(,dim = c(10,2,1000,9))
sandwich_coefs_noint <- array(,dim = c(10,2,1000,9))
bootstrap_treat_noint <- array(,dim = c(10,2,1000,9))
sandwich_treat_noint <- array(,dim = c(10,2,1000,9))
bootstrap_size_noint <- array(,dim = c(10,2,1000,6))
sandwich_size_noint <- array(,dim = c(10,2,1000,6))
coefs_noint <- array(, dim = c(10,1000,9))
treat_noint <- array(, dim = c(10,1000,9))
size_noint <- array(, dim = c(10,1000,6))

for (i in 1:9){
  load(paste0("HPC output/CI_bootstrap_coefs_PP_noint_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_coefs_PP_noint_", i, ".rda"))
  load(paste0("HPC output/CI_bootstrap_treat_PP_noint_", i, ".rda"))
  load(paste0("HPC output/CI_sandwich_treat_PP_noint_", i, ".rda"))
  load(paste0("HPC output/coefs_PP_noint_", i, ".rda"))
  load(paste0("HPC output/treat_PP_noint_", i, ".rda"))
  bootstrap_coefs_noint[,,,i] <- CI_bootstrap_coefs_PP
  sandwich_coefs_noint[,,,i] <- CI_sandwich_coefs_PP
  bootstrap_treat_noint[,,,i] <- CI_bootstrap_treat_PP
  sandwich_treat_noint[,,,i] <- CI_sandwich_treat_PP
  coefs_noint[,,i] <- coefs_PP_noint
  treat_noint[,,i] <- treat_PP_noint
}

for (i in 1:6){
  load(paste0("HPC output/CI_bootstrap_size_PP_noint_", sizes[i], ".rda"))
  load(paste0("HPC output/CI_sandwich_size_PP_noint_", sizes[i], ".rda"))
  load(paste0("HPC output/size_PP_noint_", sizes[i], ".rda"))
  bootstrap_size_noint[,,,i] <- CI_bootstrap_size_PP
  sandwich_size_noint[,,,i] <- CI_sandwich_size_PP
  size_noint[,,i] <- size_PP_noint
}

bias_point_coefs <- array(,dim = c(10,9))
bias_point_treat <- array(,dim = c(10,9))
bias_point_size <- array(,dim = c(10,6))

bias_point_coefs_noint <- 1*bias_point_coefs
bias_point_treat_noint <- 1*bias_point_coefs
bias_point_size_noint <- 1*bias_point_size

for (i in 1:9){
  bias_point_coefs[,i] <- rowMeans(coefs[,,i], na.rm = TRUE) - (true_value_conf[,1,i] - true_value_conf[,2,i])
  bias_point_treat[,i] <- rowMeans(treat[,,i], na.rm = TRUE) - (true_value_treat[,1,i] - true_value_treat[,2,i])
  bias_point_coefs_noint[,i] <- rowMeans(coefs_noint[,,i], na.rm = TRUE) - (true_value_conf[,1,i] - true_value_conf[,2,i])
  bias_point_treat_noint[,i] <- rowMeans(treat_noint[,,i], na.rm = TRUE) - (true_value_treat[,1,i] - true_value_treat[,2,i])
}

for (i in 1:6){
  bias_point_size[,i] <- rowMeans(size[,,i], na.rm = TRUE) - (true_value_conf[,1,5] - true_value_conf[,2,5])
  bias_point_size_noint[,i] <- rowMeans(size_noint[,,i], na.rm = TRUE) - (true_value_conf[,1,5] - true_value_conf[,2,5])
}

empirical_bias_mid <- array(, dim = c(10,8,9))
for (i in 1:9){
  empirical_bias_mid[,1,i] <- rowMeans((bootstrap_coefs[,2,,i]+ bootstrap_coefs[,1,,i])/2, na.rm = TRUE) - (true_value_conf[,1,i] - true_value_conf[,2,i])
  empirical_bias_mid[,2,i] <- rowMeans((sandwich_coefs[,2,,i]+ sandwich_coefs[,1,,i])/2, na.rm = TRUE) - (true_value_conf[,1,i] - true_value_conf[,2,i])
  empirical_bias_mid[,3,i] <- rowMeans((bootstrap_treat[,2,,i]+ bootstrap_treat[,1,,i])/2, na.rm = TRUE) - (true_value_treat[,1,i] - true_value_treat[,2,i])
  empirical_bias_mid[,4,i] <- rowMeans((sandwich_treat[,2,,i]+ sandwich_treat[,1,,i])/2, na.rm = TRUE) - (true_value_treat[,1,i] - true_value_treat[,2,i])
  empirical_bias_mid[,5,i] <- rowMeans((bootstrap_coefs_noint[,2,,i]+ bootstrap_coefs_noint[,1,,i])/2, na.rm = TRUE) - (true_value_conf[,1,i] - true_value_conf[,2,i])
  empirical_bias_mid[,6,i] <- rowMeans((sandwich_coefs_noint[,2,,i]+ sandwich_coefs_noint[,1,,i])/2, na.rm = TRUE) - (true_value_conf[,1,i] - true_value_conf[,2,i])
  empirical_bias_mid[,7,i] <- rowMeans((bootstrap_treat_noint[,2,,i]+ bootstrap_treat_noint[,1,,i])/2, na.rm = TRUE) - (true_value_treat[,1,i] - true_value_treat[,2,i])
  empirical_bias_mid[,8,i] <- rowMeans((sandwich_treat_noint[,2,,i]+ sandwich_treat_noint[,1,,i])/2, na.rm = TRUE) - (true_value_treat[,1,i] - true_value_treat[,2,i])
}

empirical_bias_mid_size <- array(, dim = c(10,4,6))
for (i in 1:6){
  empirical_bias_mid_size[,1,i] <- rowMeans((bootstrap_size[,2,,i]+ bootstrap_size[,1,,i])/2, na.rm = TRUE) - (true_value_conf[,1,5] - true_value_conf[,2,5])
  empirical_bias_mid_size[,2,i] <- rowMeans((sandwich_size[,2,,i]+ sandwich_size[,1,,i])/2, na.rm = TRUE) - (true_value_conf[,1,5] - true_value_conf[,2,5])
  empirical_bias_mid_size[,3,i] <- rowMeans((bootstrap_size_noint[,2,,i]+ bootstrap_size_noint[,1,,i])/2, na.rm = TRUE) - (true_value_conf[,1,5] - true_value_conf[,2,5])
  empirical_bias_mid_size[,4,i] <- rowMeans((sandwich_size_noint[,2,,i]+ sandwich_size_noint[,1,,i])/2, na.rm = TRUE) - (true_value_conf[,1,5] - true_value_conf[,2,5])
}


p13 <- ggplot()+
  geom_line(aes(x = 0:9, y = bias_point_coefs[,1], colour = "Full model")) +
  geom_point(aes(x = 0:9, y = bias_point_coefs[,1], colour = "Full model")) +
  geom_line(aes(x = 0:9, y = bias_point_coefs_noint[,1], colour = "Simple model")) +
  geom_point(aes(x = 0:9, y = bias_point_coefs_noint[,1], colour = "Simple model")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Model type", values = c("Full model"= "green", "Simple model" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Confounding = 0.1")+ 
  ylim(-0.02,0.02)+ 
  theme(plot.title = element_text(size=10))
p14 <- ggplot()+
  geom_line(aes(x = 0:9, y = bias_point_coefs[,5], colour = "Full model")) +
  geom_point(aes(x = 0:9, y = bias_point_coefs[,5], colour = "Full model")) +
  geom_line(aes(x = 0:9, y = bias_point_coefs_noint[,5], colour = "Simple model")) +
  geom_point(aes(x = 0:9, y = bias_point_coefs_noint[,5], colour = "Simple model")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Model type", values = c("Full model"= "green", "Simple model" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Confounding = 0.5")+ 
  ylim(-0.02,0.02)+ 
  theme(plot.title = element_text(size=10))
p15 <- ggplot()+
  geom_line(aes(x = 0:9, y = bias_point_coefs[,9], colour = "Full model")) +
  geom_point(aes(x = 0:9, y = bias_point_coefs[,9], colour = "Full model")) +
  geom_line(aes(x = 0:9, y = bias_point_coefs_noint[,9], colour = "Simple model")) +
  geom_point(aes(x = 0:9, y = bias_point_coefs_noint[,9], colour = "Simple model")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Model type", values = c("Full model"= "green", "Simple model" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Confounding = 0.9")+ 
  ylim(-0.02,0.02)+ 
  theme(plot.title = element_text(size=10))

figure5 <- ggarrange(p13+ rremove("ylab"),p14+ rremove("ylab"),p15+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure5 <- annotate_figure(figure5,
                           left = text_grob("Empirical bias",size = 10, rot =  90))

p16 <- ggplot()+
  geom_line(aes(x = 0:9, y = bias_point_treat[,1], colour = "Full model")) +
  geom_point(aes(x = 0:9, y = bias_point_treat[,1], colour = "Full model")) +
  geom_line(aes(x = 0:9, y = bias_point_treat_noint[,1], colour = "Simple model")) +
  geom_point(aes(x = 0:9, y = bias_point_treat_noint[,1], colour = "Simple model")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Model type", values = c("Full model"= "green", "Simple model" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Treatment prevalence = 0.1")+ 
  ylim(-0.02,0.02)+ 
  theme(plot.title = element_text(size=10))
p17 <- ggplot()+
  geom_line(aes(x = 0:9, y = bias_point_treat[,5], colour = "Full model")) +
  geom_point(aes(x = 0:9, y = bias_point_treat[,5], colour = "Full model")) +
  geom_line(aes(x = 0:9, y = bias_point_treat_noint[,5], colour = "Simple model")) +
  geom_point(aes(x = 0:9, y = bias_point_treat_noint[,5], colour = "Simple model")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Model type", values = c("Full model"= "green", "Simple model" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Treatment prevalence = 0.5")+ 
  ylim(-0.02,0.02)+ 
  theme(plot.title = element_text(size=10))
p18 <- ggplot()+
  geom_line(aes(x = 0:9, y = bias_point_treat[,9], colour = "Full model")) +
  geom_point(aes(x = 0:9, y = bias_point_treat[,9], colour = "Full model")) +
  geom_line(aes(x = 0:9, y = bias_point_treat_noint[,9], colour = "Simple model")) +
  geom_point(aes(x = 0:9, y = bias_point_treat_noint[,9], colour = "Simple model")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Model type", values = c("Full model"= "green", "Simple model" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Treatment prevalence = 0.9")+ 
  ylim(-0.02,0.02)+ 
  theme(plot.title = element_text(size=10))

figure6 <- ggarrange(p16+ rremove("ylab"),p17+ rremove("ylab"),p18+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure6 <- annotate_figure(figure6,
                           left = text_grob("Empirical bias",size = 10, rot =  90))

p19 <- ggplot()+
  geom_line(aes(x = 0:9, y = bias_point_size[,1], colour = "Full model")) +
  geom_point(aes(x = 0:9, y = bias_point_size[,1], colour = "Full model")) +
  geom_line(aes(x = 0:9, y = bias_point_size_noint[,1], colour = "Simple model")) +
  geom_point(aes(x = 0:9, y = bias_point_size_noint[,1], colour = "Simple model")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Model type", values = c("Full model"= "green", "Simple model" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Sample size = 200")+ 
  ylim(-0.02,0.02)+ 
  theme(plot.title = element_text(size=10))
p20 <- ggplot()+
  geom_line(aes(x = 0:9, y = bias_point_size[,4], colour = "Full model")) +
  geom_point(aes(x = 0:9, y = bias_point_size[,4], colour = "Full model")) +
  geom_line(aes(x = 0:9, y = bias_point_size_noint[,4], colour = "Simple model")) +
  geom_point(aes(x = 0:9, y = bias_point_size_noint[,4], colour = "Simple model")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Model type", values = c("Full model"= "green", "Simple model" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Sample size = 1000")+ 
  ylim(-0.02,0.02)+ 
  theme(plot.title = element_text(size=10))
p21 <- ggplot()+
  geom_line(aes(x = 0:9, y = bias_point_size[,6], colour = "Full model")) +
  geom_point(aes(x = 0:9, y = bias_point_size[,6], colour = "Full model")) +
  geom_line(aes(x = 0:9, y = bias_point_size_noint[,6], colour = "Simple model")) +
  geom_point(aes(x = 0:9, y = bias_point_size_noint[,6], colour = "Simple model")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Model type", values = c("Full model"= "green", "Simple model" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Sample size = 5000")+ 
  ylim(-0.02,0.02)+ 
  theme(plot.title = element_text(size=10))

figure7 <- ggarrange(p19+ rremove("ylab"),p20+ rremove("ylab"),p21+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figure7 <- annotate_figure(figure7,
                           left = text_grob("Empirical bias",size = 10, rot =  90))

ggarrange(heights = c(4, 4, 4.9),figure5, figure6, figure7, nrow = 3, ncol = 1)

p7 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,1,1], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,1,1], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,2,1], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,2,1], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Confounding = 0.1")+ ylim(-0.03,0.03) + theme(plot.title = element_text(size=10))
p8 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,1,5], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,1,5], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,2,5], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,2,5], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Confounding = 0.5")+ ylim(-0.03,0.03) + theme(plot.title = element_text(size=10))

p9 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,1,9], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,1,9], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,2,9], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,2,9], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Confounding = 0.9")+ ylim(-0.03,0.03) + theme(plot.title = element_text(size=10))

figure1 <- ggarrange(p7+ rremove("ylab"),p8+ rremove("ylab"),p9+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure1 <- annotate_figure(figure1,top = text_grob("Full model",size = 14),
                           left = text_grob("Empirical bias",size = 10, rot =  90))

p10 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,5,1], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,5,1], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,6,1], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,6,1], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Confounding = 0.1")+ ylim(-0.03,0.03) + theme(plot.title = element_text(size=10))
p11 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,5,5], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,5,5], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,6,5], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,6,5], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Confounding = 0.5") + ylim(-0.03,0.03)+ theme(plot.title = element_text(size=10))

p12 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,5,9], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,5,9], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,6,9], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,6,9], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time',
       title = "Confounding = 0.9")+ ylim(-0.03,0.03) + theme(plot.title = element_text(size=10)) 

figure2 <- ggarrange(p10 + rremove("ylab"),p11 + rremove("ylab"),p12 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T,legend = 'bottom')
figure2 <- annotate_figure(figure2,top = text_grob("Simple model",size = 14),
                           left = text_grob("Empirical bias",size = 10, rot =  90))
ggarrange(heights = c(4, 4.9),figure1, figure2, nrow = 2, ncol = 1)

p1 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,3,1], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,3,1], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,4,1], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,4,1], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Treatment prevalence = 0.1")+ ylim(-0.05,0.12) + theme(plot.title = element_text(size=10))
p2 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,3,5], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,3,5], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,4,5], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,4,5], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Treatment prevalence = 0.5")+ ylim(-0.05,0.12) + theme(plot.title = element_text(size=10))

p3 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,3,9], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,3,9], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,4,9], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,4,9], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Treatment prevalence = 0.9")+ ylim(-0.05,0.12) + theme(plot.title = element_text(size=10))

figure3 <- ggarrange(p1+ rremove("ylab"),p2+ rremove("ylab"),p3+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure3 <- annotate_figure(figure3,top = text_grob("Full model",size = 14),
                           left = text_grob("Empirical bias",size = 10, rot =  90))

p4 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,7,1], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,7,1], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,8,1], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,8,1], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Treatment prevalence = 0.1")+ ylim(-0.05,0.12) + theme(plot.title = element_text(size=10))
p5 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,7,5], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,7,5], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,8,5], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,8,5], colour = "Sandwich")) +
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Empirical bias",
       title = "Treatment prevalence = 0.5") + ylim(-0.05,0.12)+ theme(plot.title = element_text(size=10))

p6 <- ggplot() +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,7,9], colour = "Bootstrap")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,7,9], colour = "Bootstrap")) +
  geom_line(aes(x = 0:9, y = empirical_bias_mid[,8,9], colour = "Sandwich")) +
  geom_point(aes(x = 0:9, y = empirical_bias_mid[,8,9], colour = "Sandwich")) +
  
  scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time',
       title = "Treatment prevalence = 0.9")+ ylim(-0.05,0.12) + theme(plot.title = element_text(size=10)) 

figure4 <- ggarrange(p4 + rremove("ylab"),p5 + rremove("ylab"),p6 + rremove("ylab"), nrow = 1, ncol = 3, common.legend = T,legend = 'bottom')
figure4 <- annotate_figure(figure4,top = text_grob("Simple model",size = 14),
                           left = text_grob("Empirical bias",size = 10, rot =  90))
ggarrange(heights = c(4, 4.9),figure3, figure4, nrow = 2, ncol = 1)

