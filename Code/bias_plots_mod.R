library(modelr)
library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/MPhil-dissertation/Code")
source("simulate_MSM_modified.R")
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

load("true_value_conf_modified.rda")
load("true_value_treat_modified.rda")

coefs <- array(, dim = c(10,1000,9))
treat <- array(, dim = c(10,1000,9))
size <- array(, dim = c(10,1000,6))
coefs_noint <- array(, dim = c(10,1000,9))
treat_noint <- array(, dim = c(10,1000,9))
size_noint <- array(, dim = c(10,1000,6))

for (i in 1:9){
  load(paste0("HPC output/coefs_PP_mod_",i,".rda"))
  load(paste0("HPC output/treat_PP_mod_",i,".rda"))
  load(paste0("HPC output/coefs_PP_noint_mod_", i, ".rda"))
  load(paste0("HPC output/treat_PP_noint_mod_", i, ".rda"))
  coefs[,,i] <- coefs_PP_mod
  treat[,,i] <- treat_PP_mod
  coefs_noint[,,i] <- coefs_PP_noint_mod
  treat_noint[,,i] <- treat_PP_noint_mod
}

sizes <- c(2,5,8,10,25,50)
for (i in 1:6){
  load(paste0("HPC output/size_PP_mod_",sizes[i],".rda"))
  load(paste0("HPC output/size_PP_noint_mod_",sizes[i],".rda"))
  size[,,i] <- size_PP_mod
  size_noint[,,i] <- size_PP_noint_mod
}

bias_point_coefs <- array(,dim = c(10,9))
bias_point_treat <- array(,dim = c(10,9))
bias_point_size <- array(,dim = c(10,6))

bias_point_coefs_noint <- 1*bias_point_coefs
bias_point_treat_noint <- 1*bias_point_coefs
bias_point_size_noint <- 1*bias_point_size

for (i in 1:9){
  bias_point_coefs[,i] <- rowMeans(coefs[,,i], na.rm = TRUE) - (true_value_conf_modified[,1,i] - true_value_conf_modified[,2,i])
  bias_point_treat[,i] <- rowMeans(treat[,,i], na.rm = TRUE) - (true_value_treat_modified[,1,i] - true_value_treat_modified[,2,i])
  bias_point_coefs_noint[,i] <- rowMeans(coefs_noint[,,i], na.rm = TRUE) - (true_value_conf_modified[,1,i] - true_value_conf_modified[,2,i])
  bias_point_treat_noint[,i] <- rowMeans(treat_noint[,,i], na.rm = TRUE) - (true_value_treat_modified[,1,i] - true_value_treat_modified[,2,i])
}

for (i in 1:6){
  bias_point_size[,i] <- rowMeans(size[,,i], na.rm = TRUE) - (true_value_conf_modified[,1,5] - true_value_conf_modified[,2,5])
  bias_point_size_noint[,i] <- rowMeans(size_noint[,,i], na.rm = TRUE) - (true_value_conf_modified[,1,5] - true_value_conf_modified[,2,5])
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

