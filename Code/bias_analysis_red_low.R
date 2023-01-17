library(modelr)
library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/MPhil-dissertation/Code")
source("simulate_MSM_simplified.R")
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
library(matrixStats)

load("HPC output/true_value_conf_red_low.rda")
load("HPC output/true_value_treat_red_low.rda")

coefs <- array(, dim = c(5,1000,9))
treat <- array(, dim = c(5,1000,9))
size <- array(, dim = c(5,1000,7))

treat_vals <- c(-10,-8,-5,-2,0,2,5,8,10)

for (i in 1:9){
  load(paste0("HPC output/est_conf_red_low_",i,".rda"))
  load(paste0("HPC output/est_treat_red_low_",i,".rda"))
  coefs[,,i] <- est_conf
  treat[,,i] <- est_treat
}

# 
# sizes <- c(2,5,8,10,25,50,100)
# for (i in 1:7){
#   load(paste0("HPC output/est_size_red_low_",sizes[i],".rda"))
#   size[,,i] <- est_size
# }

bias_point_coefs <- array(,dim = c(5,9))
bias_point_treat <- array(,dim = c(5,9))
bias_point_size <- array(0,dim = c(5,7))

coefs_sd <- 1.0*bias_point_coefs
treat_sd <- 1.0*bias_point_coefs
size_sd <- 1.0*bias_point_size


for (i in 1:9){
  bias_point_coefs[,i] <- rowMeans(coefs[,,i], na.rm = TRUE) - (true_value_conf_red[,1,i] - true_value_conf_red[,2,i])
  bias_point_treat[,i] <- rowMeans(treat[,,i], na.rm = TRUE) - (true_value_treat_red[,1,i] - true_value_treat_red[,2,i])
  coefs_sd[,i] <- rowSds(coefs[,,i])
  treat_sd[,i] <- rowSds(treat[,,i])
}

# for (i in 1:7){
#   bias_point_size[,i] <- rowMeans(size[,,i], na.rm = TRUE) - (true_value_conf_red[,1,5] - true_value_conf_red[,2,5])
#   size_sd[,i] <- rowSds(size[,,i])
# }

################# BIAS TO SE RATIO PLOTS ##########################
p1 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_coefs[,1]/coefs_sd[,1]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_coefs[,1]/coefs_sd[,1]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Confounding = 0.1")+ 
  ylim(0,2.0)+ 
  theme(plot.title = element_text(size=10))
p2 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_coefs[,5]/coefs_sd[,5]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_coefs[,5]/coefs_sd[,5]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Confounding = 0.5")+
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
p3 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_coefs[,9]/coefs_sd[,9]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_coefs[,9]/coefs_sd[,9]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Confounding = 0.9")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
figure1 <- ggarrange(p1+ rremove("ylab"),p2+ rremove("ylab"),p3+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure1 <- annotate_figure(figure1,
                           left = text_grob("Bias/SE",size = 10, rot =  90),
                           top = text_grob("|Bias/SD| of estimates from simulations based on reduced complexity data", size = 12))
p4 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_treat[,1]/treat_sd[,1]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_treat[,1]/treat_sd[,1]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Treatment prevalence = -1")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
p5 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_treat[,5]/treat_sd[,5]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_treat[,5]/treat_sd[,5]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Treatment prevalence = 0")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
p6 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_treat[,9]/treat_sd[,9]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_treat[,9]/treat_sd[,9]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Treatment prevalence = 1")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
figure2 <- ggarrange(p4+ rremove("ylab"),p5+ rremove("ylab"),p6+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure2 <- annotate_figure(figure2,
                           left = text_grob("Bias/SE",size = 10, rot =  90))
p7 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_size[,1]/size_sd[,1]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_size[,1]/size_sd[,1]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Size = 200")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
p8 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_size[,4]/size_sd[,4]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_size[,4]/size_sd[,4]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Size = 1000")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
p9 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_size[,6]/size_sd[,6]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_size[,6]/size_sd[,6]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Size = 5000")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
figure3 <- ggarrange(p7+ rremove("ylab"),p8+ rremove("ylab"),p9+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figure3 <- annotate_figure(figure3,
                           left = text_grob("Bias/SE",size = 10, rot =  90))
ggarrange(heights = c(4.5, 4, 4.9),figure1, figure2, figure3, nrow = 3, ncol = 1)


################## BIAS PLOTS #########################
p10 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_coefs[,1], colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_coefs[,1], colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Confounding = 0.1")+ 
  ylim(-0.3,0.1)+ 
  theme(plot.title = element_text(size=10))
p11 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_coefs[,5], colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_coefs[,5], colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Confounding = 0.5")+ 
  ylim(-0.3,0.1)+ 
  theme(plot.title = element_text(size=10))
p12 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_coefs[,9], colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_coefs[,9], colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Confounding = 0.9")+ 
  ylim(-0.3,0.1)+ 
  theme(plot.title = element_text(size=10))
figure4 <- ggarrange(p10+ rremove("ylab"),p11+ rremove("ylab"),p12+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure4 <- annotate_figure(figure4,
                           left = text_grob("Empirical bias",size = 10, rot =  90),
                           top = text_grob("Empirical bias of estimates from simulations based on reduced complexity data", size = 12))
p13 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_treat[,1], colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_treat[,1], colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Treatment prevalence = -1")+ 
  ylim(-0.3,0.1)+ 
  theme(plot.title = element_text(size=10))
p14 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_treat[,5], colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_treat[,5], colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Treatment prevalence = 0")+ 
  ylim(-0.3,0.1)+ 
  theme(plot.title = element_text(size=10))
p15 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_treat[,9], colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_treat[,9], colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Treatment prevalence = 1")+ 
  ylim(-0.3,0.1)+ 
  theme(plot.title = element_text(size=10))
figure5 <- ggarrange(p13+ rremove("ylab"),p14+ rremove("ylab"),p15+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure5 <- annotate_figure(figure5,
                           left = text_grob("Empirical bias",size = 10, rot =  90))
p16 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_size[,1], colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_size[,1], colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Size = 200")+ 
  ylim(-0.3,0.1)+ 
  theme(plot.title = element_text(size=10))
p17 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_size[,4], colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_size[,4], colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Size = 1000")+ 
  ylim(-0.3,0.1)+ 
  theme(plot.title = element_text(size=10))
p18 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_size[,6], colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_size[,6], colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Size = 5000")+ 
  ylim(-0.3,0.1)+ 
  theme(plot.title = element_text(size=10))
figure6 <- ggarrange(p16+ rremove("ylab"),p17+ rremove("ylab"),p18+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figure6 <- annotate_figure(figure6,
                           left = text_grob("Empirical bias",size = 10, rot =  90))
ggarrange(heights = c(4.5, 4, 4.9),figure4, figure5, figure6, nrow = 3, ncol = 1)


############## SE PLOTS ####################
p19 <- ggplot()+
  geom_line(aes(x = 0:4, y = coefs_sd[,1])) +
  geom_point(aes(x = 0:4, y = coefs_sd[,1])) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Confounding = 0.1")+ 
  ylim(0,0.5)+ 
  theme(plot.title = element_text(size=10))
p20 <- ggplot()+
  geom_line(aes(x = 0:4, y = coefs_sd[,5])) +
  geom_point(aes(x = 0:4, y = coefs_sd[,5])) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Confounding = 0.5")+ 
  ylim(0,0.5)+ 
  theme(plot.title = element_text(size=10))
p21 <- ggplot()+
  geom_line(aes(x = 0:4, y = coefs_sd[,9])) +
  geom_point(aes(x = 0:4, y = coefs_sd[,9])) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Confounding = 0.9")+ 
  ylim(0,0.5)+ 
  theme(plot.title = element_text(size=10))
figure7 <- ggarrange(p19+ rremove("ylab"),p20+ rremove("ylab"),p21+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure7 <- annotate_figure(figure7,
                           left = text_grob("Standard deviation",size = 10, rot =  90),
                           top = text_grob("Standard deviation of estimates from simulations based on reduced complexity data", size = 12))
p22 <- ggplot()+
  geom_line(aes(x = 0:4, y = treat_sd[,1])) +
  geom_point(aes(x = 0:4, y = treat_sd[,1])) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Treatment prevalence = -1")+ 
  ylim(0,0.5)+ 
  theme(plot.title = element_text(size=10))
p23 <- ggplot()+
  geom_line(aes(x = 0:4, y = treat_sd[,5])) +
  geom_point(aes(x = 0:4, y = treat_sd[,5])) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(0,0.5)+ 
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Treatment prevalence = 0")+ 
  theme(plot.title = element_text(size=10))
p24 <- ggplot()+
  geom_line(aes(x = 0:4, y = treat_sd[,9])) +
  geom_point(aes(x = 0:4, y = treat_sd[,9])) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(0,0.5)+ 
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Treatment prevalence = 1")+ 
  theme(plot.title = element_text(size=10))

figure8 <- ggarrange(p22+ rremove("ylab"),p23+ rremove("ylab"),p24+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure8 <- annotate_figure(figure8,
                           left = text_grob("Standard deviation",size = 10, rot =  90))
p25 <- ggplot()+
  geom_line(aes(x = 0:4, y = size_sd[,1])) +
  geom_point(aes(x = 0:4, y = size_sd[,1])) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(0,0.5)+ 
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Size = 200")+ 
  theme(plot.title = element_text(size=10))
p26 <- ggplot()+
  geom_line(aes(x = 0:4, y = size_sd[,4])) +
  geom_point(aes(x = 0:4, y = size_sd[,4])) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(0,0.5)+ 
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Size = 1000")+ 
  theme(plot.title = element_text(size=10))
p27 <- ggplot()+
  geom_line(aes(x = 0:4, y = size_sd[,6])) +
  geom_point(aes(x = 0:4, y = size_sd[,6])) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(0,0.5)+ 
  labs(x = 'Follow up time', 
       y = "Bias",
       title = "Size = 5000")+ 
  theme(plot.title = element_text(size=10))

figure9 <- ggarrange(p25+ rremove("ylab"),p26+ rremove("ylab"),p27+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figure9 <- annotate_figure(figure9,
                           left = text_grob("Standard deviation",size = 10, rot =  90))
ggarrange(heights = c(4.5, 4, 4),figure7, figure8, figure9, nrow = 3, ncol = 1)

##################### MSE ##############################
p28 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_coefs[,1]^2 + coefs_sd[,1]^2, colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_coefs[,1]^2 + coefs_sd[,1]^2, colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Confounding = 0.1")+ 
  #ylim(0,2.0)+ 
  theme(plot.title = element_text(size=10))
p29 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_coefs[,5]^2 + coefs_sd[,5]^2, colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = bias_point_coefs[,5]^2 + coefs_sd[,5]^2, colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Confounding = 0.5")+
  #ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
p30 <- ggplot()+
  geom_line(aes(x = 0:4, y = bias_point_coefs[,9]^2 + coefs_sd[,9], colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_coefs[,9]^2 + coefs_sd[,9]), colour = "KM curve")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Confounding = 0.9")+ 
  #ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
figure10 <- ggarrange(p28+ rremove("ylab"),p29+ rremove("ylab"),p30+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure10 <- annotate_figure(figure1,
                            left = text_grob("Bias/SE",size = 10, rot =  90),
                            top = text_grob("MSE of estimates from simulations based on reduced complexity data", size = 12))
p31 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_treat[,1]/treat_sd[,1]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_treat[,1]/treat_sd[,1]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Treatment prevalence = -1")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
p32 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_treat[,5]/treat_sd[,5]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_treat[,5]/treat_sd[,5]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Treatment prevalence = 0")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
p33 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_treat[,9]/treat_sd[,9]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_treat[,9]/treat_sd[,9]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Treatment prevalence = 1")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
figure11 <- ggarrange(p31+ rremove("ylab"),p32+ rremove("ylab"),p33+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'none')
figure11 <- annotate_figure(figure11,
                            left = text_grob("Bias/SE",size = 10, rot =  90))
p34 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_size[,1]/size_sd[,1]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_size[,1]/size_sd[,1]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Size = 200")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
p35 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_size[,4]/size_sd[,4]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_size[,4]/size_sd[,4]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Size = 1000")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
p36 <- ggplot()+
  geom_line(aes(x = 0:4, y = abs(bias_point_size[,6]/size_sd[,6]), colour = "KM curve")) +
  geom_point(aes(x = 0:4, y = abs(bias_point_size[,6]/size_sd[,6]), colour = "KM curve")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Bias calculation method", values = c("KM curve"= "green", "Limiting value" = "purple")) +
  labs(x = 'Follow up time', 
       y = "Bias to SE ratio",
       title = "Size = 5000")+ 
  ylim(0,2.0)+
  theme(plot.title = element_text(size=10))
figure12 <- ggarrange(p34+ rremove("ylab"),p35+ rremove("ylab"),p36+ rremove("ylab"), nrow = 1, ncol = 3, common.legend = T, legend = 'bottom')
figure12 <- annotate_figure(figure12,
                            left = text_grob("Bias/SE",size = 10, rot =  90))
ggarrange(heights = c(4.5, 4, 4.9),figure1, figure2, figure3, nrow = 3, ncol = 1)
