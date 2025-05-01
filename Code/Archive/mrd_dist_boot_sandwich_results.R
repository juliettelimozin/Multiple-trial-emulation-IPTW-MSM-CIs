library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/Multiple-trial-emulation-IPTW-MSM-CIs/Code")
library(ggplot2)
library(ggpubr)
load("HPC output/true_value_red_newsimus.rda")
library(modelr)
library(tidyverse)
library(tidyr)
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
library(latex2exp)
library(grDevices)
library(xtable)

bootstrap_iter <- 500
sampling_size <- 200


hr_bootstrap <- array(,dim = c(bootstrap_iter,27))
hr_sandwich <- array(,dim = c(sampling_size,27))
mrd_bootstrap <- array(,dim = c(5,bootstrap_iter,27))
mrd_sandwich <- array(,dim = c(5,sampling_size,27))
mrd_LEF_outcome <- array(,dim = c(5,bootstrap_iter,27))
mrd_LEF_both <- array(,dim = c(5,bootstrap_iter,27))
mrd_jackknife_mvn <- array(,dim = c(5,sampling_size,9))
mrd_jackknife_wald <- array(,dim = c(5,200,9))
SE_true <- array(,dim = c(5,27))
SE_bootstrap <- array(,dim = c(5,27))
SE_sandwich <- array(,dim = c(5,27))
SE_LEF_outcome <- array(,dim = c(5,27))
SE_LEF_both <- array(,dim = c(5,27))
SE_jackknife_mvn <- array(,dim = c(5,9))
SE_jackknife_wald <- array(, dim = c(5,9))

est <- array(,dim = c(5,1000,27))


size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)

scenarios <- as.data.frame(tidyr::crossing(size,conf, treat))
true_value_red <- -true_value_red

for (i in 1:27){
  load(paste0("HPC output/hr_estimates_boot_low_",i, ".rda"))
  load(paste0("HPC output/hr_estimates_sandwich_low_",i, ".rda"))
  load(paste0("HPC output/mrd_estimates_boot_low_",i, ".rda"))
  load(paste0("HPC output/mrd_estimates_sandwich_low_",i, ".rda"))
  load(paste0("HPC output/mrd_estimates_LEF_outcome_low_",i, ".rda"))
  load(paste0("HPC output/mrd_estimates_LEF_both_low_",i, ".rda"))
  if(i %in% 1:9){
    load(paste0("HPC output/mrd_estimates_jackknife_wald_low_",i, ".rda"))
    load(paste0("HPC output/mrd_estimates_jackknife_mvn_low_",i, ".rda"))
    load(paste0("HPC output/jackknife_mrd_se_low_",i, ".rda"))
    mrd_jackknife_mvn[,,i] <- mrd_estimates_jackknife_mvn
    mrd_jackknife_wald[,,i] <- mrd_estimates_jackknife_wald
    SE_jackknife_wald[,i] <- jackknife_mrd_se
    SE_jackknife_mvn[,i] <- rowSds(mrd_estimates_jackknife_mvn, na.rm = TRUE)
  }
  est[,,i] <- -estimates
  hr_bootstrap[,i] <- hr_estimates_boot
  hr_sandwich[,i] <- hr_estimates_sandwich
  mrd_bootstrap[,,i] <- mrd_estimates_boot
  mrd_sandwich[,,i] <- mrd_estimates_sandwich
  mrd_LEF_outcome[,,i] <- mrd_estimates_LEF_outcome
  mrd_LEF_both[,,i] <- mrd_estimates_LEF_both
  
  SE_true[,i] <- rowSds(-estimates, na.rm = TRUE)
  SE_bootstrap[,i] <- rowSds(mrd_estimates_boot, na.rm = TRUE)
  if (any(!is.na(mrd_estimates_sandwich)) == T){
    SE_sandwich[,i] <- rowSds(mrd_estimates_sandwich, na.rm = TRUE)
  }
  SE_LEF_outcome[,i] <- rowSds(mrd_estimates_LEF_outcome, na.rm = TRUE)
  SE_LEF_both[,i] <- rowSds(mrd_estimates_LEF_both, na.rm = TRUE)
}


SE_ratio_low <- lapply(1:27, function(i){
  if(i %in% 1:9){
    ggplot() +
      geom_line(aes(x = 0:4, y = SE_bootstrap[,i]/SE_true[,i], colour = 'Bootstrap')) +
      geom_point(aes(x = 0:4, y = SE_bootstrap[,i]/SE_true[,i], colour = 'Bootstrap')) +
      geom_line(aes(x = 0:4, y = SE_sandwich[,i]/SE_true[,i], colour = 'Sandwich')) +
      geom_point(aes(x = 0:4, y = SE_sandwich[,i]/SE_true[,i], colour = 'Sandwich')) +
      geom_line(aes(x = 0:4, y = SE_LEF_outcome[,i]/SE_true[,i], colour = 'LEF outcome')) +
      geom_point(aes(x = 0:4, y = SE_LEF_outcome[,i]/SE_true[,i], colour = 'LEF outcome')) +
      geom_line(aes(x = 0:4, y = SE_LEF_both[,i]/SE_true[,i], colour = 'LEF both')) +
      geom_point(aes(x = 0:4, y = SE_LEF_both[,i]/SE_true[,i], colour = 'LEF both')) +
      ylim(0,22) +
       geom_line(aes(x = 0:4, y = SE_jackknife_wald[,i]/SE_true[,i], colour = 'Jackknife Wald')) +
       geom_point(aes(x = 0:4, y = SE_jackknife_wald[,i]/SE_true[,i], colour = 'Jackknife Wald')) +
      # geom_line(aes(x = 0:4, y = SE_jackknife_mvn[,i]/SE_true[,i], colour = 'Jackknife MVN')) +
      # geom_point(aes(x = 0:4, y = SE_jackknife_mvn[,i]/SE_true[,i], colour = 'Jackknife MVN')) +
      xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
      ylab(bquote(SD(hat(MRD_CI))/SD(hat(MRD))))+
      theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
      scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                        "LEF outcome" = "green", "LEF both" = "purple",
                                                      'Jackknife Wald' = 'pink', 'Jackknife MVN' = 'orange'))
  } else {
    ggplot() +
      geom_line(aes(x = 0:4, y = SE_bootstrap[,i]/SE_true[,i], colour = 'Bootstrap')) +
      geom_point(aes(x = 0:4, y = SE_bootstrap[,i]/SE_true[,i], colour = 'Bootstrap')) +
      geom_line(aes(x = 0:4, y = SE_sandwich[,i]/SE_true[,i], colour = 'Sandwich')) +
      geom_point(aes(x = 0:4, y = SE_sandwich[,i]/SE_true[,i], colour = 'Sandwich')) +
      geom_line(aes(x = 0:4, y = SE_LEF_outcome[,i]/SE_true[,i], colour = 'LEF outcome')) +
      geom_point(aes(x = 0:4, y = SE_LEF_outcome[,i]/SE_true[,i], colour = 'LEF outcome')) +
      geom_line(aes(x = 0:4, y = SE_LEF_both[,i]/SE_true[,i], colour = 'LEF both')) +
      geom_point(aes(x = 0:4, y = SE_LEF_both[,i]/SE_true[,i], colour = 'LEF both')) +
      xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
      ylab(bquote(SD(hat(MRD_CI))/SD(hat(MRD))))+
      ylim(0,3.2) +
      theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
      scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                      "LEF outcome" = "green", "LEF both" = "purple",
                                                      'Jackknife Wald' = 'pink', 'Jackknife MVN' = 'orange'))
  }
})
annotate_figure(ggarrange(plotlist = SE_ratio_low, nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

mrd_dist1_low <- lapply(c(21,24,27), function(i){
  ggplot() +
    geom_histogram(aes(x = est[1,,i], fill="True"),binwidth = 0.001)+
    geom_histogram(aes(x = mrd_bootstrap[1,,i], fill="Bootstrap"), binwidth = .001 )+
    geom_histogram(aes(x = mrd_sandwich[1,,i], fill="Sandwich"), binwidth = .001)+
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Distribution", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                      "True" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist1_low, nrow = 1, ncol = 3, common.legend = T, legend = 'bottom'))

mrd_dist2_low <- lapply(c(21,24,27), function(i){
  ggplot() +
    geom_histogram(aes(x = est[2,,i], fill="True"),binwidth = 0.001)+
    geom_histogram(aes(x = mrd_bootstrap[2,,i], fill="Bootstrap"), binwidth = .001 )+
    geom_histogram(aes(x = mrd_sandwich[2,,i], fill="Sandwich"), binwidth = .001)+
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Distribution", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                      "True" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist2_low, nrow = 1, ncol = 3, common.legend = T, legend = 'bottom'))

mrd_dist3_low <- lapply(c(21,24,27), function(i){
  ggplot() +
    geom_histogram(aes(x = est[3,,i], fill="True"),binwidth = 0.001)+
    geom_histogram(aes(x = mrd_bootstrap[3,,i], fill="Bootstrap"), binwidth = .001 )+
    geom_histogram(aes(x = mrd_sandwich[3,,i], fill="Sandwich"), binwidth = .001)+
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Distribution", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                      "True" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist3_low, nrow = 1, ncol = 3, common.legend = T, legend = 'bottom'))

mrd_dist4_low <- lapply(c(21,24,27), function(i){
  ggplot() +
    geom_histogram(aes(x = est[4,,i], fill="True"),binwidth = 0.001)+
    geom_histogram(aes(x = mrd_bootstrap[4,,i], fill="Bootstrap"), binwidth = .001 )+
    geom_histogram(aes(x = mrd_sandwich[4,,i], fill="Sandwich"), binwidth = .001)+
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Distribution", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                      "True" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist4_low, nrow = 1, ncol = 3, common.legend = T, legend = 'bottom'))

mrd_dist5_low <- lapply(c(21,24,27), function(i){
  ggplot() +
    geom_histogram(aes(x = est[5,,i], fill="True"),binwidth = 0.001)+
    geom_histogram(aes(x = mrd_bootstrap[5,,i], fill="Bootstrap"), binwidth = .001 )+
    geom_histogram(aes(x = mrd_sandwich[5,,i], fill="Sandwich"), binwidth = .001)+
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Distribution", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                      "True" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist5_low, nrow = 1, ncol = 3, common.legend = T, legend = 'bottom'))

annotate_figure(ggarrange(plotlist = c(mrd_dist1_low,mrd_dist2_low,mrd_dist3_low,mrd_dist4_low,mrd_dist5_low), nrow = 5, ncol = 3, common.legend = T, legend = 'bottom'))

