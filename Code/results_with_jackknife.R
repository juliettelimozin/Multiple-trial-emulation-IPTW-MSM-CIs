library(dplyr)
library(tidyverse)
library(tidyr)
setwd("~/rds/hpc-work/Multiple-trial-emulation-IPTW-MSM-CIs/Code")
library(ggplot2)
library(ggpubr)
load("true_value_red_newsimus.rda")
true_value_red <- -true_value_red
load("true_value_surv0.rda")
load("true_value_surv1.rda")
#load("true_value_boot_200it_200k_fixed_700it.rda")
#true_value_red <- true_value_boot_200it_200k_fixed
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

bootstrap <- array(,dim = c(5,2,iters,27,3))
LEF_outcome <- array(,dim = c(5,2,iters,27,3))
LEF_both <- array(,dim = c(5,2,iters,27,3))
sandwich <- array(,dim = c(5,2,iters,27,3))
jackknife_mvn <- array(,dim = c(5,2,iters,27,3))
jackknife_wald <- array(,dim = c(5,2,iters,27,3))
time <- array(,dim = c(6,iters,27,3))
est <- array(,dim = c(5,iters,27,3))
treat_pos <- c(-1,-0.8,-0.5,-0.2,0,0.2,0.5,0.8,1)
outcomes <- c("low", 'med', 'high')
jackknife_se <- array(,dim = c(5,iters,27,3))

size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)

scenarios <- as.data.frame(tidyr::crossing(size,conf, treat))
bias_point <- array(,dim = c(5,27,3))
bias_surv0 <- array(,dim = c(5,27,3))
bias_surv1 <- array(,dim = c(5,27,3))
sd_point <- array(,dim = c(5,27,3))
mean_time <- data.frame(matrix(,nrow = 0, ncol = 10))

na_failure_rate <- data.frame(matrix(,nrow = 0, ncol = 7))
se_ratio <- data.frame(matrix(,nrow = 0, ncol = 5))

load("~/rds/hpc-work/Project1/NewSimusJ/J_sandwich_SE.rda")
load("~/rds/hpc-work/Project1/NewSimusJ/J_bootstrap_SE.rda")
load("~/rds/hpc-work/Project1/NewSimusJ/J_LEF_outcome_SE.rda")
load("~/rds/hpc-work/Project1/NewSimusJ/J_LEF_both_SE.rda")
load("~/rds/hpc-work/Project1/NewSimusJ/J_jackknife_mvn_SE.rda")

for (i in 10:27){
  for (j in 1:3){
    load(paste0("~/rds/hpc-work/Project1/NewSimusJ/J_CI_bootstrap_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("~/rds/hpc-work/Project1/NewSimusJ/J_CI_jackknife_mvn_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("~/rds/hpc-work/Project1/NewSimusJ/J_CI_jackknife_wald_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("~/rds/hpc-work/Project1/NewSimusJ/J_CI_LEF_outcome_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("~/rds/hpc-work/Project1/NewSimusJ/J_CI_LEF_both_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("~/rds/hpc-work/Project1/NewSimusJ/J_CI_sandwich_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("~/rds/hpc-work/Project1/NewSimusJ/J_computation_time_",outcomes[j],'_', i, ".rda"))
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ/J_estimates_red_',outcomes[j],'_',i, '.rda'))
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ/J_survival_treatment_estimates_',outcomes[j],'_',i, '.rda'))
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ/J_survival_control_estimates_',outcomes[j],'_',i, '.rda'))
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ/J_jackknife_SEs_',outcomes[j],'_',i, '.rda'))
    
    bootstrap[,,,i,j] <- CI_bootstrap_PP_red
    LEF_outcome[,,,i,j] <- CI_LEF_outcome_PP_red
    LEF_both[,,,i,j] <- CI_LEF_both_PP_red
    sandwich[,,,i,j] <- CI_sandwich_PP_red
    jackknife_mvn[,,,i,j] <- CI_jackknife_mvn_PP_red
    jackknife_wald[,,,i,j] <- CI_jackknife_wald_PP_red
    jackknife_se[,,i,j] <- jackknife_SEs
    
    scenario <- i%%9
    est[,,i,j] <- estimates
    if (scenario ==0){scenario <- 9}
    bias_point[,i,j] <- rowMeans(estimates, na.rm = TRUE) - true_value_red[,scenario,j]
    bias_surv0[,i,j] <- rowMeans(survival_control_estimates, na.rm = TRUE) - surv0[,scenario,j]
    bias_surv1[,i,j] <- rowMeans(survival_treatment_estimates, na.rm = TRUE) - surv1[,scenario,j]
    
    sd_point[,i,j] <- rowSds(estimates, na.rm = TRUE)
    
    for (k in 1:5){
      se_ratio <- rbind(se_ratio, cbind(k-1,'Bootstrap',bootstrap_SE[k,,i,j]/sd_point[k,i,j],i,j))
      se_ratio <- rbind(se_ratio, cbind(k-1,'LEF outcome',LEF_outcome_SE[k,,i,j]/sd_point[k,i,j],i,j))
      
      se_ratio <- rbind(se_ratio, cbind(k-1,'LEF both',LEF_both_SE[k,,i,j]/sd_point[k,i,j],i,j))
      
      se_ratio <- rbind(se_ratio, cbind(k-1,'Sandwich',sandwich_SE[k,,i,j]/sd_point[k,i,j],i,j))
      
      se_ratio <- rbind(se_ratio, cbind(k-1,'Jackknife MVN',jackknife_mvn_SE[k,,i,j]/sd_point[k,i,j],
                                        i,j))
      
      se_ratio <- rbind(se_ratio, cbind(k-1,'Jackknife Wald',jackknife_se[k,,i,j]/sd_point[k,i,j],i,j))
    }
    
    mean_time <- rbind(mean_time, c(outcomes[j], scenarios[i,1], scenarios[i,2],scenarios[i,3], rowMeans(computation_time, na.rm = TRUE)))
    
    
  }
}

colnames(mean_time) <- c('Outcome_prevalence', 'Sample_size', 'Confounding', 'Treatment_prevalence', 'Bootstrap', 'LEF_outcome',
                         'LEF_both', 'Jackknife_Wald','Jackknife_MVN' ,'Sandwich')

colnames(se_ratio) <- c('Visit', 'CI_type','SE_ratio','i', 'j')
se_ratio$i <- as.numeric(se_ratio$i)
se_ratio$SE_ratio <- as.numeric(se_ratio$SE_ratio)
se_ratio$Visit <- as.factor(se_ratio$Visit)
se_ratio$j <- as.numeric(se_ratio$j)

mean_time <- mean_time %>% 
  dplyr::mutate(Sample_size = as.numeric(Sample_size),
                Confounding = as.numeric(Confounding),
                Treatment_prevalence = as.numeric(Treatment_prevalence),
                Bootstrap = as.numeric(Bootstrap),
                LEF_outcome = as.numeric(LEF_outcome),
                LEF_both = as.numeric(LEF_both),
                Jackknife_Wald = as.numeric(Jackknife_Wald),
                Jackknife_MVN = as.numeric(Jackknife_MVN),
                Sandwich = as.numeric(Sandwich)) %>% 
  dplyr::group_by(Outcome_prevalence,Sample_size) %>% 
  dplyr::summarise(Bootstrap = mean(Bootstrap),
                   LEF_outcome = mean(LEF_outcome),
                   LEF_both = mean(LEF_both),
                   Jackknife_Wald =mean(Jackknife_Wald),
                   Jackknife_MVN = mean(Jackknife_MVN),
                   Sandwich = mean(Sandwich))
print(xtable(mean_time), 
             type = 'latex',include.rownames=FALSE)

################ EMPIRICAL MRD SE ####################
mrd_se_quantiles_low <-  lapply(1:9, function(i){
  ggplot(se_ratio[se_ratio$i == i & se_ratio$j == 1 & se_ratio$CI_type != 'Jackknife MVN' & se_ratio$CI_type != 'Jackknife Wald',]) +
    stat_summary(
      mapping = aes(x = Visit, y = SE_ratio, colour = CI_type),
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = mean,
      size=0.3, 
      position = position_dodge(width = 0.5))+
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    #ylim(0,2.25) + 
    geom_hline(yintercept = 1, linetype = "dashed")+
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:9){
  if(i %in% 1:3){
      mrd_se_quantiles_low[[i]] <- mrd_se_quantiles_low[[i]] +
        labs(title = bquote(alpha[a] == .(scenarios[i,3])))}
  if(i %in% c(1,4,7)){
    mrd_se_quantiles_low[[i]] <- mrd_se_quantiles_low[[i]] +
      ylab(bquote(atop(alpha[c] == .(scenarios[i,2]), 'Ratio of SEs')))
  } else{mrd_se_quantiles_low[[i]] <- mrd_se_quantiles_low[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())}
  if(i %in% 7:9){
    mrd_se_quantiles_low[[i]] <- mrd_se_quantiles_low[[i]] +
      xlab('Visit')
  } else {mrd_se_quantiles_low[[i]] <- mrd_se_quantiles_low[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())}
  
}

annotate_figure(ggarrange(plotlist = mrd_se_quantiles_low, nrow = 3, ncol = 3, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.2,1,1),
                          heights = c(1.1, 0.95, 1)))

mrd_se_quantiles_low <-  lapply(1:9, function(i){
  ggplot(se_ratio[se_ratio$i == i & se_ratio$j == 2 & se_ratio$CI_type == 'Jackknife MVN',]) +
    stat_summary(
      mapping = aes(x = Visit, y = SE_ratio, colour = CI_type),
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = mean,
      size=0.3, 
      position = position_dodge(width = 0.5))+
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    ylim(0,2.25) + 
    geom_hline(yintercept = 1, linetype = "dashed")+
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:9){
  if(i %in% 1:3){
    mrd_se_quantiles_low[[i]] <- mrd_se_quantiles_low[[i]] +
      labs(title = bquote(alpha[a] == .(scenarios[i,3])))}
  if(i %in% c(1,4,7)){
    mrd_se_quantiles_low[[i]] <- mrd_se_quantiles_low[[i]] +
      ylab(bquote(atop(alpha[c] == .(scenarios[i,2]), 'Ratio of SEs')))
  } else{mrd_se_quantiles_low[[i]] <- mrd_se_quantiles_low[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())}
  if(i %in% 7:9){
    mrd_se_quantiles_low[[i]] <- mrd_se_quantiles_low[[i]] +
      xlab('Visit')
  } else {mrd_se_quantiles_low[[i]] <- mrd_se_quantiles_low[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())}
  
}

annotate_figure(ggarrange(plotlist = mrd_se_quantiles_low, nrow = 3, ncol = 3, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.2,1,1),
                          heights = c(1.1, 0.95, 1)))

mrd_se_quantiles_med <-  lapply(10:27, function(i){
  ggplot(se_ratio[se_ratio$i == i & se_ratio$j == 1 & se_ratio$CI_type != 'Jackknife MVN' & se_ratio$CI_type != 'Jackknife Wald',]) +
    stat_summary(
      mapping = aes(x = Visit, y = SE_ratio, colour = CI_type),
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = mean,
      size=0.3, 
      position = position_dodge(width = 0.5))+
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    #ylim(0,2.25) + 
    geom_hline(yintercept = 1, linetype = "dashed")+
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:18){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      mrd_se_quantiles_med[[i]] <- mrd_se_quantiles_med[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          mrd_se_quantiles_med[[i]] <- mrd_se_quantiles_med[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10)){
    mrd_se_quantiles_med[[i]] <- mrd_se_quantiles_med[[i]] +
      ylab(bquote(atop(n == .(scenarios[i+9,1]), 'Ratio of SEs')))
  } else{mrd_se_quantiles_med[[i]] <- mrd_se_quantiles_med[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 10:18){
    mrd_se_quantiles_med[[i]] <- mrd_se_quantiles_med[[i]] +
      xlab('Visit')
  } else {mrd_se_quantiles_med[[i]] <- mrd_se_quantiles_med[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = mrd_se_quantiles_med, nrow = 2, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 1)))

mrd_se_quantiles_med <-  lapply(1:9, function(i){
  ggplot(se_ratio[se_ratio$i == i & se_ratio$j == 2 & se_ratio$CI_type == 'Jackknife MVN',]) +
    stat_summary(
      mapping = aes(x = Visit, y = SE_ratio, colour = CI_type),
      fun.min = function(z) { quantile(z,0.25, na.rm = TRUE) },
      fun.max = function(z) { quantile(z,0.75, na.rm = TRUE) },
      fun = function(z) { mean(z, na.rm = TRUE) },
      size=0.3, 
      position = position_dodge(width = 0.5),
      na.rm = TRUE)+
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    ylim(0,2.25) + 
    geom_hline(yintercept = 1, linetype = "dashed")+
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:9){
  if(i %in% 1:3){
    mrd_se_quantiles_med[[i]] <- mrd_se_quantiles_med[[i]] +
      labs(title = bquote(alpha[a] == .(scenarios[i,3])))}
  if(i %in% c(1,4,7)){
    mrd_se_quantiles_med[[i]] <- mrd_se_quantiles_med[[i]] +
      ylab(bquote(atop(alpha[c] == .(scenarios[i,2]), 'Ratio of SEs')))
  } else{mrd_se_quantiles_med[[i]] <- mrd_se_quantiles_med[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())}
}

annotate_figure(ggarrange(plotlist = mrd_se_quantiles_med, nrow = 3, ncol = 3, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.2,1,1),
                          heights = c(1.1, 0.95, 1)))

mrd_se_quantiles_high <-  lapply(10:27, function(i){
  ggplot(se_ratio[se_ratio$i == i & se_ratio$j == 3 & se_ratio$CI_type != 'Jackknife MVN',]) +
    stat_summary(
      mapping = aes(x = Visit, y = SE_ratio, colour = CI_type),
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = function(z) { mean(z, na.rm = TRUE) },
      size=0.3, 
      position = position_dodge(width = 0.5))+
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Ratio of SEs") +  
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    ylim(0,2.25) + 
    geom_hline(yintercept = 1, linetype = "dashed")+
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14))
}
)
annotate_figure(ggarrange(plotlist = mrd_se_quantiles_high, nrow = 2, ncol = 9, common.legend = T,
                          legend = 'bottom'))
mrd_se_quantiles_high <-  lapply(1:9, function(i){
  ggplot(se_ratio[se_ratio$i == i & se_ratio$j == 3 & se_ratio$CI_type == 'Jackknife MVN',]) +
    stat_summary(
      mapping = aes(x = Visit, y = SE_ratio, colour = CI_type),
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = function(z) { mean(z, na.rm = TRUE) },
      size=0.3, 
      position = position_dodge(width = 0.5),
      na.rm = TRUE)+
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Ratio of SEs") +  
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    ylim(0,15) + 
    geom_hline(yintercept = 1, linetype = "dashed")+
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14))
}
)
annotate_figure(ggarrange(plotlist = mrd_se_quantiles_high, nrow = 3, ncol = 3, common.legend = T,
                          legend = 'bottom'))

##################### MRD DIST ##############################
mrd_dist1_low <- lapply(1:27, function(i){
  scenario <- i%%9
  if (scenario ==0){scenario <- 9}
  ggplot() +
    geom_histogram(aes(x = est[1,,i,3], fill="High"), binwidth = .001)+
    geom_histogram(aes(x = est[1,,i,2], fill="Medium"), binwidth = .001 )+
    geom_histogram(aes(x = est[1,,i,1], fill="Low"),binwidth = 0.001)+
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical MRD dist.") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                      "High" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist1_low, nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

mrd_dist2_low <- lapply(1:27, function(i){
  scenario <- i%%9
  if (scenario ==0){scenario <- 9}
  ggplot() +
    geom_histogram(aes(x = est[2,,i,3], fill="High"), binwidth = .001)+
    
    geom_histogram(aes(x = est[2,,i,2], fill="Medium"), binwidth = .001 )+
    geom_histogram(aes(x = est[2,,i,1], fill="Low"),binwidth = 0.001)+
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical MRD dist.") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                      "High" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist2_low, nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

mrd_dist3_low <- lapply(1:27, function(i){
  scenario <- i%%9
  if (scenario ==0){scenario <- 9}
  ggplot() +
    geom_histogram(aes(x = est[3,,i,3], fill="High"), binwidth = .001)+
    geom_histogram(aes(x = est[3,,i,2], fill="Medium"), binwidth = .001 )+
    geom_histogram(aes(x = est[3,,i,1], fill="Low"),binwidth = 0.001)+
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical MRD dist.") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                      "High" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist3_low, nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

mrd_dist4_low <- lapply(1:27, function(i){
  scenario <- i%%9
  if (scenario ==0){scenario <- 9}
  ggplot() +
    geom_histogram(aes(x = est[4,,i,3], fill="High"), binwidth = .001)+
    geom_histogram(aes(x = est[4,,i,2], fill="Medium"), binwidth = .001 )+
    geom_histogram(aes(x = est[4,,i,1], fill="Low"),binwidth = 0.001)+
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical MRD dist.") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                      "High" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist4_low, nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

mrd_dist5_low <- lapply(1:27, function(i){
  scenario <- i%%9
  if (scenario ==0){scenario <- 9}
  ggplot() +
    geom_histogram(aes(x = est[5,,i,3], fill="High"), binwidth = .001)+
    geom_histogram(aes(x = est[5,,i,2], fill="Medium"), binwidth = .001 )+
    geom_histogram(aes(x = est[5,,i,1], fill="Low"),binwidth = 0.001)+
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical MRD dist.") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                      "High" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist5_low, nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

################BIAS, SD, MSE PLOTS ###################
bias_plots <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_point[,i,1], colour = 'Low')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,1], colour = 'Low')) +
    geom_line(aes(x = 0:4, y = bias_point[,i,2], colour = 'Medium')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,2], colour = 'Medium')) +
    geom_line(aes(x = 0:4, y = bias_point[,i,3], colour = 'High')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,3], colour = 'High')) +
    scale_color_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                       "High" = "green")) +
    ylim(-0.1,0.1) +
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      bias_plots[[i]] <- bias_plots[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          bias_plots[[i]] <- bias_plots[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    bias_plots[[i]] <- bias_plots[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'Empirical Bias')))
  } else{bias_plots[[i]] <- bias_plots[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    bias_plots[[i]] <- bias_plots[[i]] +
      xlab('Visit')
  } else {bias_plots[[i]] <- bias_plots[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = bias_plots, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))

bias_plots_surv0 <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_surv0[,i,1], colour = 'Low')) +
    geom_point(aes(x = 0:4, y = bias_surv0[,i,1], colour = 'Low')) +
    geom_line(aes(x = 0:4, y = bias_surv0[,i,2], colour = 'Medium')) +
    geom_point(aes(x = 0:4, y = bias_surv0[,i,2], colour = 'Medium')) +
    geom_line(aes(x = 0:4, y = bias_surv0[,i,3], colour = 'High')) +
    geom_point(aes(x = 0:4, y = bias_surv0[,i,3], colour = 'High')) +
    #xlab(bquote(paste('N = ',.(scenarios[i,1]),", ", alpha[c] == .(scenarios[i,2]),', ',alpha[a] == .(scenarios[i,3])))) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    
    ylab("Empirical bias") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_color_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                       "High" = "green")) +
    ylim(-0.1,0.1)
})
annotate_figure(ggarrange(plotlist = bias_plots_surv0[1:27], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

bias_plots_surv1 <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_surv1[,i,1], colour = 'Low')) +
    geom_point(aes(x = 0:4, y = bias_surv1[,i,1], colour = 'Low')) +
    geom_line(aes(x = 0:4, y = bias_surv1[,i,2], colour = 'Medium')) +
    geom_point(aes(x = 0:4, y = bias_surv1[,i,2], colour = 'Medium')) +
    geom_line(aes(x = 0:4, y = bias_surv1[,i,3], colour = 'High')) +
    geom_point(aes(x = 0:4, y = bias_surv1[,i,3], colour = 'High')) +
    #xlab(bquote(paste('N = ',.(scenarios[i,1]),", ", alpha[c] == .(scenarios[i,2]),', ',alpha[a] == .(scenarios[i,3])))) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    
    ylab("Empirical bias") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_color_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                       "High" = "green")) +
    ylim(-0.1,0.1)
})
annotate_figure(ggarrange(plotlist = bias_plots_surv1[1:27], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

sd_plots <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = sd_point[,i,1], colour = 'Low')) +
    geom_point(aes(x = 0:4, y = sd_point[,i,1], colour = 'Low')) +
    geom_line(aes(x = 0:4, y = sd_point[,i,2], colour = 'Medium')) +
    geom_point(aes(x = 0:4, y = sd_point[,i,2], colour = 'Medium')) +
    geom_line(aes(x = 0:4, y = sd_point[,i,3], colour = 'High')) +
    geom_point(aes(x = 0:4, y = sd_point[,i,3], colour = 'High')) +
    scale_color_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                       "High" = "green")) +
    ylim(0,0.28) +
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      sd_plots[[i]] <- sd_plots[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          sd_plots[[i]] <- sd_plots[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    sd_plots[[i]] <- sd_plots[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'Empirical SD')))
  } else{sd_plots[[i]] <- sd_plots[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    sd_plots[[i]] <- sd_plots[[i]] +
      xlab('Visit')
  } else {sd_plots[[i]] <- sd_plots[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = sd_plots, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))
bias.se_plots <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_point[,i,1]/sd_point[,i,1], colour = 'Low')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,1]/sd_point[,i,1], colour = 'Low')) +
    geom_line(aes(x = 0:4, y = bias_point[,i,2]/sd_point[,i,2], colour = 'Medium')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,2]/sd_point[,i,2], colour = 'Medium')) +
    geom_line(aes(x = 0:4, y = bias_point[,i,3]/sd_point[,i,3], colour = 'High')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,3]/sd_point[,i,3], colour = 'High')) +
    #xlab(bquote(paste('N = ',.(scenarios[i,1]),", ", alpha[c] == .(scenarios[i,2]),', ',alpha[a] == .(scenarios[i,3])))) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    
    ylab("Bias/SD") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_color_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                       "High" = "green")) +
    ylim(-0.25,0.35)
})
annotate_figure(ggarrange(plotlist = bias.se_plots[1:27], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

mse_plots <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = sqrt(bias_point[,i,1]^2 + sd_point[,i,1]^2), colour = 'Low')) +
    geom_point(aes(x = 0:4, y = sqrt(bias_point[,i,1]^2 + sd_point[,i,1]^2), colour = 'Low')) +
    geom_line(aes(x = 0:4, y = sqrt(bias_point[,i,2]^2 + sd_point[,i,2]^2), colour = 'Medium')) +
    geom_point(aes(x = 0:4, y = sqrt(bias_point[,i,2]^2 + sd_point[,i,2]^2), colour = 'Medium')) +
    geom_line(aes(x = 0:4, y = sqrt(bias_point[,i,3]^2 + sd_point[,i,3]^2), colour = 'High')) +
    geom_point(aes(x = 0:4, y = sqrt(bias_point[,i,3]^2 + sd_point[,i,3]^2), colour = 'High')) +
    scale_color_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                       "High" = "green")) +
    ylim(0,0.3) +
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      mse_plots[[i]] <- mse_plots[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          mse_plots[[i]] <- mse_plots[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    mse_plots[[i]] <- mse_plots[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'Empirical root-MSE')))
  } else{mse_plots[[i]] <- mse_plots[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    mse_plots[[i]] <- mse_plots[[i]] +
      xlab('Visit')
  } else {mse_plots[[i]] <- mse_plots[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = mse_plots, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))

# ############# COVERAGE #####################
# 
# coverage_ind <- array(0,dim = c(4,5,27,3))
# success <- array(0,dim = c(4,5,27,3))
# 
# for (i in 1:1000){
#   for (k in 1:5){
#     for (j in 1:27){
#       for (l in 1:3){
#         scenario <- j%%9
#         if (scenario ==0){scenario <- 9}
# 
#         if (is.na(bootstrap[k,1,i,j,l]) == F){
#           success[1,k,j,l] <- success[1,k,j,l] + 1
#           if (all(bootstrap[k,1,i,j,l] <= true_value_red[k,scenario, l])
#               & all(bootstrap[k,2,i,j,l] >= true_value_red[k,scenario, l])){
#             coverage_ind[1,k,j,l] <- coverage_ind[1,k,j,l] + 1
#           }
#         }
# 
#         if (is.na(LEF_outcome[k,1,i,j,l]) == F){
#           success[2,k,j,l] <- success[2,k,j,l] + 1
#           if (all(LEF_outcome[k,1,i,j,l] <= true_value_red[k,scenario, l])
#               & all(LEF_outcome[k,2,i,j,l] >= true_value_red[k,scenario, l])){
#             coverage_ind[2,k,j,l] <- coverage_ind[2,k,j,l] + 1
#           }
#         }
# 
#         if (is.na(LEF_both[k,1,i,j,l]) == F){
#           success[3,k,j,l] <- success[3,k,j,l] + 1
#           if (all(LEF_both[k,1,i,j,l] <= true_value_red[k,scenario, l])
#               & all(LEF_both[k,2,i,j,l] >= true_value_red[k,scenario, l])){
#             coverage_ind[3,k,j,l] <- coverage_ind[3,k,j,l] + 1
#           }
#         }
# 
#         if (all(is.na(sandwich[k,1,i,j,l])) == F){
#           success[4,k,j,l] <- success[4,k,j,l] + 1
#           if (all(sandwich[k,1,i,j,l] <= true_value_red[k,scenario, l])
#               & all(sandwich[k,2,i,j,l] >= true_value_red[k,scenario, l])){
#             coverage_ind[4,k,j,l]<- coverage_ind[4,k,j,l]+ 1
#           }
#         }
#       }
#     }
#   }
# }
# 
# coverage_ind <- coverage_ind/success
# 
# bias_elim_coverage_ind <- array(0,dim = c(4,5,9,3))
# bias_elim_success <- array(0,dim = c(4,5,9,3))
# 
# for (i in 1:1000){
#   for (k in 1:5){
#     for (j in 1:9){
#       for (l in 1:3){
#         scenario <- j%%9
#         if (scenario ==0){scenario <- 9}
#         
#         if (is.na(bootstrap[k,1,i,j,l]) == F){
#           bias_elim_success[1,k,j,l] <- bias_elim_success[1,k,j,l] + 1
#           if (all(bootstrap[k,1,i,j,l] <= bias_point[k,j,l] + true_value_red[k,scenario, l])
#               & all(bootstrap[k,2,i,j,l] >= bias_point[k,j,l] + true_value_red[k,scenario, l])){
#             bias_elim_coverage_ind[1,k,j,l] <- bias_elim_coverage_ind[1,k,j,l] + 1
#           }
#         }
#         
#         if (is.na(LEF_outcome[k,1,i,j,l]) == F){
#           bias_elim_success[2,k,j,l] <- bias_elim_success[2,k,j,l] + 1
#           if (all(LEF_outcome[k,1,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
#               & all(LEF_outcome[k,2,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
#             bias_elim_coverage_ind[2,k,j,l] <- bias_elim_coverage_ind[2,k,j,l] + 1
#           }
#         }
#         
#         if (is.na(LEF_both[k,1,i,j,l]) == F){
#           bias_elim_success[3,k,j,l] <- bias_elim_success[3,k,j,l] + 1
#           if (all(LEF_both[k,1,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
#               & all(LEF_both[k,2,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
#             bias_elim_coverage_ind[3,k,j,l] <- bias_elim_coverage_ind[3,k,j,l] + 1
#           }
#         }
#         
#         if (all(is.na(sandwich[k,1,i,j,l])) == F){
#           bias_elim_success[4,k,j,l] <- bias_elim_success[4,k,j,l] + 1
#           if (all(sandwich[k,1,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
#               & all(sandwich[k,2,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
#             bias_elim_coverage_ind[4,k,j,l]<- bias_elim_coverage_ind[4,k,j,l]+ 1
#           }
#         }
#       }
#     }
#   }
# }
# 
# bias_elim_coverage_ind <- bias_elim_coverage_ind/bias_elim_success
# 
# 

##############
pivot_coverage_ind <- array(0,dim = c(6,5,27,3))
pivot_success <- array(0,dim = c(6,5,27,3))

for (i in 1:iters){
  for (k in 1:5){
    for (j in 1:9){
      for (l in 1:3){
        scenario <- j%%9
        if (scenario ==0){scenario <- 9}
        
        if (is.na(bootstrap[k,1,i,j,l]) == F){
          pivot_success[1,k,j,l] <- pivot_success[1,k,j,l] + 1
          if (all(2*est[k,i,j,l] - bootstrap[k,2,i,j,l] <=true_value_red[k,scenario, l])
              & all(2*est[k,i,j,l] - bootstrap[k,1,i,j,l] >= true_value_red[k,scenario, l])){
            pivot_coverage_ind[1,k,j,l] <- pivot_coverage_ind[1,k,j,l] + 1
          }
        }
        
        if (is.na(LEF_outcome[k,1,i,j,l]) == F){
          pivot_success[2,k,j,l] <- pivot_success[2,k,j,l] + 1
          if (all(2*est[k,i,j,l] - LEF_outcome[k,2,i,j,l] <= true_value_red[k,scenario, l])
              & all(2*est[k,i,j,l] - LEF_outcome[k,1,i,j,l] >=true_value_red[k,scenario, l])){
            pivot_coverage_ind[2,k,j,l] <- pivot_coverage_ind[2,k,j,l] + 1
          }
        }
        
        if (is.na(LEF_both[k,1,i,j,l]) == F){
          pivot_success[3,k,j,l] <- pivot_success[3,k,j,l] + 1
          if (all(2*est[k,i,j,l] - LEF_both[k,2,i,j,l] <= true_value_red[k,scenario, l])
              & all(2*est[k,i,j,l] - LEF_both[k,1,i,j,l] >= true_value_red[k,scenario, l])){
            pivot_coverage_ind[3,k,j,l] <- pivot_coverage_ind[3,k,j,l] + 1
          }
        }
        
        if (all(is.na(sandwich[k,1,i,j,l])) == F){
          pivot_success[4,k,j,l] <- pivot_success[4,k,j,l] + 1
          if (all(sandwich[k,1,i,j,l] <= true_value_red[k,scenario, l])
              & all(sandwich[k,2,i,j,l] >= true_value_red[k,scenario, l])){
            pivot_coverage_ind[4,k,j,l]<- pivot_coverage_ind[4,k,j,l]+ 1
          }
        }
        if (all(is.na(jackknife_mvn[k,1,i,j,l])) == F){
          pivot_success[5,k,j,l] <- pivot_success[5,k,j,l] + 1
          if (all(jackknife_mvn[k,1,i,j,l] <= true_value_red[k,scenario, l])
              & all(jackknife_mvn[k,2,i,j,l] >= true_value_red[k,scenario, l])){
            pivot_coverage_ind[5,k,j,l]<- pivot_coverage_ind[5,k,j,l]+ 1
          }
        }
        if (all(is.na(jackknife_wald[k,1,i,j,l])) == F){
          pivot_success[6,k,j,l] <- pivot_success[6,k,j,l] + 1
          if (all(jackknife_wald[k,1,i,j,l] <= true_value_red[k,scenario, l])
              & all(jackknife_wald[k,2,i,j,l] >= true_value_red[k,scenario, l])){
            pivot_coverage_ind[6,k,j,l]<- pivot_coverage_ind[6,k,j,l]+ 1
          }
        }
      }
    }
  }
}

pivot_coverage_ind <- pivot_coverage_ind/pivot_success

bias_elim_pivot_coverage_ind <- array(0,dim = c(6,5,27,3))
bias_elim_pivot_success <- array(0,dim = c(6,5,27,3))

for (i in 1:iters){
  for (k in 1:5){
    for (j in 1:27){
      for (l in 1:3){
        scenario <- j%%9
        if (scenario ==0){scenario <- 9}
        
        if (is.na(bootstrap[k,1,i,j,l]) == F){
          bias_elim_pivot_success[1,k,j,l] <- bias_elim_pivot_success[1,k,j,l] + 1
          if (all(2*est[k,i,j,l] - bootstrap[k,2,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
              & all(2*est[k,i,j,l] - bootstrap[k,1,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
            bias_elim_pivot_coverage_ind[1,k,j,l] <- bias_elim_pivot_coverage_ind[1,k,j,l] + 1
          }
        }
        
        if (is.na(LEF_outcome[k,1,i,j,l]) == F){
          bias_elim_pivot_success[2,k,j,l] <- bias_elim_pivot_success[2,k,j,l] + 1
          if (all(2*est[k,i,j,l] - LEF_outcome[k,2,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
              & all(2*est[k,i,j,l] - LEF_outcome[k,1,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
            bias_elim_pivot_coverage_ind[2,k,j,l] <- bias_elim_pivot_coverage_ind[2,k,j,l] + 1
          }
        }
        
        if (is.na(LEF_both[k,1,i,j,l]) == F){
          bias_elim_pivot_success[3,k,j,l] <- bias_elim_pivot_success[3,k,j,l] + 1
          if (all(2*est[k,i,j,l] - LEF_both[k,2,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
              & all(2*est[k,i,j,l] - LEF_both[k,1,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
            bias_elim_pivot_coverage_ind[3,k,j,l] <- bias_elim_pivot_coverage_ind[3,k,j,l] + 1
          }
        }
        
        if (all(is.na(sandwich[k,1,i,j,l])) == F){
          bias_elim_pivot_success[4,k,j,l] <- bias_elim_pivot_success[4,k,j,l] + 1
          if (all(sandwich[k,1,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
              & all(sandwich[k,2,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
            bias_elim_pivot_coverage_ind[4,k,j,l]<- bias_elim_pivot_coverage_ind[4,k,j,l]+ 1
          }
        }
        if (all(is.na(jackknife_mvn[k,1,i,j,l])) == F){
          bias_elim_pivot_success[5,k,j,l] <- bias_elim_pivot_success[5,k,j,l] + 1
          if (all(jackknife_mvn[k,1,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
              & all(jackknife_mvn[k,2,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
            bias_elim_pivot_coverage_ind[5,k,j,l]<- bias_elim_pivot_coverage_ind[5,k,j,l]+ 1
          }
        }
        if (all(is.na(jackknife_wald[k,1,i,j,l])) == F){
          bias_elim_pivot_success[6,k,j,l] <- bias_elim_pivot_success[6,k,j,l] + 1
          if (all(jackknife_wald[k,1,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
              & all(jackknife_wald[k,2,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
            bias_elim_pivot_coverage_ind[6,k,j,l]<- bias_elim_pivot_coverage_ind[6,k,j,l]+ 1
          }
        }
      }
    }
  }
}

bias_elim_pivot_coverage_ind <- bias_elim_pivot_coverage_ind/bias_elim_pivot_success


###############COVERAGE PLOTS ####################
# coverage_low <-  lapply(1:9, function(i){
#   ggplot() +
#   geom_line(aes(x = 0:4, y = coverage_ind[1,,i,1], colour = "Bootstrap")) +
#   geom_point(aes(x = 0:4, y = coverage_ind[1,,i,1], colour = "Bootstrap")) +
#   geom_line(aes(x = 0:4, y = coverage_ind[2,,i,1], colour = "LEF outcome")) +
#   geom_point(aes(x = 0:4, y = coverage_ind[2,,i,1], colour = "LEF outcome")) +
#   geom_line(aes(x = 0:4, y = coverage_ind[3,,i,1], colour = "LEF both")) +
#   geom_point(aes(x = 0:4, y = coverage_ind[3,,i,1], colour = "LEF both")) +
#   geom_line(aes(x = 0:4, y = coverage_ind[4,,i,1], colour = "Sandwich")) +
#   geom_point(aes(x = 0:4, y = coverage_ind[4,,i,1], colour = "Sandwich")) +
#   scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
#                                                   "LEF outcome" = "green", "LEF both" = "purple")) +
#   geom_hline(yintercept = 0.95, linetype = "dashed") +
#   labs(x = 'Follow up time',
#        y = "Empirical coverage rate",
#        title = paste("N =", scenarios[i,1],
#                      '\nConfounding =',scenarios[i,2],
#                      '\nTreat. prev. =', scenarios[i,3]))+ ylim(0.3,1) + 
#     theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
#   }
# )
# annotate_figure(ggarrange(plotlist = coverage_low[1:9], nrow = 3, ncol = 9, common.legend = T,
#                           legend = 'bottom'), top = 'Low event rate')
# 
# coverage_med <-  lapply(1:9, function(i){
#   ggplot() +
#     geom_line(aes(x = 0:4, y = coverage_ind[1,,i,2], colour = "Bootstrap")) +
#     geom_point(aes(x = 0:4, y = coverage_ind[1,,i,2], colour = "Bootstrap")) +
#     geom_line(aes(x = 0:4, y = coverage_ind[2,,i,2], colour = "LEF outcome")) +
#     geom_point(aes(x = 0:4, y = coverage_ind[2,,i,2], colour = "LEF outcome")) +
#     geom_line(aes(x = 0:4, y = coverage_ind[3,,i,2], colour = "LEF both")) +
#     geom_point(aes(x = 0:4, y = coverage_ind[3,,i,2], colour = "LEF both")) +
#     geom_line(aes(x = 0:4, y = coverage_ind[4,,i,2], colour = "Sandwich")) +
#     geom_point(aes(x = 0:4, y = coverage_ind[4,,i,2], colour = "Sandwich")) +
#     scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
#                                                     "LEF outcome" = "green", "LEF both" = "purple")) +
#     geom_hline(yintercept = 0.95, linetype = "dashed") +
#     labs(x = 'Follow up time',
#          y = "Empirical coverage rate",
#          title = paste("N =", scenarios[i,1],
#                        '\nConfounding =',scenarios[i,2],
#                        '\nTreat. prev. =', scenarios[i,3]))+ ylim(0.3,1) + 
#     theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
# }
# )
# annotate_figure(ggarrange(plotlist = coverage_med[1:9], nrow = 3, ncol = 9, common.legend = T,
#                           legend = 'bottom'), top = 'Medium event rate')
# 
# coverage_high <-  lapply(1:9, function(i){
#   ggplot() +
#     geom_line(aes(x = 0:4, y = coverage_ind[1,,i,3], colour = "Bootstrap")) +
#     geom_point(aes(x = 0:4, y = coverage_ind[1,,i,3], colour = "Bootstrap")) +
#     geom_line(aes(x = 0:4, y = coverage_ind[2,,i,3], colour = "LEF outcome")) +
#     geom_point(aes(x = 0:4, y = coverage_ind[2,,i,3], colour = "LEF outcome")) +
#     geom_line(aes(x = 0:4, y = coverage_ind[3,,i,3], colour = "LEF both")) +
#     geom_point(aes(x = 0:4, y = coverage_ind[3,,i,3], colour = "LEF both")) +
#     geom_line(aes(x = 0:4, y = coverage_ind[4,,i,3], colour = "Sandwich")) +
#     geom_point(aes(x = 0:4, y = coverage_ind[4,,i,3], colour = "Sandwich")) +
#     scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
#                                                     "LEF outcome" = "green", "LEF both" = "purple")) +
#     geom_hline(yintercept = 0.95, linetype = "dashed") +
#     labs(x = 'Follow up time',
#          y = "Empirical coverage rate",
#          title = paste("N =", scenarios[i,1],
#                        '\nConfounding =',scenarios[i,2],
#                        '\nTreat. prev. =', scenarios[i,3]))+ ylim(0.3,1) + 
#     theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
# }
# )
# annotate_figure(ggarrange(plotlist = coverage_high[1:9], nrow = 3, ncol = 9, common.legend = T,
#                           legend = 'bottom'), top = 'High event rate')
# 
# ######################BIAS ELIM COVERAGE PLOTS #################
# coverage_low <-  lapply(1:9, function(i){
#   ggplot() +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[1,,i,1], colour = "Bootstrap")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[1,,i,1], colour = "Bootstrap")) +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[2,,i,1], colour = "LEF outcome")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[2,,i,1], colour = "LEF outcome")) +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[3,,i,1], colour = "LEF both")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[3,,i,1], colour = "LEF both")) +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[4,,i,1], colour = "Sandwich")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[4,,i,1], colour = "Sandwich")) +
#     scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
#                                                     "LEF outcome" = "green", "LEF both" = "purple")) +
#     geom_hline(yintercept = 0.95, linetype = "dashed") +
#     labs(x = 'Follow up time',
#          y = "Bias-eliminated coverage",
#          title = paste("N =", scenarios[i,1],
#                        '\nConfounding =',scenarios[i,2],
#                        '\nTreat. prev. =', scenarios[i,3]))+ ylim(0.3,1) + 
#     theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
# }
# )
# annotate_figure(ggarrange(plotlist = coverage_low[1:9], nrow = 3, ncol = 9, common.legend = T,
#                           legend = 'bottom'), top = 'Low event rate')
# 
# coverage_med <-  lapply(1:9, function(i){
#   ggplot() +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[1,,i,2], colour = "Bootstrap")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[1,,i,2], colour = "Bootstrap")) +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[2,,i,2], colour = "LEF outcome")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[2,,i,2], colour = "LEF outcome")) +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[3,,i,2], colour = "LEF both")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[3,,i,2], colour = "LEF both")) +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[4,,i,2], colour = "Sandwich")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[4,,i,2], colour = "Sandwich")) +
#     scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
#                                                     "LEF outcome" = "green", "LEF both" = "purple")) +
#     geom_hline(yintercept = 0.95, linetype = "dashed") +
#     labs(x = 'Follow up time',
#          y = "Bias-eliminated coverage",
#          title = paste("N =", scenarios[i,1],
#                        '\nConfounding =',scenarios[i,2],
#                        '\nTreat. prev. =', scenarios[i,3]))+ ylim(0.3,1) + 
#     theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
# }
# )
# annotate_figure(ggarrange(plotlist = coverage_med[1:9], nrow = 3, ncol = 9, common.legend = T,
#                           legend = 'bottom'), top = 'Medium event rate')
# 
# coverage_high <-  lapply(1:9, function(i){
#   ggplot() +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[1,,i,3], colour = "Bootstrap")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[1,,i,3], colour = "Bootstrap")) +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[2,,i,3], colour = "LEF outcome")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[2,,i,3], colour = "LEF outcome")) +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[3,,i,3], colour = "LEF both")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[3,,i,3], colour = "LEF both")) +
#     geom_line(aes(x = 0:4, y = bias_elim_coverage_ind[4,,i,3], colour = "Sandwich")) +
#     geom_point(aes(x = 0:4, y = bias_elim_coverage_ind[4,,i,3], colour = "Sandwich")) +
#     scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
#                                                     "LEF outcome" = "green", "LEF both" = "purple")) +
#     geom_hline(yintercept = 0.95, linetype = "dashed") +
#     labs(x = 'Follow up time',
#          y = "Bias-eliminated coverage",
#          title = paste("N =", scenarios[i,1],
#                        '\nConfounding =',scenarios[i,2],
#                        '\nTreat. prev. =', scenarios[i,3]))+ ylim(0.3,1) + 
#     theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
# }
# )
# annotate_figure(ggarrange(plotlist = coverage_high[1:9], nrow = 3, ncol = 9, common.legend = T,
#                           legend = 'bottom'), top = 'High event rate')
# ############### PIVOT COVERAGE #######################
coverage_low <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[1,,i,1], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[1,,i,1], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[2,,i,1], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[2,,i,1], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[3,,i,1], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[3,,i,1], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[4,,i,1], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[4,,i,1], colour = "Sandwich")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[5,,i,1], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[5,,i,1], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[6,,i,1], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[6,,i,1], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    ylim(0.4,1) + 
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      coverage_low[[i]] <- coverage_low[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          coverage_low[[i]] <- coverage_low[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
    }
  }
  if(i %in% c(1,10,19)){
    coverage_low[[i]] <- coverage_low[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'Coverage')))
  } else{coverage_low[[i]] <- coverage_low[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    coverage_low[[i]] <- coverage_low[[i]] +
      xlab('Visit')
  } else {coverage_low[[i]] <- coverage_low[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = coverage_low, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))

coverage_med <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[1,,i,2], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[1,,i,2], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[2,,i,2], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[2,,i,2], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[3,,i,2], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[3,,i,2], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[4,,i,2], colour = "Sandwich")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[5,,i,2], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[5,,i,2], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[6,,i,2], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[6,,i,2], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    ylim(0.4,1) + 
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      coverage_med[[i]] <- coverage_med[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          coverage_med[[i]] <- coverage_med[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    coverage_med[[i]] <- coverage_med[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'Coverage')))
  } else{coverage_med[[i]] <- coverage_med[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    coverage_med[[i]] <- coverage_med[[i]] +
      xlab('Visit')
  } else {coverage_med[[i]] <- coverage_med[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = coverage_med, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))

coverage_high <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[4,,i,3], colour = "Sandwich")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[5,,i,3], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[5,,i,3], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[6,,i,3], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[6,,i,3], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    ylim(0.4,1) + 
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      coverage_high[[i]] <- coverage_high[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          coverage_high[[i]] <- coverage_high[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    coverage_high[[i]] <- coverage_high[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'Coverage')))
  } else{coverage_high[[i]] <- coverage_high[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    coverage_high[[i]] <- coverage_high[[i]] +
      xlab('Visit')
  } else {coverage_high[[i]] <- coverage_high[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = coverage_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))

######################   BIAS ELIM PIVOT COVERAGE ########################
coverage_low <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,1], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,1], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,1], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,1], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,1], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,1], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,1], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,1], colour = "Sandwich")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[5,,i,1], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[5,,i,1], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[6,,i,1], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[6,,i,1], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    ylim(0.4,1) + 
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      coverage_low[[i]] <- coverage_low[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          coverage_low[[i]] <- coverage_low[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    coverage_low[[i]] <- coverage_low[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'Coverage')))
  } else{coverage_low[[i]] <- coverage_low[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    coverage_low[[i]] <- coverage_low[[i]] +
      xlab('Visit')
  } else {coverage_low[[i]] <- coverage_low[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = coverage_low, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))

coverage_med <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,2], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,2], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,2], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,2], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,2], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,2], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,2], colour = "Sandwich")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[5,,i,2], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[5,,i,2], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[6,,i,2], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[6,,i,2], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    ylim(0.4,1) + 
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      coverage_med[[i]] <- coverage_med[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          coverage_med[[i]] <- coverage_med[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    coverage_med[[i]] <- coverage_med[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'Coverage')))
  } else{coverage_med[[i]] <- coverage_med[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    coverage_med[[i]] <- coverage_med[[i]] +
      xlab('Visit')
  } else {coverage_med[[i]] <- coverage_med[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = coverage_med, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))

coverage_high <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,3], colour = "Sandwich")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[5,,i,3], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[5,,i,3], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[6,,i,3], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[6,,i,3], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    ylim(0.4,1) + 
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      coverage_high[[i]] <- coverage_high[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          coverage_high[[i]] <- coverage_high[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    coverage_high[[i]] <- coverage_high[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'Coverage')))
  } else{coverage_high[[i]] <- coverage_high[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    coverage_high[[i]] <- coverage_high[[i]] +
      xlab('Visit')
  } else {coverage_high[[i]] <- coverage_high[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = coverage_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))

###################### CI WIDTH ##############################

CI_width <- array(, dim = c(6,5,27,3))

CI_width <- data.frame(matrix(,nrow = 0, ncol = 7))
for (j in 1:27){
  for (l in 1:3){
    CI_width <-rbind(CI_width,cbind(0:4, 'Bootstrap', rowMeans(abs(bootstrap[,2,,j,l] - bootstrap[,1,,j,l]), na.rm = TRUE), 
                                    rowMeans(abs(bootstrap[,2,,j,l] - bootstrap[,1,,j,l]), na.rm = TRUE) - rowSds(abs(bootstrap[,2,,j,l] - bootstrap[,1,,j,l]), na.rm = TRUE),
                                    rowMeans(abs(bootstrap[,2,,j,l] - bootstrap[,1,,j,l]), na.rm = TRUE) + rowSds(abs(bootstrap[,2,,j,l] - bootstrap[,1,,j,l]), na.rm = TRUE),
                                    j,l))
    CI_width <-rbind(CI_width,cbind(0:4, 'LEF outcome', rowMeans(abs(LEF_outcome[,2,,j,l] - LEF_outcome[,1,,j,l]), na.rm = TRUE), 
                                    rowMeans(abs(LEF_outcome[,2,,j,l] - LEF_outcome[,1,,j,l]), na.rm = TRUE) - rowSds(abs(LEF_outcome[,2,,j,l] - LEF_outcome[,1,,j,l]), na.rm = TRUE),
                                    rowMeans(abs(LEF_outcome[,2,,j,l] - LEF_outcome[,1,,j,l]), na.rm = TRUE) + rowSds(abs(LEF_outcome[,2,,j,l] - LEF_outcome[,1,,j,l]), na.rm = TRUE),
                                    l,j))
    CI_width <-rbind(CI_width,cbind(0:4, 'LEF both', rowMeans(abs(LEF_both[,2,,j,l] - LEF_both[,1,,j,l]), na.rm = TRUE), 
                                    rowMeans(abs(LEF_both[,2,,j,l] - LEF_both[,1,,j,l]), na.rm = TRUE) - rowSds(abs(LEF_both[,2,,j,l] - LEF_both[,1,,j,l]), na.rm = TRUE),
                                    rowMeans(abs(LEF_both[,2,,j,l] - LEF_both[,1,,j,l]), na.rm = TRUE) + rowSds(abs(LEF_both[,2,,j,l] - LEF_both[,1,,j,l]), na.rm = TRUE),
                                    l,j))
    CI_width <-rbind(CI_width,cbind(0:4, 'Sandwich', rowMeans(abs(sandwich[,2,,j,l] - sandwich[,1,,j,l]), na.rm = TRUE), 
                                    rowMeans(abs(sandwich[,2,,j,l] - sandwich[,1,,j,l]), na.rm = TRUE) - rowSds(abs(sandwich[,2,,j,l] - sandwich[,1,,j,l]), na.rm = TRUE),
                                    rowMeans(abs(sandwich[,2,,j,l] - sandwich[,1,,j,l]), na.rm = TRUE) + rowSds(abs(sandwich[,2,,j,l] - sandwich[,1,,j,l]), na.rm = TRUE),
                                    l,j))
    CI_width <-rbind(CI_width,cbind(0:4, 'Jackknife MVN', rowMeans(abs(jackknife_mvn[,2,,j,l] - jackknife_mvn[,1,,j,l]), na.rm = TRUE), 
                                    rowMeans(abs(jackknife_mvn[,2,,j,l] - jackknife_mvn[,1,,j,l]), na.rm = TRUE) - rowSds(abs(jackknife_mvn[,2,,j,l] - jackknife_mvn[,1,,j,l]), na.rm = TRUE),
                                    rowMeans(abs(jackknife_mvn[,2,,j,l] - jackknife_mvn[,1,,j,l]), na.rm = TRUE) + rowSds(abs(jackknife_mvn[,2,,j,l] - jackknife_mvn[,1,,j,l]), na.rm = TRUE),
                                    l,j))
    CI_width <-rbind(CI_width,cbind(0:4, 'Jackknife Wald', rowMeans(abs(jackknife_wald[,2,,j,l] - jackknife_wald[,1,,j,l]), na.rm = TRUE), 
                                    rowMeans(abs(jackknife_wald[,2,,j,l] - jackknife_wald[,1,,j,l]), na.rm = TRUE) - rowSds(abs(jackknife_wald[,2,,j,l] - jackknife_wald[,1,,j,l]), na.rm = TRUE),
                                    rowMeans(abs(jackknife_wald[,2,,j,l] - jackknife_wald[,1,,j,l]), na.rm = TRUE) + rowSds(abs(jackknife_wald[,2,,j,l] - jackknife_wald[,1,,j,l]), na.rm = TRUE),
                                    l,j))
  }}

colnames(CI_width) <- c('Visit', 'CI_type', "Mean", "Min", "Max", 'i', 'j')
CI_width$Visit <- as.numeric(CI_width$Visit)
CI_width$Mean <- as.numeric(CI_width$Mean)
CI_width$Min <- as.numeric(CI_width$Min)
CI_width$Max <- as.numeric(CI_width$Max)
CI_width$i <- as.numeric(CI_width$i)
CI_width$j <- as.numeric(CI_width$j)


width_low <-  lapply(1:27, function(i){
    ggplot(CI_width[CI_width$i == i & CI_width$j == 1,],aes( x = Visit, y = Mean, group = CI_type, colour = CI_type)) +
    geom_point(position=position_dodge(width=0.5), na.rm = TRUE) +    
    geom_errorbar(aes(x = Visit, ymin = Min, ymax = Max),position=position_dodge(width=0.5), na.rm = TRUE) +
  
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("CI Width") +  ylim(0,1.5) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))
}
)
annotate_figure(ggarrange(plotlist = width_low, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

width_med <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width[1,,i,2], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width[1,,i,2], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width[2,,i,2], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width[2,,i,2], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width[3,,i,2], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width[3,,i,2], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width[4,,i,2], colour = "Sandwich")) + 
    geom_line(aes(x = 0:4, y = CI_width[5,,i,2], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = CI_width[5,,i,2], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = CI_width[6,,i,2], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = CI_width[6,,i,2], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Average Width") +  ylim(0,1.5) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_med, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

width_high <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width[4,,i,3], colour = "Sandwich")) +
    geom_line(aes(x = 0:4, y = CI_width[5,,i,3], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = CI_width[5,,i,3], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = CI_width[6,,i,3], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = CI_width[6,,i,3], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Average Width") +  ylim(0,1.5) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

############################# CI WIDTH SD###############################
CI_width_sd <- array(, dim = c(4,5,9,3))
for (k in 1:5){
  for (j in 1:27){
    for (l in 1:3){
      CI_width_sd[1,k,j,l] <- sd(abs(bootstrap[k,2,,j,l] - bootstrap[k,1,,j,l]), na.rm = TRUE)
      CI_width_sd[2,k,j,l] <- sd(abs(LEF_outcome[k,2,,j,l] - LEF_outcome[k,1,,j,l]), na.rm = TRUE)
      CI_width_sd[3,k,j,l] <- sd(abs(LEF_both[k,2,,j,l] - LEF_both[k,1,,j,l]), na.rm = TRUE)
      CI_width_sd[4,k,j,l] <- sd(abs(sandwich[k,2,,j,l] - sandwich[k,1,,j,l]), na.rm = TRUE)
    }}}

width_low <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width_sd[1,,i,1], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[1,,i,1], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width_sd[2,,i,1], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[2,,i,1], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width_sd[3,,i,1], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[3,,i,1], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width_sd[4,,i,1], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[4,,i,1], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width SD") +  ylim(0,0.7) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))
}
)
annotate_figure(ggarrange(plotlist = width_low, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))


width_med <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width_sd[1,,i,2], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[1,,i,2], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width_sd[2,,i,2], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[2,,i,2], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width_sd[3,,i,2], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[3,,i,2], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width_sd[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[4,,i,2], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width SD") +  ylim(0,0.7) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_med, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

width_high <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width_sd[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width_sd[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width_sd[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width_sd[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width_sd[4,,i,3], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width SD") +  ylim(0,0.7) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

########################################### CI WIDTH MEDIAN ##########################
CI_width_median <- array(, dim = c(4,5,9,3))
for (k in 1:5){
  for (j in 1:27){
    for (l in 1:3){
      CI_width_median[1,k,j,l] <- median(abs(bootstrap[k,2,,j,l] - bootstrap[k,1,,j,l]), na.rm = TRUE)
      CI_width_median[2,k,j,l] <- median(abs(LEF_outcome[k,2,,j,l] - LEF_outcome[k,1,,j,l]), na.rm = TRUE)
      CI_width_median[3,k,j,l] <- median(abs(LEF_both[k,2,,j,l] - LEF_both[k,1,,j,l]), na.rm = TRUE)
      CI_width_median[4,k,j,l] <- median(abs(sandwich[k,2,,j,l] - sandwich[k,1,,j,l]), na.rm = TRUE)
    }}}

width_low <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width_median[1,,i,1], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width_median[1,,i,1], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width_median[2,,i,1], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width_median[2,,i,1], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width_median[3,,i,1], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width_median[3,,i,1], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width_median[4,,i,1], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width_median[4,,i,1], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width median") +  ylim(0,0.8) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))
}
)
annotate_figure(ggarrange(plotlist = width_low, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))


width_med <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width_median[1,,i,2], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width_median[1,,i,2], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width_median[2,,i,2], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width_median[2,,i,2], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width_median[3,,i,2], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width_median[3,,i,2], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width_median[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width_median[4,,i,2], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width median") +  ylim(0,0.8) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_med, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

width_high <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width_median[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width_median[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width_median[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width_median[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width_median[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width_median[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width_median[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width_median[4,,i,3], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    xlab(bquote(atop(n ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width median") +  ylim(0,0.8) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

######################PIVOT FAILURE RATE #######################
failure_table <- data.frame(matrix(,nrow = 0, ncol = 8))
for (i in 1:27){
  for (j in 1:3){
    failure_table <- rbind(failure_table, c(outcomes[j], scenarios[i,1], scenarios[i,2],scenarios[i,3], 
                                            (iters - pivot_success[,1,i,j])/iters))
  }
}
colnames(failure_table) <- c('Outcome_prevalence', 'Sample_size', 'Confounding', 'Treatment_prevalence', 'Bootstrap', 'LEF_outcome',
                             'LEF_both', 'Sandwich', 'Jackknife_MVN', 'Jackknife_Wald')
failure_table$Sample_size <- as.numeric(failure_table$Sample_size)
failure_table$Confounding <- as.numeric(failure_table$Confounding)
failure_table$Treatment_prevalence <- as.numeric(failure_table$Treatment_prevalence)
failure_table$Bootstrap <- as.numeric(failure_table$Bootstrap)
failure_table$Jackknife_MVN <- as.numeric(failure_table$Jackknife_MVN)
failure_table$LEF_outcome <- as.numeric(failure_table$LEF_outcome)
failure_table$LEF_both <- as.numeric(failure_table$LEF_both)
failure_table$Sandwich <- as.numeric(failure_table$Sandwich)
failure_table$Jackknife_Wald <- as.numeric(failure_table$Jackknife_Wald)
failure_table <- failure_table %>% 
  dplyr::group_by(Outcome_prevalence, Sample_size, Treatment_prevalence) %>% 
  dplyr::summarise(Bootstrap = mean(as.numeric(Bootstrap)),
                   LEF_outcome = mean(as.numeric(LEF_outcome)),
                   LEF_both = mean(as.numeric(LEF_both)),
                   Sandwich = mean(as.numeric(Sandwich)),
                   Jackknife_MVN = ifelse(Sample_size == 200, mean(as.numeric(Jackknife_MVN)), NA),
                   Jackknife_Wald = ifelse(Sample_size == 200, mean(as.numeric(Jackknife_Wald)), NA))
print(xtable(failure_table, type = 'latex'), include.rownames = FALSE) 
no_na_frequency <- data.frame(matrix(,nrow = 0, ncol = 7))
for (i in 1:27){
  for (j in 1:3){
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ/coeff_dim_',outcomes[j],'_',i, '.rda'))
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ/bootstrap_nas_',outcomes[j],'_',i, '.rda'))
    load(paste0('~/rds/hpc-work/Project1/NewSimusJ/jackknife_nas_',outcomes[j],'_',i, '.rda'))
    no_na_frequency <- rbind(no_na_frequency, cbind(outcomes[j], scenarios[i,1], scenarios[i,2],scenarios[i,3],
                                            coeff_dim, bootstrap_nas, jackknife_nas))
  }
}
colnames(no_na_frequency) <- c('Outcome_prevalence', 'Sample_size', 'Confounding', 'Treatment_prevalence', 'Coeff_dim', 'Bootstrap_NAs',
                             'Jackknife_NAs')

no_na_frequency$Sample_size <- as.numeric(no_na_frequency$Sample_size)
no_na_frequency$Confounding <- as.numeric(no_na_frequency$Confounding)
no_na_frequency$Treatment_prevalence <- as.numeric(no_na_frequency$Treatment_prevalence)

no_na_frequency$Coeff_dim <- as.numeric(no_na_frequency$Coeff_dim)
no_na_frequency$Bootstrap_NAs <- as.numeric(no_na_frequency$Bootstrap_NAs)
no_na_frequency$Jackknife_NAs <- as.numeric(no_na_frequency$Jackknife_NAs)


frenquency_contingency <- no_na_frequency %>% 
  dplyr::mutate(MSM_coeffs = ifelse(Coeff_dim == 20, 'No NAs', 'Some NAs'),
                Bootstrap = ifelse(Bootstrap_NAs == 0, 'No NAs', 'Some NAs'),
                Jackknife = ifelse(Jackknife_NAs == 0, 'No NAs', 'Some NAs'))

con4 <- xtabs(~MSM_coeffs+Jackknife+ Outcome_prevalence + Confounding + Treatment_prevalence, 
              data = frenquency_contingency[frenquency_contingency$Sample_size == 200,])
ftable(con4, row.vars = c('Outcome_prevalence', 'Confounding', 'Treatment_prevalence', 'MSM_coeffs'))

con4 <- xtabs(~MSM_coeffs+Bootstrap+ Outcome_prevalence + Confounding + Treatment_prevalence, 
              data = frenquency_contingency[frenquency_contingency$Sample_size == 200,])
ftable(con4, row.vars = c('Outcome_prevalence', 'Confounding', 'Treatment_prevalence', 'MSM_coeffs'))
 ################ MC SE PLOTS #######################
MC_SE_pivot<- array(, dim = c(6,5,27,3))
for (i in 1:27){
  for (j in 1:3){
    MC_SE_pivot[,,i,j] <- sqrt((pivot_coverage_ind[,,i,j]*(1-pivot_coverage_ind[,,i,j]))/iters)
  }
}


MCSE_low <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[1,,i,1], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[1,,i,1], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[2,,i,1], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[2,,i,1], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[3,,i,1], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[3,,i,1], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[4,,i,1], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[4,,i,1], colour = "Sandwich")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[5,,i,1], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[5,,i,1], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[6,,i,1], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[6,,i,1], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    ylim(0,0.02) + 
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      MCSE_low[[i]] <- MCSE_low[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          MCSE_low[[i]] <- MCSE_low[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    MCSE_low[[i]] <- MCSE_low[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'MC SE')))
  } else{MCSE_low[[i]] <- MCSE_low[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    MCSE_low[[i]] <- MCSE_low[[i]] +
      xlab('Visit')
  } else {MCSE_low[[i]] <- MCSE_low[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = MCSE_low, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))

MCSE_med <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[1,,i,2], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[1,,i,2], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[2,,i,2], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[2,,i,2], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[3,,i,2], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[3,,i,2], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[4,,i,2], colour = "Sandwich")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[5,,i,2], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[5,,i,2], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[6,,i,2], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[6,,i,2], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    ylim(0,0.02) + 
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      MCSE_med[[i]] <- MCSE_med[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          MCSE_med[[i]] <- MCSE_med[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    MCSE_med[[i]] <- MCSE_med[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'MC SE')))
  } else{MCSE_med[[i]] <- MCSE_med[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    MCSE_med[[i]] <- MCSE_med[[i]] +
      xlab('Visit')
  } else {MCSE_med[[i]] <- MCSE_med[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = MCSE_med, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))

MCSE_high <-  lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[4,,i,3], colour = "Sandwich")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[5,,i,3], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[5,,i,3], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[6,,i,3], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[6,,i,3], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'deepskyblue' )) +
    ylim(0,0.02) + 
    theme(legend.text = element_text(size=14))
}
)
for(i in 1:27){
  if(i %in% 1:9){
    if(i %in% c(2, 5, 8)){
      MCSE_high[[i]] <- MCSE_high[[i]] +
        labs(title = bquote(atop(alpha[c] == .(scenarios[i,2]), alpha[a] == .(scenarios[i,3]))))} else{
          MCSE_high[[i]] <- MCSE_high[[i]] +
            labs(title = bquote(atop(phantom(3),alpha[a] == .(scenarios[i,3]))))
        }
  }
  if(i %in% c(1,10,19)){
    MCSE_high[[i]] <- MCSE_high[[i]] +
      ylab(bquote(atop(n == .(scenarios[i,1]), 'MC SE')))
  } else{MCSE_high[[i]] <- MCSE_high[[i]] +
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())
  }
  if(i %in% 19:27){
    MCSE_high[[i]] <- MCSE_high[[i]] +
      xlab('Visit')
  } else {MCSE_high[[i]] <- MCSE_high[[i]] +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank())
  }
}

annotate_figure(ggarrange(plotlist = MCSE_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom',
                          widths = c(1.4,1,1,1,1,1,1,1,1),
                          heights = c(1.1, 0.95, 1)))
################WEIGHTS########################
weights <- data.frame(matrix(,nrow = 0, ncol = 7))
for(i in 1:27){
  simdata_censored<-DATA_GEN_censored_reduced(as.numeric(scenarios[i,1]), 5, 
                                              conf = as.numeric(scenarios[i,2]), 
                                              treat_prev = as.numeric(scenarios[i,3]),
                                              outcome_prev = -4.7,
                                              censor = F)
  PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                              eligible ='eligible',
                                              switch_d_cov = ~X2 + X4,
                                              outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                              use_weight=T, use_censor=T, quiet = T,
                                              save_weight_models = F,
                                              data_dir = data_direction)
  
  weights <- rbind(weights, c(-4.7, as.numeric(scenarios[i,1]), as.numeric(scenarios[i,3]),as.numeric(scenarios[i,2]),
                              mean(PP_prep$data$weight), sd(PP_prep$data$weight), max(PP_prep$data$weight)))
}

for(i in 1:27){
  simdata_censored<-DATA_GEN_censored_reduced(as.numeric(scenarios[i,1]), 5, 
                                              conf = as.numeric(scenarios[i,2]), 
                                              treat_prev = as.numeric(scenarios[i,3]),
                                              outcome_prev = -3.8,
                                              censor = F)
  PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                              eligible ='eligible',
                                              switch_d_cov = ~X2 + X4,
                                              outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                              use_weight=T, use_censor=T, quiet = T,
                                              save_weight_models = F,
                                              data_dir = data_direction)
  weights <- rbind(weights, c(-3.8, as.numeric(scenarios[i,1]),as.numeric(scenarios[i,3]),as.numeric(scenarios[i,2]),
                              mean(PP_prep$data$weight), sd(PP_prep$data$weight), max(PP_prep$data$weight)))
}

for(i in 1:27){
  simdata_censored<-DATA_GEN_censored_reduced(as.numeric(scenarios[i,1]), 5, 
                                              conf = as.numeric(scenarios[i,2]), 
                                              treat_prev = as.numeric(scenarios[i,3]),
                                              outcome_prev = -3,
                                              censor = F)
  PP_prep <- TrialEmulation::data_preparation(simdata_censored, id='ID', period='t', treatment='A', outcome='Y', 
                                              eligible ='eligible',
                                              switch_d_cov = ~X2 + X4,
                                              outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                              use_weight=T, use_censor=T, quiet = T,
                                              save_weight_models = F,
                                              data_dir = data_direction)
  weights <- rbind(weights, c(-3, as.numeric(scenarios[i,1]), as.numeric(scenarios[i,3]),as.numeric(scenarios[i,2]),
                              mean(PP_prep$data$weight), sd(PP_prep$data$weight), max(PP_prep$data$weight)))
}

colnames(weights) <- c('Outcome_prevalence','Sample_size', 'Treatment_prevalence', 'Confounding', 'Mean', 'SD',
                       'Max')
weights$Mean <- as.numeric(weights$Mean)
weights$SD <- as.numeric(weights$SD)
weights$Max <- as.numeric(weights$Max)

print(xtable(weights %>% filter(Sample_size == 5000) %>% select(-Sample_size) %>% arrange(Outcome_prevalence,Treatment_prevalence, Confounding ), 
             type = 'latex', digits = c(0,2,2,2,2,2,2)),include.rownames=FALSE)



##################### DATA TABULATION ###################
l <- 1
data_med <-DATA_GEN_censored_reduced(100000, 5, 
                                                              conf = as.numeric(scenarios[l,2]), 
                                                              treat_prev = as.numeric(scenarios[l,3]),
                                                              outcome_prev = -4.7,
                                                              censor = F)
PP_prep <- TrialEmulation::data_preparation(data_med, id='ID', period='t', treatment='A', outcome='Y',cense = 'C', 
                                            eligible ='eligible',
                                            switch_d_cov = ~X2 + X4,
                                            outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                            estimand_type = 'PP', quiet = F,
                                            save_weight_models = F,
                                            data_dir = data_direction)
switch_data <- PP_prep$data %>% 
  dplyr::mutate(t_1 = ifelse(followup_time == 1,1,0),
                t_2 = ifelse(followup_time == 2,1,0),
                t_3 = ifelse(followup_time == 3,1,0),
                t_4 = ifelse(followup_time == 4,1,0),
                t_1A = t_1*assigned_treatment,
                t_2A = t_2*assigned_treatment,
                t_3A = t_3*assigned_treatment,
                t_4A = t_4*assigned_treatment,
                t_1X2 = t_1*X2,
                t_2X2 = t_2*X2,
                t_3X2 = t_3*X2,
                t_4X2 = t_4*X2,
                t_1X4 = t_1*X4,
                t_2X4 = t_2*X4,
                t_3X4 = t_3*X4,
                t_4X4 = t_4*X4)

my_covariates <- ~ X2 + X4+ assigned_treatment+
  t_1 + t_2 + t_3 + t_4 +
  t_1A + t_2A + t_3A + t_4A + 
  t_1X2 + t_2X2 + t_3X2 + t_4X2 + 
  t_1X4 + t_2X4 + t_3X4 + t_4X4

PP <- TrialEmulation::trial_msm(data = switch_data,
                                outcome_cov = my_covariates,
                                model_var = c('assigned_treatment'),
                                glm_function = 'glm',
                                include_trial_period = ~1, include_followup_time = ~1,
                                estimand_type = "PP", quiet = F, use_sample_weights =  F)
big_data_med <- PP_prep$data %>% 
  dplyr::mutate(sample_size = as.numeric(scenarios[l,1]), conf = as.numeric(scenarios[l,2]), treat = as.numeric(scenarios[l,3]))
for (l in 20:27){
  data_med <-DATA_GEN_censored_reduced(as.numeric(scenarios[l,1]), 5, 
                                       conf = as.numeric(scenarios[l,2]), 
                                       treat_prev = as.numeric(scenarios[l,3]),
                                       outcome_prev = -3.0,
                                       censor = F)
  PP_prep <- TrialEmulation::data_preparation(data_med, id='ID', period='t', treatment='A', outcome='Y',cense = 'C', 
                                              eligible ='eligible',
                                              switch_d_cov = ~X2 + X4,
                                              outcome_cov = ~X2 + X4, model_var = c('assigned_treatment'),
                                              estimand_type = 'PP', quiet = T,
                                              save_weight_models = F,
                                              data_dir = data_direction)
  big_data_med_loop <- PP_prep$data %>% 
    dplyr::mutate(sample_size = as.numeric(scenarios[l,1]), conf = as.numeric(scenarios[l,2]), treat = as.numeric(scenarios[l,3]))
  
  big_data_med <- rbind(big_data_med, big_data_med_loop)
}  

con4<-xtabs(~sample_size+conf+treat+assigned_treatment + outcome + followup_time, data=big_data_med)
ftable(con4)

xftbl <- xtableFtable(ftable(con4), method = "compact")
print.xtableFtable(xftbl, booktabs = T) 

############ Table of SE for n = 5000 ###############


se_5000_table <- rbind(cbind(0:4,sd_point[,19:27,1]),
                       cbind(0:4,sd_point[,19:27,2]),
                       cbind(0:4,sd_point[,19:27,3]))

print(xtable(se_5000_table,
             type = 'latex', digits = c(0,0,3,3,3,3,3,3,3,3,3)),include.rownames=FALSE)
