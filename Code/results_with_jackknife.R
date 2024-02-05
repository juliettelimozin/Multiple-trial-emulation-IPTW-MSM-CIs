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

iters <- 1000

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


size <- c(200,1000,5000)
treat <- c(-1,0,1)
conf <- c(0.1,0.5,0.9)

scenarios <- as.data.frame(tidyr::crossing(size,conf, treat))
bias_point <- array(,dim = c(5,27,3))
sd_point <- array(,dim = c(5,27,3))
mean_time <- data.frame(matrix(,nrow = 0, ncol = 10))
true_value_red <- -true_value_red
na_failure_rate <- data.frame(matrix(,nrow = 0, ncol = 7))
for (i in 1:27){
  for (j in 1:3){
    load(paste0("HPC output/NewSimus/J_CI_bootstrap_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("HPC output/NewSimus/J_CI_jackknife_mvn_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("HPC output/NewSimus/J_CI_jackknife_wald_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("HPC output/NewSimus/J_CI_LEF_outcome_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("HPC output/NewSimus/J_CI_LEF_both_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("HPC output/NewSimus/J_CI_sandwich_PP_red_",outcomes[j],'_', i, ".rda"))
    load(paste0("HPC output/NewSimus/J_computation_time_",outcomes[j],'_', i, ".rda"))
    load(paste0('HPC output/NewSimus/J_estimates_red_',outcomes[j],'_',i, '.rda'))
    load(paste0('HPC output/NewSimus/coeff_dim_',outcomes[j],'_',i, '.rda'))
    load(paste0('HPC output/NewSimus/bootstrap_nas_',outcomes[j],'_',i, '.rda'))
    load(paste0('HPC output/NewSimus/jackknife_nas_',outcomes[j],'_',i, '.rda'))
    
    bootstrap[,,,i,j] <- CI_bootstrap_PP_red
    LEF_outcome[,,,i,j] <- CI_LEF_outcome_PP_red
    LEF_both[,,,i,j] <- CI_LEF_both_PP_red
    sandwich[,,,i,j] <- CI_sandwich_PP_red
    jackknife_mvn[,,,i,j] <- CI_jackknife_mvn_PP_red
    jackknife_wald[,,,i,j] <- CI_jackknife_wald_PP_red
    
    na_failure_rate <- rbind(na_failure_rate, c(outcomes[j], scenarios[i,1], scenarios[i,2],scenarios[i,3], 
                                          mean(coeff_dim, na.rm = TRUE),
                                          mean(bootstrap_nas, na.rm = TRUE),
                                          mean(jackknife_nas, na.rm = TRUE)
                                          ))
    
    mean_time <- rbind(mean_time, c(outcomes[j], scenarios[i,1], scenarios[i,2],scenarios[i,3], rowMeans(computation_time, na.rm = TRUE)))
    scenario <- i%%9
    est[,,i,j] <- estimates
    if (scenario ==0){scenario <- 9}
    bias_point[,i,j] <- rowMeans(estimates, na.rm = TRUE) - true_value_red[,scenario,j]
    sd_point[,i,j] <- rowSds(estimates, na.rm = TRUE)
  }
}

colnames(mean_time) <- c('Outcome_prevalence', 'Sample_size', 'Confounding', 'Treatment_prevalence', 'Bootstrap', 'LEF_outcome',
                         'LEF_both', 'Sandwich')
mean_time <- mean_time %>% 
  dplyr::mutate(Sample_size = as.numeric(Sample_size),
                Confounding = as.numeric(Confounding),
                Treatment_prevalence = as.numeric(Treatment_prevalence),
                Bootstrap = as.numeric(Bootstrap)/as.numeric(Sandwich),
                LEF_outcome = as.numeric(LEF_outcome)/as.numeric(Sandwich),
                LEF_both = as.numeric(LEF_both)/as.numeric(Sandwich),
                Sandwich = 1) %>% 
  dplyr::group_by(Outcome_prevalence,Sample_size) %>% 
  dplyr::summarise(Bootstrap = mean(Bootstrap),
                   LEF_outcome = mean(LEF_outcome),
                   LEF_both = mean(LEF_both),
                   sandwich = 1)

mrd_dist1_low <- lapply(1:9, function(i){
  scenario <- i%%9
  if (scenario ==0){scenario <- 9}
  ggplot() +
    geom_histogram(aes(x = est[1,,i,1], fill="Low"),binwidth = 0.001)+
    geom_histogram(aes(x = est[1,,i,2], fill="Medium"), binwidth = .001 )+
    geom_histogram(aes(x = est[1,,i,3], fill="High"), binwidth = .001)+
    geom_vline(xintercept = true_value_red[1,scenario,1]) +
    geom_vline(xintercept = true_value_red[1,scenario,2]) +
    geom_vline(xintercept = true_value_red[1,scenario,3]) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical MRD dist.") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                      "High" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist1_low[1:9], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))
mrd_dist2_low <- lapply(1:9, function(i){
  scenario <- i%%9
  if (scenario ==0){scenario <- 9}
  ggplot() +
    geom_histogram(aes(x = est[2,,i,1], fill="Low"),binwidth = 0.001)+
    geom_histogram(aes(x = est[2,,i,2], fill="Medium"), binwidth = .001 )+
    geom_histogram(aes(x = est[2,,i,3], fill="High"), binwidth = .001)+
    geom_vline(xintercept = true_value_red[2,scenario,1]) +
    geom_vline(xintercept = true_value_red[2,scenario,2]) +
    geom_vline(xintercept = true_value_red[2,scenario,3]) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical MRD dist.") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                      "High" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist2_low[1:9], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

mrd_dist3_low <- lapply(1:9, function(i){
  scenario <- i%%9
  if (scenario ==0){scenario <- 9}
  ggplot() +
    geom_histogram(aes(x = est[3,,i,1], fill="Low"),binwidth = 0.001)+
    geom_histogram(aes(x = est[3,,i,2], fill="Medium"), binwidth = .001 )+
    geom_histogram(aes(x = est[3,,i,3], fill="High"), binwidth = .001)+
    geom_vline(xintercept = true_value_red[3,scenario,1]) +
    geom_vline(xintercept = true_value_red[3,scenario,2]) +
    geom_vline(xintercept = true_value_red[3,scenario,3]) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical MRD dist.") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                      "High" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist3_low[1:9], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

mrd_dist4_low <- lapply(1:9, function(i){
  scenario <- i%%9
  if (scenario ==0){scenario <- 9}
  ggplot() +
    geom_histogram(aes(x = est[4,,i,1], fill="Low"),binwidth = 0.001)+
    geom_histogram(aes(x = est[4,,i,2], fill="Medium"), binwidth = .001 )+
    geom_histogram(aes(x = est[4,,i,3], fill="High"), binwidth = .001)+
    geom_vline(xintercept = true_value_red[4,scenario,1]) +
    geom_vline(xintercept = true_value_red[4,scenario,2]) +
    geom_vline(xintercept = true_value_red[4,scenario,3]) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical MRD dist.") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                      "High" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist4_low[1:9], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

mrd_dist5_low <- lapply(1:9, function(i){
  scenario <- i%%9
  if (scenario ==0){scenario <- 9}
  ggplot() +
    geom_histogram(aes(x = est[5,,i,1], fill="Low"),binwidth = 0.001)+
    geom_histogram(aes(x = est[5,,i,2], fill="Medium"), binwidth = .001 )+
    geom_histogram(aes(x = est[5,,i,3], fill="High"), binwidth = .001)+
    geom_vline(xintercept = true_value_red[5,scenario,1]) +
    geom_vline(xintercept = true_value_red[5,scenario,2]) +
    geom_vline(xintercept = true_value_red[5,scenario,3]) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical MRD dist.") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_fill_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                      "High" = "green"))
})
annotate_figure(ggarrange(plotlist = mrd_dist5_low[1:9], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

################BIAS, SD, MSE PLOTS ###################
bias_plots <- lapply(1:27, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_point[,i,1], colour = 'Low')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,1], colour = 'Low')) +
    geom_line(aes(x = 0:4, y = bias_point[,i,2], colour = 'Medium')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,2], colour = 'Medium')) +
    geom_line(aes(x = 0:4, y = bias_point[,i,3], colour = 'High')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,3], colour = 'High')) +
    #xlab(bquote(paste('N = ',.(scenarios[i,1]),", ", alpha[c] == .(scenarios[i,2]),', ',alpha[a] == .(scenarios[i,3])))) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    
    ylab("Empirical bias") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_color_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                       "High" = "green")) +
    ylim(-0.1,0.1)
})
annotate_figure(ggarrange(plotlist = bias_plots[1:27], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

sd_plots <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = sd_point[,i,1], colour = 'Low')) +
    geom_point(aes(x = 0:4, y = sd_point[,i,1], colour = 'Low')) +
    geom_line(aes(x = 0:4, y = sd_point[,i,2], colour = 'Medium')) +
    geom_point(aes(x = 0:4, y = sd_point[,i,2], colour = 'Medium')) +
    geom_line(aes(x = 0:4, y = sd_point[,i,3], colour = 'High')) +
    geom_point(aes(x = 0:4, y = sd_point[,i,3], colour = 'High')) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical SD") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), ) +  
    scale_color_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                       "High" = "green")) +
    ylim(0,0.28)
})
annotate_figure(ggarrange(plotlist = sd_plots[1:9], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

bias.se_plots <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_point[,i,1]/sd_point[,i,1], colour = 'Low')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,1]/sd_point[,i,1], colour = 'Low')) +
    geom_line(aes(x = 0:4, y = bias_point[,i,2]/sd_point[,i,2], colour = 'Medium')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,2]/sd_point[,i,2], colour = 'Medium')) +
    geom_line(aes(x = 0:4, y = bias_point[,i,3]/sd_point[,i,3], colour = 'High')) +
    geom_point(aes(x = 0:4, y = bias_point[,i,3]/sd_point[,i,3], colour = 'High')) +
    #xlab(bquote(paste('N = ',.(scenarios[i,1]),", ", alpha[c] == .(scenarios[i,2]),', ',alpha[a] == .(scenarios[i,3])))) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    
    ylab("Empirical bias") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14)) +  
    scale_color_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                       "High" = "green")) +
    ylim(-0.1,0.1)
})
annotate_figure(ggarrange(plotlist = bias.se_plots[1:9], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))

mse_plots <- lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = sqrt(bias_point[,i,1]^2 + sd_point[,i,1]^2), colour = 'Low')) +
    geom_point(aes(x = 0:4, y = sqrt(bias_point[,i,1]^2 + sd_point[,i,1]^2), colour = 'Low')) +
    geom_line(aes(x = 0:4, y = sqrt(bias_point[,i,2]^2 + sd_point[,i,2]^2), colour = 'Medium')) +
    geom_point(aes(x = 0:4, y = sqrt(bias_point[,i,2]^2 + sd_point[,i,2]^2), colour = 'Medium')) +
    geom_line(aes(x = 0:4, y = sqrt(bias_point[,i,3]^2 + sd_point[,i,3]^2), colour = 'High')) +
    geom_point(aes(x = 0:4, y = sqrt(bias_point[,i,3]^2 + sd_point[,i,3]^2), colour = 'High')) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Empirical root-MSE") + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), ) +  
    scale_color_manual(name = "Event rate", values = c("Low"= "red", "Medium" = "blue",
                                                       "High" = "green")) +
    ylim(0,0.3)
})
annotate_figure(ggarrange(plotlist = mse_plots[1:9], nrow = 3, ncol = 9, common.legend = T, legend = 'bottom'))


# ############# COVERAGE #####################
# 
# coverage_ind <- array(0,dim = c(4,5,9,3))
# success <- array(0,dim = c(4,5,9,3))
# 
# for (i in 1:1000){
#   for (k in 1:5){
#     for (j in 1:9){
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

pivot_coverage_ind <- array(0,dim = c(6,5,27,3))
pivot_success <- array(0,dim = c(6,5,27,3))

for (i in 1:iters){
  for (k in 1:5){
    for (j in 1:27){
      for (l in 1){
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

bias_elim_pivot_coverage_ind <- array(0,dim = c(4,5,27,3))
bias_elim_pivot_success <- array(0,dim = c(4,5,27,3))

for (i in 1:iters){
  for (k in 1:5){
    for (j in 1:27){
      for (l in 1:3){
        scenario <- j%%9
        if (scenario ==0){scenario <- 9}
        
        if (is.na(bootstrap[k,1,i,j,l]) == F){
          bias_elim_pivot_success[1,k,j,l] <- bias_elim_pivot_success[1,k,j,l] + 1
          if (all(2*est[k,i,j,l] + bootstrap[k,1,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
              & all(2*est[k,i,j,l] + bootstrap[k,2,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
            bias_elim_pivot_coverage_ind[1,k,j,l] <- bias_elim_pivot_coverage_ind[1,k,j,l] + 1
          }
        }
        
        if (is.na(LEF_outcome[k,1,i,j,l]) == F){
          bias_elim_pivot_success[2,k,j,l] <- bias_elim_pivot_success[2,k,j,l] + 1
          if (all(2*est[k,i,j,l] + LEF_outcome[k,1,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
              & all(2*est[k,i,j,l] + LEF_outcome[k,2,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
            bias_elim_pivot_coverage_ind[2,k,j,l] <- bias_elim_pivot_coverage_ind[2,k,j,l] + 1
          }
        }
        
        if (is.na(LEF_both[k,1,i,j,l]) == F){
          bias_elim_pivot_success[3,k,j,l] <- bias_elim_pivot_success[3,k,j,l] + 1
          if (all(2*est[k,i,j,l] + LEF_both[k,1,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
              & all(2*est[k,i,j,l] + LEF_both[k,2,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
            bias_elim_pivot_coverage_ind[3,k,j,l] <- bias_elim_pivot_coverage_ind[3,k,j,l] + 1
          }
        }
        
        if (all(is.na(sandwich[k,1,i,j,l])) == F){
          bias_elim_pivot_success[4,k,j,l] <- bias_elim_pivot_success[4,k,j,l] + 1
          if (all(-sandwich[k,2,i,j,l] <= bias_point[k,j,l] +true_value_red[k,scenario, l])
              & all(-sandwich[k,1,i,j,l] >= bias_point[k,j,l] +true_value_red[k,scenario, l])){
            bias_elim_pivot_coverage_ind[4,k,j,l]<- bias_elim_pivot_coverage_ind[4,k,j,l]+ 1
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
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'pink' )) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Coverage") +  ylim(0.2,1) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14))
}
)
annotate_figure(ggarrange(plotlist = coverage_low, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

coverage_med <-  lapply(1:9, function(i){
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
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'pink' )) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Coverage") +  ylim(0.2,1) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14))
}
)
annotate_figure(ggarrange(plotlist = coverage_med[1:9], nrow = 3, ncol = 3, common.legend = T,
                          legend = 'bottom'))

coverage_high <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = pivot_coverage_ind[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = pivot_coverage_ind[4,,i,3], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Coverage") +  ylim(0.2,1) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14))
}
)
annotate_figure(ggarrange(plotlist = coverage_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

######################   BIAS ELIM PIVOT COVERAGE ########################
coverage_low <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,1], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,1], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,1], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,1], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,1], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,1], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,1], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,1], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Bias-eliminated coverage") +  ylim(0.2,1) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))
}
)
annotate_figure(ggarrange(plotlist = coverage_low, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

coverage_med <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,2], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,2], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,2], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,2], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,2], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,2], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,2], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Bias-eliminated coverage") +  ylim(0.2,1) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = coverage_med, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

coverage_high <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = bias_elim_pivot_coverage_ind[4,,i,3], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Bias-eliminated coverage") +  ylim(0.2,1) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = coverage_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

###################### CI WIDTH ##############################

CI_width <- array(, dim = c(6,5,9,3))
for (k in 1:5){
  for (j in 1:9){
    for (l in 1:3){
      CI_width[1,k,j,l] <- mean(abs(bootstrap[k,2,,j,l] - bootstrap[k,1,,j,l]), na.rm = TRUE)
      CI_width[2,k,j,l] <- mean(abs(LEF_outcome[k,2,,j,l] - LEF_outcome[k,1,,j,l]), na.rm = TRUE)
      CI_width[3,k,j,l] <- mean(abs(LEF_both[k,2,,j,l] - LEF_both[k,1,,j,l]), na.rm = TRUE)
      CI_width[4,k,j,l] <- mean(abs(sandwich[k,2,,j,l] - sandwich[k,1,,j,l]), na.rm = TRUE)
      CI_width[5,k,j,l] <- mean(abs(jackknife_mvn[k,2,,j,l] - jackknife_mvn[k,1,,j,l]), na.rm = TRUE)
      CI_width[6,k,j,l] <- mean(abs(jackknife_wald[k,2,,j,l] - jackknife_wald[k,1,,j,l]), na.rm = TRUE)
    }}}

width_low <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width[1,,i,1], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width[1,,i,1], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width[2,,i,1], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width[2,,i,1], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width[3,,i,1], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width[3,,i,1], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width[4,,i,1], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width[4,,i,1], colour = "Sandwich")) + 
    geom_line(aes(x = 0:4, y = CI_width[5,,i,1], colour = "Jackknife MVN")) +
    geom_point(aes(x = 0:4, y = CI_width[5,,i,1], colour = "Jackknife MVN")) +
    geom_line(aes(x = 0:4, y = CI_width[6,,i,1], colour = "Jackknife Wald")) +
    geom_point(aes(x = 0:4, y = CI_width[6,,i,1], colour = "Jackknife Wald")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple",
                                                    "Jackknife MVN" = 'orange',"Jackknife Wald" = 'pink' )) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Average Width") +  ylim(0,1.5) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))
}
)
annotate_figure(ggarrange(plotlist = width_low, nrow = 3, ncol = 3, common.legend = T,
                          legend = 'bottom'))

width_med <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width[1,,i,2], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width[1,,i,2], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width[2,,i,2], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width[2,,i,2], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width[3,,i,2], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width[3,,i,2], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width[4,,i,2], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Average Width") +  ylim(0,0.7) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_med, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

width_high <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = CI_width[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = CI_width[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = CI_width[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = CI_width[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = CI_width[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = CI_width[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = CI_width[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = CI_width[4,,i,3], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Average Width") +  ylim(0,0.7) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

############################# CI WIDTH SD###############################
CI_width_sd <- array(, dim = c(4,5,9,3))
for (k in 1:5){
  for (j in 1:9){
    for (l in 1:3){
      CI_width_sd[1,k,j,l] <- sd(abs(bootstrap[k,2,,j,l] - bootstrap[k,1,,j,l]), na.rm = TRUE)
      CI_width_sd[2,k,j,l] <- sd(abs(LEF_outcome[k,2,,j,l] - LEF_outcome[k,1,,j,l]), na.rm = TRUE)
      CI_width_sd[3,k,j,l] <- sd(abs(LEF_both[k,2,,j,l] - LEF_both[k,1,,j,l]), na.rm = TRUE)
      CI_width_sd[4,k,j,l] <- sd(abs(sandwich[k,2,,j,l] - sandwich[k,1,,j,l]), na.rm = TRUE)
    }}}

width_low <-  lapply(1:9, function(i){
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
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width SD") +  ylim(0,0.7) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))
}
)
annotate_figure(ggarrange(plotlist = width_low, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))


width_med <-  lapply(1:9, function(i){
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
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width SD") +  ylim(0,0.7) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_med, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

width_high <-  lapply(1:9, function(i){
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
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width SD") +  ylim(0,0.7) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

########################################### CI WIDTH MEDIAN ##########################
CI_width_median <- array(, dim = c(4,5,9,3))
for (k in 1:5){
  for (j in 1:9){
    for (l in 1:3){
      CI_width_median[1,k,j,l] <- median(abs(bootstrap[k,2,,j,l] - bootstrap[k,1,,j,l]), na.rm = TRUE)
      CI_width_median[2,k,j,l] <- median(abs(LEF_outcome[k,2,,j,l] - LEF_outcome[k,1,,j,l]), na.rm = TRUE)
      CI_width_median[3,k,j,l] <- median(abs(LEF_both[k,2,,j,l] - LEF_both[k,1,,j,l]), na.rm = TRUE)
      CI_width_median[4,k,j,l] <- median(abs(sandwich[k,2,,j,l] - sandwich[k,1,,j,l]), na.rm = TRUE)
    }}}

width_low <-  lapply(1:9, function(i){
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
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width median") +  ylim(0,0.8) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))
}
)
annotate_figure(ggarrange(plotlist = width_low, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))


width_med <-  lapply(1:9, function(i){
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
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width median") +  ylim(0,0.8) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_med, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

width_high <-  lapply(1:9, function(i){
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
    xlab(bquote(atop(N ==.(scenarios[i,1]), alpha[c] ==.(scenarios[i,2])~', '~alpha[a] == .(scenarios[i,3])))) +
    ylab("Width median") +  ylim(0,0.8) + 
    theme(title=element_text(size=15), axis.text = element_text(size=10), legend.text = element_text(size=14), axis.title.y = element_text(size = 12))}
)
annotate_figure(ggarrange(plotlist = width_high, nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'))

######################PIVOT FAILURE RATE #######################
failure_table <- data.frame(matrix(,nrow = 0, ncol = 8))
for (i in 1:9){
  for (j in 1:3){
    failure_table <- rbind(failure_table, c(outcomes[j], scenarios[i,1], scenarios[i,2],scenarios[i,3], 
                                            (iters - pivot_success[,1,i,j])/iters))
  }
}
colnames(failure_table) <- c('Outcome_prevalence', 'Sample_size', 'Confounding', 'Treatment_prevalence', 'Bootstrap', 'LEF_outcome',
                             'LEF_both', 'Sandwich')
failure_table$Sample_size <- as.numeric(failure_table$Sample_size)
failure_table <- failure_table %>% 
  dplyr::group_by(Outcome_prevalence, Sample_size, Confounding, Treatment_prevalence) %>% 
  dplyr::summarise(Bootstrap = mean(as.numeric(Bootstrap)),
                   LEF_outcome = mean(as.numeric(LEF_outcome)),
                   LEF_both = mean(as.numeric(LEF_both)),
                   Sandwich = mean(as.numeric(Sandwich)))

################ MC SE PLOTS #######################
MC_SE<- array(, dim = c(4,5,9,3))
for (i in 1:9){
  for (j in 1:3){
    MC_SE[,,i,j] <- sqrt((coverage_ind[,,i,j]*(1-coverage_ind[,,i,j]))/iters)
  }
}

MC_SE_pivot<- array(, dim = c(4,5,9,3))
for (i in 1:9){
  for (j in 1:3){
    MC_SE_pivot[,,i,j] <- sqrt((pivot_coverage_ind[,,i,j]*(1-pivot_coverage_ind[,,i,j]))/iters)
  }
}

MCSE_low <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = MC_SE[1,,i,1], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = MC_SE[1,,i,1], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = MC_SE[2,,i,1], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = MC_SE[2,,i,1], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = MC_SE[3,,i,1], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = MC_SE[3,,i,1], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = MC_SE[4,,i,1], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = MC_SE[4,,i,1], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    labs(x = 'Follow up time',
         y = "MC SE",
         title = paste("N =", scenarios[i,1],
                       '\nConfounding =',scenarios[i,2],
                       '\nTreat. prev. =', scenarios[i,3]))+ ylim(0,0.03) + 
    theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
}
)
annotate_figure(ggarrange(plotlist = MCSE_low[1:9], nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'), top = 'Low event rate')

MCSE_med <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = MC_SE[1,,i,2], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = MC_SE[1,,i,2], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = MC_SE[2,,i,2], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = MC_SE[2,,i,2], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = MC_SE[3,,i,2], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = MC_SE[3,,i,2], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = MC_SE[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = MC_SE[4,,i,2], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    labs(x = 'Follow up time',
         y = "MC SE",
         title = paste("N =", scenarios[i,1],
                       '\nConfounding =',scenarios[i,2],
                       '\nTreat. prev. =', scenarios[i,3]))+ ylim(0,0.03) + 
    theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
}
)
annotate_figure(ggarrange(plotlist = MCSE_med[1:9], nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'), top = 'Medium event rate')

MCSE_high <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = MC_SE[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = MC_SE[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = MC_SE[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = MC_SE[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = MC_SE[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = MC_SE[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = MC_SE[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = MC_SE[4,,i,3], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    labs(x = 'Follow up time',
         y = "MC SE",
         title = paste("N =", scenarios[i,1],
                       '\nConfounding =',scenarios[i,2],
                       '\nTreat. prev. =', scenarios[i,3]))+ ylim(0,0.03) + 
    theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
}
)
annotate_figure(ggarrange(plotlist = MCSE_high[1:9], nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'), top = 'High event rate')

MCSE_low <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[1,,i,1], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[1,,i,1], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[2,,i,1], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[2,,i,1], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[3,,i,1], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[3,,i,1], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[4,,i,1], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[4,,i,1], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    labs(x = 'Follow up time',
         y = "MC SE",
         title = paste("N =", scenarios[i,1],
                       '\nConfounding =',scenarios[i,2],
                       '\nTreat. prev. =', scenarios[i,3]))+ ylim(0,0.03) + 
    theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
}
)
annotate_figure(ggarrange(plotlist = MCSE_low[1:9], nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'), top = 'Low event rate')

MCSE_med <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[1,,i,2], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[1,,i,2], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[2,,i,2], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[2,,i,2], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[3,,i,2], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[3,,i,2], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[4,,i,2], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[4,,i,2], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    labs(x = 'Follow up time',
         y = "MC SE",
         title = paste("N =", scenarios[i,1],
                       '\nConfounding =',scenarios[i,2],
                       '\nTreat. prev. =', scenarios[i,3])) + 
    theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
}
)
annotate_figure(ggarrange(plotlist = MCSE_med[1:9], nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'), top = 'Medium event rate')

MCSE_high <-  lapply(1:9, function(i){
  ggplot() +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[1,,i,3], colour = "Bootstrap")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[1,,i,3], colour = "Bootstrap")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[2,,i,3], colour = "LEF outcome")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[2,,i,3], colour = "LEF outcome")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[3,,i,3], colour = "LEF both")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[3,,i,3], colour = "LEF both")) +
    geom_line(aes(x = 0:4, y = MC_SE_pivot[4,,i,3], colour = "Sandwich")) +
    geom_point(aes(x = 0:4, y = MC_SE_pivot[4,,i,3], colour = "Sandwich")) +
    scale_color_manual(name = "CI type", values = c("Bootstrap"= "red", "Sandwich" = "blue",
                                                    "LEF outcome" = "green", "LEF both" = "purple")) +
    labs(x = 'Follow up time',
         y = "MC SE",
         title = paste("N =", scenarios[i,1],
                       '\nConfounding =',scenarios[i,2],
                       '\nTreat. prev. =', scenarios[i,3]))+ ylim(0,0.03) + 
    theme(plot.title = element_text(size=10))+ theme(aspect.ratio = 1, axis.title = element_text(size = 10))
}
)
annotate_figure(ggarrange(plotlist = MCSE_high[1:9], nrow = 3, ncol = 9, common.legend = T,
                          legend = 'bottom'), top = 'High event rate')

################WEIGHTS########################
weights <- data.frame(matrix(,nrow = 0, ncol = 7))
for(i in 1:9){
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

for(i in 1:9){
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

for(i in 1:9){
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
