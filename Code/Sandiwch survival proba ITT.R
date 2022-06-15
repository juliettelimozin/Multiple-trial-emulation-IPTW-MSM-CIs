## simulate data for testing TrialEmulation package, using the algorithm in Young and Tchetgen Tchetgen (2014) 
library(modelr)
library(reshape2)
library(tidyverse)
library(tidyr)
setwd("/Users/juliette/Documents/MPhil PHS 21-22/MPhil-dissertation/Code")
source("simulate_MSM.R")
set.seed(20222022)
library(MASS)
library(RandomisedTrialsEmulation, lib.loc='/home/li/lib/R/R_LIBS/')


#RUN BOOTSTRAP FILES TO WORK ON SAME SIMULATED DATA AND MODEL FIT
load("ITT_fit.rda")
predicted_probas_ITT <- read.csv("predicted_probas_ITT")
fitting_data_treatment <- read.csv("fitting_data_treatment_ITT.csv")
fitting_data_control <- read.csv("fitting ")
coeffs_sample <- mvrnorm(2000,ITT$model$model$coefficients, vcov(ITT$model$model))

for (i in 1:2000){
  
}

  
  
  
  