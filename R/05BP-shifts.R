# Migration application for multiple shifts

rm(list=ls())
# Packages ----
setwd("C:/Users/na1d22/Dropbox/Surrey/IVSelection") # Adjust.
library(haven)
library(AER)
library(ivpack)
library(gtools)
library(glmnet)
library(dplyr)

# Load data ----
analysis_data <- read_dta("Stata/data/analysis_shifts_data.dta")
analysis_data <- na.omit(analysis_data)

attach(analysis_data)

# Functions ---- 
source("R/Code/01Functions_060723Newtry.R")

# Sample restriction
analysis_data <- filter(analysis_data, survey >= 1990)

# Variables ----
# Instruments
dataprep(wFB="^shs")

# Analyses ---- P = 1
adamAR(Y1, D1)
alpha.ad
l
which(a.end == 1)
adamAR(Y2, D1)
l
which(a.end == 1)
adamAR(Y3, D1)
l
which(a.end == 1)

# Analyses ---- P = 2
adamAR(Y1, D)
alpha.ad
l
which(a.end == 1)
adamAR(Y2, D)
l
which(a.end == 1)
adamAR(Y3, D)
l
which(a.end == 1)

# Ada Multiple, HS

adamHS(Y1, D, tau = 0.1)
alpha.ad
l
which(a.end == 1)
adamHS(Y2, D, tau = 0.1)
l
which(a.end == 1)
adamHS(Y3, D, tau = 0.1)
l
which(a.end == 1)

detach(analysis_data)
