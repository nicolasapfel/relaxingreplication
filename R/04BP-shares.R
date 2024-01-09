# This produces the results for Anderson-Rubin selection
# Single and multiple

rm(list=ls())
# Packages ----
setwd("C:/Users/na1d22/Dropbox/Surrey/IVSelection/R/Code") # Adjust

source("01Functions_repo.R")

# Load data ----
analysis_data <- read_dta("../../Stata/data/analysis_data.dta")
analysis_data <- na.omit(analysis_data)

attach(analysis_data)
# Sample restriction
analysis_data <- filter(analysis_data, survey >= 1990)

dataprep(wFB="^share70_")
# weighted versions of variables in wY1, wZ etc.
# weighted and FWL-ed version of variables in Y1, Z, etc. 

# Compare FWLed and weighted with weighted only
# IV regressions
ss <- ivreg(Y1 ~ D1 - 1 | Z - 1 , subset = survey>=1990)
summary(ss)
ss <- ivreg(wY1 ~ wD1 + wyear3 + wyear4 + sqrt(lpop) - 1 | wZ + wyear3 + wyear4 + sqrt(lpop) - 1, subset = survey>=1990)
summary(ss)
# Compared with STATA - use FWLed and weighted - gives same results as STATA

# Anderson-Rubin tests 
# I need and.rub(), because ivmodel() can only handle single regressor
and.rub(Y=wY1, D=wD1, Z=wZ[,c(6,11)], X=cbind(wX, Z[,1]))
AR
che <- ivmodel(Y=as.numeric(Y1), D=D1, Z=Z[,c(6,11)], X=Z[,1], intercept = F, 
               beta0 = 0, alpha = 0.05, k = c(0, 1), 
               heteroSE = T, clusterID = NULL,
               deltarange = NULL, na.action = na.omit)
bl <- che$LIML[1]$point.est[1,1]
(kappa <- che$LIML$k)
2166*log(kappa)
# ivmodel and my and.rub test give same result with a single regressor, with all IVs
# when using weighted data but no FWL with and.rub()
# and weights and FWL with ivmodel()
# Same results also when using fewer IVs
# At this point, I'm confident that and.rub is doing the correct thing
# and will do the correct thing also when using multiple endog. regressors

# Analyses single ----
#rm(list=setdiff(ls(),c("analysis_data", "fnMatSqrtInverse", "and.rub", "dataprep", "adam" )))

#dlweekly
adamAR(Y=Y1, D=D1)
l
ivmod <- ivreg(Y1 ~ D1 + Z[,which(a.end==1)] -1 | Z -1 )
summary(ivmod) # SEs not valid because of FWL
which(a.end==1)
# Go to Stata (01Analysis-Shares.do) and estimate, using these as invalid

#dlweekly_hskill
adamAR(Y=Y2, D=D1)
l
sum(a.end)
ivmod <- ivreg(Y2 ~ D1 + Z[,which(a.end==1)] -1 | Z -1 )
summary(ivmod)
which(a.end==1)

#dlweekly_lskill
adamAR(Y=Y3,D1)
l
ivmod <- ivreg(Y3 ~ D1 + Z[,which(a.end==1)] -1 | Z -1 )
summary(ivmod)
which(a.end==1)

# Analyses multiple ----
analysis_data <- read_dta("../../Stata/data/analysis_data.dta")
analysis_data <- na.omit(analysis_data)
analysis_data <- filter(analysis_data, survey >= 1990)
dataprep(wFB="^share70|^share80")
colnames(Z)

# adam with HS test
# dlweekly
adamHS(Y=Y1, D=D, tau=0.1)
l
which(a.end == 1)
# None selected as invalid

# dlweekly_hskill
adamHS(Y=Y2, D, tau=0.1/log(n))
l
which(a.end == 1)
adamHS(Y=Y2, D, tau=0.1)
l
which(a.end == 1)
ivmod <- ivreg(Y2 ~ D + Z[,which(a.end==1)] -1 | Z )
cluster.robust.se(ivmod, czone)

# dlweekly_lskill
adamHS(Y=Y3, D=D, tau=0.1)
l
which(a.end == 1)
# None selected as invalid

# adam with AR testing procedure
adamAR(Y=Y1, D=D)
which(a.end==1)

adamAR(Y=Y2, D=D)
which(a.end==1)

adamAR(Y=Y3, D=D)
which(a.end==1)
# Go to STATA, 01Analysis-shares.do, to estimate

detach(analysis_data)
