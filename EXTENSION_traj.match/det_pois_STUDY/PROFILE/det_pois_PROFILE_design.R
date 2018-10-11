rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows
# close graphics windows
library(pomp)
library(magrittr)
library(plyr)
library(reshape2)
library("rootSolve")
library(ggplot2)
library(scales)
library(foreach)
library(doParallel)
library("MASS")
library("quantreg")
registerDoParallel()
stopifnot(packageVersion("pomp") >= "1.7")

load("~/Dropbox/AAPAPER/R_code/traj.match/det_pois_model_study.rda")
##PROFILE DESIGNS
profileDesign(
  beta1=seq(from=12.625,to=12.675,length=10),
  lower=c(beta2=0.2, beta3=0.4, beta11=0.15,phi=0.02,sir_fixed_params ),
  upper=c(beta2=0.3, beta3=0.5, beta11=0.16,phi=0.04,sir_fixed_params ),
  nprof=10
) -> pd_beta1


pairs(~beta1+beta2+beta3+phi+delta1,data=pd_beta1)

profileDesign(
  beta2=seq(from=0.234,to=0.24,length=10),
  lower=c(beta1=12, beta3=0.4, beta11=0.15,phi=0.02,sir_fixed_params ),
  upper=c(beta1=13, beta3=0.5, beta11=0.16,phi=0.03,sir_fixed_params ),
  nprof=10
) -> pd_beta2


profileDesign(
  beta3=seq(from=0.419,to=0.422,length=10),
  lower=c(beta1=12, beta2=0.2, beta11=0.15,phi=0.03,sir_fixed_params ),
  upper=c(beta1=13, beta2=0.3, beta11=0.16,phi=0.04,sir_fixed_params ),
  nprof=10
) -> pd_beta3
dim(pd_beta3)

profileDesign(
  beta11=seq(from=0.1496,to=0.1504,length=20),
  lower=c(beta1=12, beta2=0.2, beta3=0.4,phi=0.03,sir_fixed_params ),
  upper=c(beta1=13, beta2=0.3, beta3=0.5,phi=0.04,sir_fixed_params ),
  nprof=5
) -> pd_beta11

profileDesign(
  phi=seq(from=0.09,to=0.11,length=10),
  #lower=coef(mifs_global[[best]])[-1],upper=coef(mifs_global[[best]])[-1],
  lower=c(beta1=12, beta2=0.2, beta11=0.15,beta3=0.4,sir_fixed_params ),
  upper=c(beta1=13, beta2=0.3, beta11=0.16,beta3=0.5,sir_fixed_params ),
  nprof=5
) -> pd_phi


save(pd_beta1, pd_beta2,pd_beta3 , pd_beta11 ,pd_phi,file = "det_pois_PROFILE_study_design.rda")

