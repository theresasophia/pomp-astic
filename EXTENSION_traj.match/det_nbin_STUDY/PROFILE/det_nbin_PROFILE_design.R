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

load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_model_study.rda")
##PROFILE DESIGNS
profileDesign(
  beta1=seq(from=12,to=13,length=20),
  lower=c(beta2=0.2, beta3=0.4, beta11=0.15,phi=0.02,od=0.4,sir_fixed_params ),
  upper=c(beta2=0.3, beta3=0.5, beta11=0.16,phi=0.04,od=0.6,sir_fixed_params ),
  nprof=5
) -> pd_beta1


pairs(~beta1+beta2+beta3+phi+delta1,data=pd_beta1)

profileDesign(
  beta2=seq(from=0.2,to=0.3,length=20),
  lower=c(beta1=12, beta3=0.4, beta11=0.15,phi=0.02,od=0.4,sir_fixed_params ),
  upper=c(beta1=13, beta3=0.5, beta11=0.16,phi=0.03,od=0.6,sir_fixed_params ),
  nprof=5
) -> pd_beta2


profileDesign(
  beta3=seq(from=0.35,to=0.5,length=20),
  lower=c(beta1=12, beta2=0.2, beta11=0.15,phi=0.03,od=0.4,sir_fixed_params ),
  upper=c(beta1=13, beta2=0.3, beta11=0.16,phi=0.04,od=0.6,sir_fixed_params ),
  nprof=5
) -> pd_beta3
dim(pd_beta3)

profileDesign(
  beta11=seq(from=0.14,to=0.17,length=30),
  lower=c(beta1=12, beta2=0.2, beta3=0.4,phi=0.03,od=0.4,sir_fixed_params ),
  upper=c(beta1=13, beta2=0.3, beta3=0.5,phi=0.04,od=0.6,sir_fixed_params ),
  nprof=5
) -> pd_beta11

profileDesign(
  phi=seq(from=0.05,to=0.2,length=30),
  #lower=coef(mifs_global[[best]])[-1],upper=coef(mifs_global[[best]])[-1],
  lower=c(beta1=12, beta2=0.2, beta11=0.15,beta3=0.4,od=0.4,sir_fixed_params ),
  upper=c(beta1=13, beta2=0.3, beta11=0.16,beta3=0.5,od=0.6,sir_fixed_params ),
  nprof=5
) -> pd_phi

profileDesign(
  od=seq(from=0.4,to=0.6,length=20),
  #lower=coef(mifs_global[[best]])[-1],upper=coef(mifs_global[[best]])[-1],
  lower=c(beta1=12, beta2=0.2, beta11=0.15,beta3=0.4,phi=0.04,sir_fixed_params ),
  upper=c(beta1=13, beta2=0.3, beta11=0.16,beta3=0.5,phi=0.04,sir_fixed_params ),
  nprof=5
) -> pd_od


save(pd_beta1, pd_beta2,pd_beta3 , pd_beta11 ,pd_phi, pd_od, file = "det_nbin_PROFILE_study_design.rda")

