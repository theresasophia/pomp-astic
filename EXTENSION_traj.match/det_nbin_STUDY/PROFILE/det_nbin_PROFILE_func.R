rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows
# close graphics windows
library(pomp)
library(magrittr)
library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
library(foreach)
library(doParallel)
library("MASS")
library("quantreg")
registerDoParallel()
stopifnot(packageVersion("pomp") >= "1.7")

## FUNKTIONS
appl<- function(pd,f){apply(x(pd),2,FUN=f)}

h<- function(x) {summary(mf[[x]])$convergence}
g <- function(x) {coef(mf[[x]])}
f <- function(x) {summary(mf[[x]])$loglik}
x <- function(x){matrix(c(seq(from=1,to=dim(x)[1],by=1)),nrow=1,ncol=dim(x)[1])}

poly <- function(x,coef_quant){coef_quant[1]+coef_quant[2]*x+coef_quant[3]*x^2}

chi_half <- 1.92

res<- function(liks,coef){data.frame(logLik=liks,beta1=coef["beta1",],beta2=coef["beta2",],beta3=coef["beta3",],
           beta11=coef["beta11",],phi=coef["phi",],od=coef["od",])}
qant_reg <- function(x){ quantreg::rq( logLik ~ 1 +x+ I(x^2), tau=0.9,data=results_global)}

max_x <- function(coef_quant){-coef_quant[2]/(2*coef_quant[3])}


lower <- function(max_fun,chi_half,coef_quant){-sqrt((max_fun-chi_half-coef_quant[1])/coef_quant[3]+
                                             coef_quant[2]^2/(coef_quant[3]^2*4))-coef_quant[2]/coef_quant[3]*1/2}
upper <- function(max_fun,chi_half,coef_quant){sqrt((max_fun-chi_half-coef_quant[1])/coef_quant[3]+
                                                       coef_quant[2]^2/(coef_quant[3]^2*4))-coef_quant[2]/coef_quant[3]*1/2}




save(h,f,g,x,poly,chi_half,qant_reg,appl,res,lower,max_x, upper,file = "det_nbin_PROFILE_study_func.rda")
