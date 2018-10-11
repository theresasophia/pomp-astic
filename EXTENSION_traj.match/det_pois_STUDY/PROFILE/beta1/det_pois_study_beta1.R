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


load("~/Dropbox/AAPAPER/R_code/traj.match/det_pois_model_study.rda")
load("~/Dropbox/AAPAPER/R_code/traj.match/det_pois_STUDY/PROFILE/det_pois_PROFILE_study_design.rda")
load("~/Dropbox/AAPAPER/R_code/traj.match/det_pois_STUDY/PROFILE/det_pois_PROFILE_study_func.rda")
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)


###################BETA1###################################################
pd <- pd_beta1
dim(pd)
stew(file="prof_gam_beta1_study_first-%d.rda",{
  
  t_global <- system.time({
    mf <- foreach(p=iter(pd,"row"),
                  .packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
                    traj.match(study,
                               start=unlist(p),
                               est=c("beta2","beta3","beta11","phi"),
                               method="Nelder-Mead",reltol=1e-6,maxit=2000, transform=TRUE)
                    
                  }
  })
},seed=1270401374,kind="L'Ecuyer")

appl<- function(pd,f){apply(x(pd),2,FUN=f)}
liks <- appl(pd,f)
coef <- appl(pd,g)
conv <- appl(pd,h)
stopifnot(conv==0)
length(which(conv==0))
trans <-t(coef)

results_global <- res(liks, coef)
pairs(~logLik+beta1,data=subset(results_global,logLik>max(logLik)-10))

stew(file="prof_gam_beta1_study_second-%d.rda",{
  
  t_global <- system.time({
    guesses <- as.data.frame(trans)
    mf <- foreach(guess=iter(guesses,"row"),.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      traj.match(study,
                 start=unlist(guess),
                 est=c("beta2","beta3","beta11","phi"),
                 method="Nelder-Mead",reltol=1e-5,maxit=2000, transform=TRUE)
    }
  })
},seed=1270401374,kind="L'Ecuyer")


appl<- function(pd,f){apply(x(pd),2,FUN=f)}
liks <- appl(pd,f)
coef <- appl(pd,g)
conv <- appl(pd,h)
stopifnot(conv==0)
length(which(conv==0))

results_global <- res(liks, coef)
pairs(~logLik+beta1,data=results_global)

results_global%>% 
  mutate(beta1=signif(beta1,8)) %>%
  ddply(~beta1,subset,logLik==max(logLik)) %>%
  ggplot(aes(x=beta1,y=logLik))+geom_point() +
  geom_quantile(formula=y~1+x+I(x^2), quantiles=0.9,col="blue")

# linear interpolation
interp <- with(results_global,approxfun(beta1,logLik))
uni_fun<- function(x){ interp(x) - (max(results_global$logLik)-1.92)}
x<- seq(10,15,0.001)
plot(x,uni_fun(x),type="l")
library("rootSolve")
int_beta1 <- sprintf("[%.3f, %.3f]", uniroot.all(uni_fun, lower=12, upper=13)[1],uniroot.all(uni_fun, lower=12, upper=13)[2])
save(int_beta1,file = "~/Dropbox/AAPAPER/R_code/traj.match/det_pois_STUDY/beta1_pois.rda")
