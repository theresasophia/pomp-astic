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
library(plyr)
library(tidyr)
registerDoParallel()
stopifnot(packageVersion("pomp") >= "1.7")

load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_model_data.rda")


require(doParallel)
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

N <- 100 # number of loops
coef_list <- vector("list", N)

for(i in 1:N){
  
  seed<-i
  set.seed(seed)

  
  sim <- trajectory(sir,as.data.frame=TRUE)
  data.frame(time=seq(1:416),cases1=rpois(416,sim$H1[-1]),cases2=rpois(416,sim$H2[-1]),cases3=rpois(416,sim$H3[-1]))%>%
    rbind(data.frame(time=0,cases1=NA,cases2=NA,cases3=NA)) %>%
    arrange(time) -> sim_study
  
  
  pomp(data = sim_study,
       times="time",
       t0=1-6*52,
       dmeasure = dmeas,
       rmeasure = rmeas,
       rprocess = euler.sim(step.fun = sir.step, delta.t = 1/10),
       statenames = c("S1", "I1", "R1", "H1", "S2", "I2", "R2", "H2","S3","I3", "R3", "H3"),
       paramnames = names(params),
       skeleton=vectorfield(Csnippet(sir.skel)),
       zeronames=c("H1", "H2", "H3"),
       initializer= init,
       toEstimationScale=toEst,
       fromEstimationScale=fromEst,
       params = params
  ) -> study

#box for starting values
sir_box <- rbind(
  beta1=c(12,16),
  beta2=c(0.2,0.4),
  beta3=c(0.3,0.5),
  beta11=c(0.12,0.16),
  phi=c(0,0.2)
)

# number of daws from box
n<- 100

#save starting values and solution

  
  t_global <- system.time({
    guesses <- as.data.frame(apply(sir_box,1,function(x)runif(n,x[1],x[2])))
    det_pois_study_global <- foreach(guess=iter(guesses,"row"),.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      traj.match(study,
                 start=c(unlist(guess),sir_fixed_params),
                 est=c("beta1","beta2","beta3","beta11","phi"),
                 method="Nelder-Mead",reltol=1e-6,maxit=2000, transform=TRUE)
      
    }
  })


f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef(det_pois_study_global[[best]])

guesses[best,]
coef(sir)
summary(det_pois_study_global[[best]])$loglik
h<- function(x) {summary(det_pois_study_global[[x]])$convergence}
conv <- apply(x,2,FUN=h)
conv
length(which(conv==0))
which(conv==0)
g <- function(x) {coef(det_pois_study_global[[x]])}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)

#re-use these values as staring values for another try 
coef <- apply(x,2,FUN=g)
g<-t(coef[,which(conv==0)])



  
  t_global <- system.time({
    guesses <- as.data.frame(g)
    det_pois_study_global <- foreach(guess=iter(guesses,"row"),.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      traj.match(study,
                 start=unlist(guess),
                 est=c("beta1","beta2","beta3","beta11","phi"),
                 method="Nelder-Mead",reltol=1e-5,maxit=2000, transform=TRUE)
    }
  })


n<- ncol(coef[,which(conv==0)])
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle <-coef(det_pois_study_global[[best]])

guesses[best,]
coef(sir)
loglik_mle <- summary(det_pois_study_global[[best]])$loglik
h<- function(x) {summary(det_pois_study_global[[x]])$convergence}
conv <- apply(x,2,FUN=h)
conv
length(which(conv==0))
which(conv==0)
g <- function(x) {coef(det_pois_study_global[[x]])}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
coef <- apply(x,2,FUN=g)


results_global_all <- data.frame(logLik=liks,beta1=coef["beta1",],beta2=coef["beta2",],beta3=coef["beta3",],
                                 beta11=coef["beta11",],phi=coef["phi",])
results_global <- results_global_all[which(conv==0),]
guesses_con <- guesses[which(conv==0),] 
quan<-0.1
# starting values and corresponding estimates which are in a 90 % quantile of the loglik
all <- ldply( list(guess=guesses_con[which(results_global$logLik > quantile(results_global$logLik,probs = quan)),], 
                   results_global=subset(results_global, logLik > quantile(results_global$logLik,probs = quan)) ), .id="type")
plot(liks[which(conv==0)])
pairs(~logLik+beta1+beta2+beta3+beta11+phi,data=all, col=ifelse(all$type=="guess", grey(0.5), "red"), pch=16)
sum(liks < quantile(results_global$logLik,probs = quan))



#plot the results
sim <- trajectory(det_pois_study_global[[best]],as.data.frame=TRUE)
periods <- nrow(dat[-1,])/52
axis.spots <- (0:(periods))*52+2;  
axis.labels <- as.character((2001):(2001+periods)); 
nfu<-c(qpois(0.975,sim$H1)[-1],qpois(0.975,sim$H2)[-1],qpois(0.975,sim$H3)[-1])
nfl<-c(qpois(0.025,sim$H1)[-1],qpois(0.025,sim$H2)[-1],qpois(0.025,sim$H3)[-1])
df2 <- data.frame(lower= nfl, upper=nfu)
df2$age <- factor(rep(c("children", "adult", "elderly"), each = nrow(dat)-1),
                  levels = c("children", "adult", "elderly"))
df2$time <- rep(seq.int( nrow(dat)-1), 3)
df2$type <- "95% prediction interval"
df2$type <- factor(df2$type, levels = c("95% prediction interval"))
df2$labelmedian <- "Mean"
df2$labelmedian <- factor(df2$labelmedian, levels = c("Mean"))
dfcases <-melt(sim_study[-1,], id.vars = c("time"), variable.name = "cases")
df2$cases <- dfcases$value
df2$mean <- c(sim$H1[-1],sim$H2[-1],sim$H3[-1])
df2$labelcases <- "Data"
df2$labelcases <- factor(df2$labelcases, levels = c("Data"))


#pdf(file="~/Dropbox/AAPAPER/Latex/det_pois_study.pdf", width=17, height=4)
ggplot(df2) + 
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill=type), alpha = 0.6) + 
  geom_line(aes(x = time, y = cases, color = labelcases)) +
  geom_line(aes(x = time, y = mean, color = labelmedian)) +
  facet_wrap( ~age, ncol=1, scales =  "free_y") + 
  scale_x_continuous( breaks = axis.spots,labels = axis.labels)+ ylab("Number of new cases")+ 
  xlab("Time(weeks)")+  labs(color="")+theme_bw()
#dev.off()


#check  in and out of 95% confidence interval
a <- which(df2$lower>df2$cases)
b <- which(df2$upper<df2$cases)
length(a)
length(b)
out.of.bounds_mle <- (length(a)+length(b))/(3*(nrow(dat)-1))
out.of.bounds_mle
coef_true <- coef(sir)
#save(coef_true,out.of.bounds_mle, coef_mle, loglik_mle,file = "~/Dropbox/AAPAPER/R_code/traj.match/det_pois_STUDY/res_mle_study_pois.rda")


out.of.bounds_mle <- (length(a)+length(b))/(3*(nrow(sim_study)-1))
aic<- 2*5- 2*loglik_mle
res<- c(coef_mle,loglik=loglik_mle,oob = out.of.bounds_mle,aic=aic)

coef_list[[i]] <-  res
}

coef_fr <- t(sapply(coef_list, c))

save(coef_fr,file = paste("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/coef_traj_pois",N,".rda"))
#plot this as a boxplot (maybe has to be a matrix though!)






