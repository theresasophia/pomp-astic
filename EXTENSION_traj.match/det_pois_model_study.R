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
registerDoParallel()
stopifnot(packageVersion("pomp") >= "1.7")

load("det_pois_model_data.rda")
################################################################################################
################ STUDY #########################################################################
################################################################################################
set.seed(1234567891)
#simulation study
coef(sir)
sim <- trajectory(sir,as.data.frame=TRUE)
data.frame(time=seq(1:416),cases1=rpois(416,sim$H1[-1]),cases2=rpois(416,sim$H2[-1]),cases3=rpois(416,sim$H3[-1]))%>%
  rbind(data.frame(time=0,cases1=NA,cases2=NA,cases3=NA)) %>%
  arrange(time) -> sim_study


 write.table(sim_study[-1,], file = "~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/rota_study10.txt" )


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


save(study,sir_fixed_params,sim_study,file = "det_pois_model_study10.rda")
