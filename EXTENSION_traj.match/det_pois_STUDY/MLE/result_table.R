
load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_STUDY/MLE/det_pois_MLE_study_second_1-%d.rda")
n<- 100
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle1 <-coef(det_pois_study_global[[best]])

load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_STUDY/MLE/det_pois_MLE_study_second_2-%d.rda")
n<- 100
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle2 <-coef(det_pois_study_global[[best]])


load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_STUDY/MLE/det_pois_MLE_study_second_3-%d.rda")
n<- 100
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle3 <-coef(det_pois_study_global[[best]])

load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_STUDY/MLE/det_pois_MLE_study_second_4-%d.rda")
n<- 100
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle4 <-coef(det_pois_study_global[[best]])

load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_STUDY/MLE/det_pois_MLE_study_second-%d.rda")
n<- 100
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle5 <-coef(det_pois_study_global[[best]])

load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_STUDY/MLE/det_pois_MLE_study_second_6-%d.rda")
n<- 100
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle6 <-coef(det_pois_study_global[[best]])

load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_STUDY/MLE/det_pois_MLE_study_second_7-%d.rda")
n<- 100
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle7 <-coef(det_pois_study_global[[best]])

load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_STUDY/MLE/det_pois_MLE_study_second_8-%d.rda")
n<- 100
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle8 <-coef(det_pois_study_global[[best]])

load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_STUDY/MLE/det_pois_MLE_study_second_9-%d.rda")
n<- 100
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle9 <-coef(det_pois_study_global[[best]])

load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_pois_STUDY/MLE/det_pois_MLE_study_second_10-%d.rda")
n<- 100
f <- function(x) {summary(det_pois_study_global[[x]])$loglik}
x <- matrix(c(seq(from=1,to=n,by=1)),nrow=1,ncol=n)
liks <- apply(x,2,FUN=f)
best <- which.max(liks)
coef_mle10 <-coef(det_pois_study_global[[best]])


table<- round(rbind(coef_mle1,coef_mle2,coef_mle3,coef_mle4,coef_mle5,coef_mle6,coef_mle7,coef_mle8,coef_mle9,coef_mle10),4)
table_ext <- rbind(table, mean=colMeans(table))

library(xtable)
xtable(as.data.frame(table_ext)[,1:5],digits=4)
