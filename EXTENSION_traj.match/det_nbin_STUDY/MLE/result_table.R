load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_STUDY/res_mle_study_neg1.rda")
coef_mle1 <- coef_mle
load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_STUDY/res_mle_study_neg2.rda")
coef_mle2 <- coef_mle
load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_STUDY/res_mle_study_neg3.rda")
coef_mle3 <- coef_mle
load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_STUDY/res_mle_study_neg4.rda")
coef_mle4 <- coef_mle
load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_STUDY/res_mle_study_neg.rda")
coef_mle5 <- coef_mle
load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_STUDY/res_mle_study_neg6.rda")
coef_mle6 <- coef_mle
load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_STUDY/res_mle_study_neg7.rda")
coef_mle7 <- coef_mle
load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_STUDY/res_mle_study_neg8.rda")
coef_mle8 <- coef_mle
load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_STUDY/res_mle_study_neg9.rda")
coef_mle9 <- coef_mle
load("~/Dropbox/AAPAPER/R_code/traj.match/det_nbin_STUDY/res_mle_study_neg10.rda")
coef_mle10 <- coef_mle



table<- round(rbind(coef_mle1,coef_mle2,coef_mle3,coef_mle4,coef_mle5,coef_mle6,coef_mle7,coef_mle8,coef_mle9,coef_mle10),4)
table_ext <- round(rbind(table, mean=colMeans(table)),3)
table_ext[,1:6]
library(xtable)
xtable(as.data.frame(table_ext)[,1:6],digits=4)
