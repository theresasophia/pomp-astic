
########creating a boxplot for the 100 samples of the traj_nbin simulation study matching###



rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows
library(pomp)
library(magrittr)
library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
library(foreach)
library(doParallel)
library(plyr)
library("tidyr")
library("ggplot2")
library(dplyr)

library("grDevices")
registerDoParallel()
stopifnot(packageVersion("pomp") >= "1.7")
########creating a boxplot for the 100 samples of the traj_pois simulation study matching###



load(file="~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/coef_traj_nbin 100 .rda")
load("~/Dropbox/AAPAPER/R_code/EXTENSION_traj.match/det_nbin_model_data.rda")
true <- coef(sir)
coef_true <- as.data.frame(true, row.names = NULL)
coef_true$variable <- names(true)

coef_true %>%
  mutate(variable = factor(variable))->mle_true

class(coef_fr)
coef_fr%>%
  subset(select=c("beta1" , "beta2" , "beta3",  "beta11", "phi" ,   "od")) %>%
  melt(id=c("sample"))  %>%
  mutate(Var2= factor(Var2))-> coef_fr

new_data <- merge(coef_fr,mle_true,  by.x = "Var2", by.y = "variable")


cairo_ps(file="~/Dropbox/AAPAPER/Revision_Biostatistics/ShareLaTeX/Revision_biostatistics/boxplot_nbin.eps", width=17, height=10)
new_data%>%
  mutate(Var2 = recode(Var2, phi = "phi")) %>%
  mutate(Var2 = recode(Var2, beta11 = "rho")) %>%
  mutate(Var2 = recode(Var2, beta1 = "beta[1]")) %>%
  mutate(Var2 = recode(Var2, beta2 = "beta[2]")) %>%
  mutate(Var2 = recode(Var2, beta3 = "beta[3]")) %>%
  mutate(Var2 = recode(Var2, od = "theta")) %>%
ggplot(aes(x=Var2, y=value))+  geom_boxplot() +facet_wrap( ~Var2, ncol=3, scales =  "free",labeller= label_parsed) +scale_x_discrete(breaks=NULL) +
  geom_hline(aes(yintercept = true), color="red")+xlab("") +ylab("Value")+theme(text = element_text(size=30))
dev.off()
