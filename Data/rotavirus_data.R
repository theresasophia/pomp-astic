###############################################
# In this code we calculate the weekly number of new rotavirus cases in Germany
# between the years 2001-2008 stratified by age, namely age groups 0-4, 5-59 and 60+, 
# scaled up the the underreporting rates as inferred in Weidemann et al. 2013. 
# The data loaded is available at the following Github repository:
# https://github.com/weidemannf/thesis_code/tree/master/Mode_inference_paperI
# Author: Theresa Stocks
# Date: June 28, 2017
#################################################




# init
rm(list = ls())
#setwd("~/Dropbox/two-age-strata/Felix_data/thesis_code-master/Mode_inference_paperI")

EFSdata <- function(n){
  EFSdat<-numeric(0)
  for (i in 1:n){
    direct<-paste("RotadataEFS200",i,".txt")
    rotadata <- read.table(file=direct)
    for (j in 1:min(52,length(rotadata))){EFSdat<-cbind(EFSdat,rotadata[[j]])}
  }
  return(EFSdat)
}
# number of age groups 1:5 children, 6:8 adults, 9:10 elderly 
nrow(EFSdata(8))


# we consider data form 2001-2008
year <- 8
child <- seq(1,5,by=1)
adult <- seq(6,8,by=1)
elderly <- seq(9,10, by=1)

EFS <- matrix(nrow=3,ncol= ncol(EFSdata(8)))
EFS[1,] <- colSums(EFSdata(year)[child,])
EFS[2,]<- colSums(EFSdata(year)[adult,])
EFS[3,] <- colSums(EFSdata(year)[elderly,])
EFS
#check if all data is used

e <- colSums(EFSdata(8))
e_sum <- colSums(EFS)
e-e_sum

WFSdata <- function(n){
  WFSdat<-numeric(0)
  for (i in 1:n){
    direct<-paste("RotadataWFS200",i,".txt")
    rotadata <- read.table(file=direct)
    for (j in 1:min(52,length(rotadata))){WFSdat<-cbind(WFSdat,rotadata[[j]])}
  }
  return(WFSdat)
}

WFS <- matrix(nrow=3,ncol= ncol(EFSdata(8)))
WFS[1,] <- colSums(WFSdata(year)[child,])
WFS[2,] <- colSums(WFSdata(year)[adult,])
WFS[3,] <- colSums(WFSdata(year)[elderly,])
WFS

#check if all data is used
w <- colSums(WFSdata(8))
w_sum <- colSums(WFS)
w-w_sum

# change of reporting behaviour in beginning of 2005, so from beg 2001- end 2004 
# constant and from beg 2005- end 2008
years_till_change <- 4
time_unit <- 52
time_change <- years_till_change* time_unit
before <- seq(1,time_change,by=1)
after <- seq(time_change+1, ncol(EFSdata(8)), by=1)

# aveaged posterior density of fitted underreportings rates from Weidemann et al. 2013 before and after 2005
EFS_rep_before <- 0.19
EFS_rep_after <- 0.241
WFS_rep_before <- 0.043
WFS_rep_after <- 0.063

# cases of children without underreporting
child_no_rep <- matrix(nrow=2,ncol= ncol(EFSdata(8)))
child_no_rep[1,] <- round(c(EFS[1,][before]/EFS_rep_before,EFS[1,][after]/EFS_rep_after))
child_no_rep[2,]<- round(c(WFS[1,][before]/WFS_rep_before,WFS[1,][after]/WFS_rep_after))
colSums(child_no_rep)

# cases of adult without underreporting
adult_no_rep <- matrix(nrow=2,ncol= ncol(EFSdata(8)))
adult_no_rep[1,] <- round(c(EFS[2,][before]/EFS_rep_before,EFS[2,][after]/EFS_rep_after))
adult_no_rep[2,] <- round(c(WFS[2,][before]/WFS_rep_before,WFS[2,][after]/WFS_rep_after))
colSums(adult_no_rep)

# cases of elderly without underreporting
elderly_no_rep <- matrix(nrow=2,ncol= ncol(EFSdata(8)))
elderly_no_rep[1,] <- round(c(EFS[3,][before]/EFS_rep_before,EFS[3,][after]/EFS_rep_after))
elderly_no_rep[2,] <- round(c(WFS[3,][before]/WFS_rep_before,WFS[3,][after]/WFS_rep_after))
colSums(elderly_no_rep)


rotavirus <- data.frame(time = seq(1:ncol(EFSdata(8))),
                  cases1 = colSums(child_no_rep),
                  cases2 = colSums(adult_no_rep),
                  cases3 = colSums(elderly_no_rep))


rotavirus$cases1 +  rotavirus$cases2 + rotavirus$cases3 -> total 
total  <- data.frame(time = seq(1:ncol(EFSdata(8))),
                        cases = total,
                     cases1 = colSums(child_no_rep),
                     cases2 = colSums(adult_no_rep),
                     cases3 = colSums(elderly_no_rep))
plot(total, type='l')
# plot the data
par(mfrow=c(1,1))
matplot(rotavirus$time,rotavirus[,paste0("cases",1:3)],type='l',lty=1,lwd=2,xlab= "Time", ylab="Cases")
legend(300,36000 , c("Children","Adults","Elderly"),bty="n", lty=c(1,1,1),lwd=c(2,2,2),
       col=c("black", "red", "green"))

write.table(rotavirus,file="rotavirus.txt")
write.table(total,file="total_rota.txt")
