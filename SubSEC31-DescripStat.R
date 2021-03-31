###################################################################################
# PAPER: A new regression model for the COVID-19 mortality rate in the U.S. States
# SUBSECTION: 4.1. Descriptive statistical analysis
# GOAL: Doing a descriptive statistical analysis of the data set
# AUTHOR: Tatiane Fontana Ribeiro
# LAST UPDATE: August 28, 2020
###################################################################################
rm(list = ls())
library(tidyverse)
library(corrplot) 
library(corrgram)
library(fBasics)
library(Hmisc)

##Directory
#setwd("Please, define here your directory")
setwd("/media/tatiane/239fb286-e175-4578-b0d3-48c750946446/RUBXII_regres_c19/GitHub_C19")

# COVID-19 data set
data.set <- read_csv("C19_MR_cov_USAbyState.csv")

View(data.set)
data.set.select <- data.set%>%   
  dplyr::select(
    MR,
    PD, 
    GINI,
    BEDS,
    SR,
    PR,
    LE
  )
M = data.set.select


###############################
# Descriptive Statistic table
###############################
# Response for TIME = 30,90,180, 
# and 180 days
MR_30_data <-  data.set%>%dplyr::filter(TIME==30)%>%
  dplyr::select(MR)
MR_30 = MR_30_data$MR

MR_90_data <-  data.set%>%dplyr::filter(TIME==90)%>%
  dplyr::select(MR)
MR_90 = MR_90_data$MR

MR_180_data <-  data.set%>%dplyr::filter(TIME==180)%>%
  dplyr::select(MR)
MR_180 = MR_180_data$MR

# Covariates without dummies
cov_data <-  M%>%dplyr::select(-MR)

PD = cov_data$PD 
GINI = cov_data$GINI
BEDS = cov_data$BEDS
SR = cov_data$SR
PR = cov_data$PR
LE = cov_data$LE

# Table
sum_qt_cov <- matrix(NA, 9, 7)
rownames(sum_qt_cov) <- c("MR_30","MR_90", "MR_180", 
                          "PD",
                          "GINI",
                          "BEDS",
                          "SR",
                          "PR",
                          "LE")
colnames(sum_qt_cov) <- c("Mean","Median","Skewness",
                          "Kurtosis","Min.","Max.","CV")
sum_qt_cov[1,] <- c(mean(MR_30),median(MR_30),
                    skewness(MR_30)[1],kurtosis(MR_30)[1],
                    min(MR_30),max(MR_30),
                    sd(MR_30)/mean(MR_30)*100)

sum_qt_cov[2,] <- c(mean(MR_90),median(MR_90),
                    skewness(MR_90)[1],kurtosis(MR_90)[1],
                    min(MR_90),max(MR_90),
                    sd(MR_90)/mean(MR_90)*100)

sum_qt_cov[3,] <- c(mean(MR_180),median(MR_180), 
                    skewness(MR_180)[1],kurtosis(MR_180)[1],
                    min(MR_180),max(MR_180),
                    sd(MR_180)/mean(MR_180)*100)

sum_qt_cov[4,] <- c(mean(PD),median(PD), 
                    skewness(PD)[1],kurtosis(PD)[1],
                    min(PD),max(PD),
                    sd(PD)/mean(PD)*100)

sum_qt_cov[5,] <- c(mean(GINI),median(GINI), 
                    skewness(GINI)[1],kurtosis(GINI)[1],
                    min(GINI),max(GINI),
                    sd(GINI)/mean(GINI)*100)

sum_qt_cov[6,] <- c(mean(BEDS),median(BEDS), 
                    skewness(BEDS)[1],kurtosis(BEDS)[1],
                    min(BEDS),max(BEDS),
                    sd(BEDS)/mean(BEDS)*100)

sum_qt_cov[7,] <- c(mean(SR),median(SR), 
                    skewness(SR)[1],kurtosis(SR)[1],
                    min(SR),max(SR),
                    sd(SR)/mean(SR)*100)

sum_qt_cov[8,] <- c(mean(PR),median(PR), 
                    skewness(PR)[1],kurtosis(PR)[1],
                    min(PR),max(PR),
                    sd(PR)/mean(PR)*100)

sum_qt_cov[9,] <- c(mean(LE),median(LE), 
                     skewness(LE)[1],kurtosis(LE)[1],
                     min(LE),max(LE),
                     sd(LE)/mean(LE)*100)

round(sum_qt_cov,4)
stargazer::stargazer(sum_qt_cov,digits = 4)

######################
# Response histogram
#####################
par(mfrow=c(1,2))
hist(data.set$MR,freq = F,breaks = 20, main = "",
     xlab = "MR", col="grey")
boxplot(data.set$MR~data.set$TIME, col =  "grey", xlab = "TIME",ylab = "MR")

####################
#Correlation matrix
####################
par(ask=F)
par(mfrow=c(1,1))
corrplot(corrgram(M, cor.method = "spearman"), method = "circle",type = "upper", 
         addCoef.col = "black", diag = F) #plot matrix

##################
# Dispersion plots
##################
MR = data.set$MR
par(mfrow=c(2,3))
plot(MR~PD, pch = 20)
plot(MR~GINI, pch = 20)
plot(MR~BEDS, pch = 20)
plot(MR~SR, pch = 20)
plot(MR~PR, pch = 20)
plot(MR~LE, pch = 20)

###################
#Correlation tests
##################
data_cor <- data.set.select
rcorr(as.matrix(data_cor), type = c("spearman"))[3]
stargazer::stargazer(rcorr(as.matrix(data_cor), type = c("spearman"))[3],digits=4)


#For tex
setwd("/media/tatiane/239fb286-e175-4578-b0d3-48c750946446/RUBXII_regres_c19/Sub_200929_AMA-stix/ama/figures")
w1<-7
h2<-5
postscript(file = "response_hist.eps",horizontal=F,
           pointsize = 14,
           paper="special",
           width = w1, height = h2,
           family = "Times")
{
  par(mfrow=c(1,2))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  par(mfrow=c(1,2))
  hist(data.set$MR,freq = F,breaks = 20, main = "",
       xlab = "MR", col="grey")
  boxplot(data.set$MR~data.set$TIME, col =  "grey", xlab = "TIME",ylab = "MR")
}
dev.off()


w1<-7
h2<-5
postscript(file = "corr_matrix.eps",horizontal=F,
           pointsize = 14,
           paper="special",
           width = w1, height = h2,
           family = "Times")
{
  
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  par(ask=F)
  
  corrplot(corrgram(M, cor.method = "spearman"), method = "circle",type = "upper", 
           addCoef.col = "black", diag = F) #plot matrix
  
}
dev.off()


w1<-7
h2<-5
postscript(file = "disp_plots.eps",horizontal=F,
           pointsize = 14,
           paper="special",
           width = w1, height = h2,
           family = "Times")
{
  par(mfrow=c(2,3))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  
  MR = data.set$MR
  plot(MR~PD, pch = 20)
  plot(MR~GINI, pch = 20)
  plot(MR~BEDS, pch = 20)
  plot(MR~SR, pch = 20)
  plot(MR~PR, pch = 20)
  plot(MR~LE, pch = 20)
}
dev.off()



