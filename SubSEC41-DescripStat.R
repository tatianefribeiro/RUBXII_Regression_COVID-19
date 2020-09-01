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
setwd("Please, define here your directory")

# COVID-19 data set
data.set <- read_csv("COVID19_US_states.csv")

data.set.select <- data.set%>%   
  dplyr::select(
    MR,
    PD,
    HDI,
    GINI,
    BEDS,
    PR,
    SR,
    AT,
    MA
  )
M = data.set.select


###############################
# Descriptive Statistic table
###############################
# Response for TIME = 30,60,90, 
# and 120 days
MR_30_data <-  data.set%>%dplyr::filter(TIME==30)%>%
  dplyr::select(MR)
MR_30 = MR_30_data$MR

MR_60_data <-  data.set%>%dplyr::filter(TIME==60)%>%
  dplyr::select(MR)
MR_60 = MR_60_data$MR

MR_90_data <-  data.set%>%dplyr::filter(TIME==90)%>%
  dplyr::select(MR)
MR_90 = MR_90_data$MR

MR_120_data <-  data.set%>%dplyr::filter(TIME==120)%>%
  dplyr::select(MR)
MR_120 = MR_120_data$MR

# Covariates without dummies
cov_data <-  M%>%dplyr::select(-MR)

PD = cov_data$PD 
HDI = cov_data$HDI
GINI = cov_data$GINI
BEDS = cov_data$BEDS
PR = cov_data$PR
SR = cov_data$SR
AT = cov_data$AT
MA = cov_data$MA

# Table
sum_qt_cov <- matrix(NA, 12, 7)
rownames(sum_qt_cov) <- c("MR_30","MR_60","MR_90", "MR_120", 
                          "PD",
                          "HDI",
                          "GINI",
                          "BEDS",
                          "PR",
                          "SR",
                          "AT",
                          "MA")
colnames(sum_qt_cov) <- c("Mean","Median","Skewness",
                          "Kurtosis","Min.","Max.","CV")
sum_qt_cov[1,] <- c(mean(MR_30),median(MR_30),
                    skewness(MR_30)[1],kurtosis(MR_30)[1],
                    min(MR_30),max(MR_30),
                    sd(MR_30)/mean(MR_30)*100)

sum_qt_cov[2,] <- c(mean(MR_60),median(MR_60),
                    skewness(MR_60)[1],kurtosis(MR_60)[1],
                    min(MR_60),max(MR_60),
                    sd(MR_60)/mean(MR_60)*100)

sum_qt_cov[3,] <- c(mean(MR_90),median(MR_90), 
                    skewness(MR_90)[1],kurtosis(MR_90)[1],
                    min(MR_90),max(MR_90),
                    sd(MR_90)/mean(MR_90)*100)

sum_qt_cov[4,] <- c(mean(MR_120),median(MR_120), 
                    skewness(MR_120)[1],kurtosis(MR_120)[1],
                    min(MR_120),max(MR_120),
                    sd(MR_120)/mean(MR_120)*100)

sum_qt_cov[5,] <- c(mean(PD),median(PD), 
                    skewness(PD)[1],kurtosis(PD)[1],
                    min(PD),max(PD),
                    sd(PD)/mean(PD)*100)

sum_qt_cov[6,] <- c(mean(HDI),median(HDI), 
                     skewness(HDI)[1],kurtosis(HDI)[1],
                     min(HDI),max(HDI),
                     sd(HDI)/mean(HDI)*100)

sum_qt_cov[7,] <- c(mean(GINI),median(GINI), 
                    skewness(GINI)[1],kurtosis(GINI)[1],
                    min(GINI),max(GINI),
                    sd(GINI)/mean(GINI)*100)

sum_qt_cov[8,] <- c(mean(BEDS),median(BEDS), 
                    skewness(BEDS)[1],kurtosis(BEDS)[1],
                    min(BEDS),max(BEDS),
                    sd(BEDS)/mean(BEDS)*100)

sum_qt_cov[9,] <- c(mean(PR),median(PR), 
                    skewness(PR)[1],kurtosis(PR)[1],
                    min(PR),max(PR),
                    sd(PR)/mean(PR)*100)

sum_qt_cov[10,] <- c(mean(SR),median(SR), 
                    skewness(SR)[1],kurtosis(SR)[1],
                    min(SR),max(SR),
                    sd(SR)/mean(SR)*100)

sum_qt_cov[11,] <- c(mean(AT),median(AT), 
                     skewness(AT)[1],kurtosis(AT)[1],
                     min(AT),max(AT),
                     sd(AT)/mean(AT)*100)

sum_qt_cov[12,] <- c(mean(MA),median(MA), 
                    skewness(MA)[1],kurtosis(MA)[1],
                    min(MA),max(MA),
                    sd(MA)/mean(MA)*100)

round(sum_qt_cov,4)

######################
# Response histogram
#####################
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
par(mfrow=c(2,4))
plot(MR~PD, pch = 20)
plot(MR~HDI, pch = 20)
plot(MR~GINI, pch = 20)
plot(MR~BEDS, pch = 20)
plot(MR~PR, pch = 20)
plot(MR~SR, pch = 20)
plot(MR~AT, pch = 20)
plot(MR~MA, pch = 20)

###################
#Correlation tests
##################
data_cor <- data.set.select
rcorr(as.matrix(data_cor), type = c("spearman"))[3]
stargazer(rcorr(as.matrix(data_cor), type = c("spearman"))[3],digits=4)


