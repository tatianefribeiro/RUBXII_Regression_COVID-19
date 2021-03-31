###################################################################################
# PAPER: A new regression model for the COVID-19 mortality rate in the U.S. States
# SUBSECTION: 4.2. Fitted regressions
# GOAL: Fiting the RUBXII, KW, and UW regression models
# AUTHOR: Tatiane Fontana Ribeiro
# LAST UPDATE: March 10, 2021
###################################################################################
rm(list = ls())
library(tidyverse)

##Directory
setwd("/media/tatiane/239fb286-e175-4578-b0d3-48c750946446/RUBXII_regres_c19/GitHub_C19")


#####################              AUXILIARY FUNCTIONS          ############################### 
source("unit_regressions_fit.R") # It fits a unit reg model, provides LOOCV measure, and plots.
source("RUBXII_Kw_UW_functions.R")    #link function and others.

#COVID-19 data set
data.set <- read_csv("C19_MR_cov_USAbyState.csv")

#View(data.set)
# Response and sample size
z <- data.set$MR
n <- length(z)

#*****************************************************************************************
#                             RUBXII Regression Model
#*****************************************************************************************
#Selected covariates
data.set.select <- data.set%>%  
  dplyr::select(
    MR,
    PD, 
    GINI,
    BEDS,
    SR,
   # LE,
  # PR,
    T90,
    T180)
  
    # Covariates matrix
X <- model.matrix(MR~.,data = data.set.select)  

# FITTED RUBXII REGRESSION MODEL
RUBXII_reg = UnitReg.fit(z,X,regression = "RUBXII")
round(RUBXII_reg$summary,4) 
RUBXII_reg$diagnosticMEASURES


# FITTED Kw REGRESSION MODEL
Kw_reg = UnitReg.fit(z,X,regression = "Kw")
round(Kw_reg$summary,4)                       
Kw_reg$diagnosticMEASURES

#RESET-type test
unit_reg <- UnitReg.fit(z,X,regression = "RUBXII")
q_hat2 <- (unit_reg$q.fv)^2
q_hat3 <- (unit_reg$q.fv)^3
data_aug_mod <- data.frame(z,X[,-1],q_hat2,q_hat3)
X_aug <- model.matrix(z~.,data = data_aug_mod)  
X <- X_aug
mod_aug <- UnitReg.fit(z,X,regression = "RUBXII")
ell_rest <- unit_reg$max_ell
ell_unr <- mod_aug$max_ell
(w_LR <- 2*(ell_unr-ell_rest))
(p_value_LR <- 1-pchisq(w_LR,2)) 

#Goodness-of-fit measures
RUBXII_goodfit <- matrix(NA,nrow = 1,ncol = 4)
colnames(RUBXII_goodfit) <- c("LL","R_G","p-val(AD)","p-val(RES)")
RUBXII_goodfit[1,] <- c(RUBXII_reg$max_ell, RUBXII_reg$diagnosticMEASURES[2], 
                        RUBXII_reg$diagnosticMEASURES[,3],
                        as.numeric(p_value_LR))

# Plots
plots_quantres(RUBXII_reg$residuals,RUBXII_reg$q.fv,n)

#*****************************************************************************************
#                             Kw Regression Model
#*****************************************************************************************
# Covariates matrix
X <- model.matrix(MR~.,data = data.set.select)  

# FITTED RUBXII REGRESSION MODEL
Kw_reg = UnitReg.fit(z,X,regression = "Kw")
round(Kw_reg$summary,4)                       
Kw_reg$diagnosticMEASURES


#RESET-type test
unit_reg <- UnitReg.fit(z,X,regression = "Kw")
omega_hat2 <- (unit_reg$q.fv)^2
omega_hat3 <- (unit_reg$q.fv)^3
data_aug_mod <- data.frame(z,X[,-1],omega_hat2,omega_hat3)
X_aug <- model.matrix(z~.,data = data_aug_mod)  
X <- X_aug
mod_aug <- UnitReg.fit(z,X,regression = "Kw")
ell_rest <- unit_reg$max_ell
ell_unr <- mod_aug$max_ell
(w_LR <- 2*(ell_unr-ell_rest))
(p_value_LR <- 1-pchisq(w_LR,2)) 


#Goodness-of-fit measures
Kw_goodfit <- matrix(NA,nrow = 1,ncol = 4)
colnames(Kw_goodfit) <- c("LL","R_G","p-val(AD)","p-val(RES)")
Kw_goodfit[1,] <- c(Kw_reg$max_ell, Kw_reg$diagnosticMEASURES[2], 
                        Kw_reg$diagnosticMEASURES[,3],
                        as.numeric(p_value_LR))

# Plots
plots_quantres(Kw_reg$residuals,Kw_reg$mu.fv,n)


#*****************************************************************************************
#                             UW Regression Model
#*****************************************************************************************
# Covariates matrix
X <- model.matrix(MR~.,data = data.set.select)  

# FITTED UW REGRESSION MODEL
UW_reg = UnitReg.fit(z,X,regression = "UW")
round(UW_reg$summary,4)                       
UW_reg$diagnosticMEASURES

#RESET-type test
unit_reg <- UnitReg.fit(z,X,regression = "UW")
omega_hat2 <- (unit_reg$mu.fv)^2
omega_hat3 <- (unit_reg$mu.fv)^3
data_aug_mod <- data.frame(z,X[,-1],omega_hat2,omega_hat3)
X_aug <- model.matrix(z~.,data = data_aug_mod)  
X <- X_aug
mod_aug <- UnitReg.fit(z,X,regression = "UW")
ell_rest <- unit_reg$max_ell
ell_unr <- mod_aug$max_ell
(w_LR <- 2*(ell_unr-ell_rest))
(p_value_LR <- 1-pchisq(w_LR,2)) 

#Goodness-of-fit measures
UW_goodfit <- matrix(NA,nrow = 1,ncol = 4)
colnames(UW_goodfit) <- c("LL","R_G","p-val(AD)","p-val(RES)")
UW_goodfit[1,] <- c(UW_reg$max_ell, UW_reg$diagnosticMEASURES[2], 
                    UW_reg$diagnosticMEASURES[,3],
                    as.numeric(p_value_LR))

# Plots
plots_quantres(UW_reg$residuals,UW_reg$mu.fv,n)


#For tex (Tables)
cbind(RUBXII_reg$summary[,c(1,4)],
      Kw_reg$summary[,c(1,4)],
      UW_reg$summary[,c(1,4)])

stargazer::stargazer(cbind(RUBXII_reg$summary[,c(1,4)],
                           Kw_reg$summary[,c(1,4)],
                           UW_reg$summary[,c(1,4)]), digits = 4
)

round(rbind(RUBXII_goodfit,Kw_goodfit,UW_goodfit),4)
stargazer::stargazer(rbind(RUBXII_goodfit,Kw_goodfit,UW_goodfit),digits = 4)

# For tex (plots)
setwd("/media/tatiane/239fb286-e175-4578-b0d3-48c750946446/RUBXII_regres_c19/Sub_200929_AMA-stix/ama/figures")
w1<-7
h2<-3
postscript(file = "residuals_plot.eps",horizontal=F,
           pointsize = 14,
           paper="special",
           width = w1, height = h2,
           family = "Times")
{
  par(mfrow=c(1,3))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  
  quant_residuals = RUBXII_reg$residuals
  qqnorm(quant_residuals,pch=1,frame=T, main = "RUBXII",
         make.plot = T, lwd=1)
  qqline(quant_residuals,col="red")
  
  quant_residuals = Kw_reg$residuals
  qqnorm(quant_residuals,pch=1,frame=T, main = "Kw",
         make.plot = T, lwd=1)
  qqline(quant_residuals,col="red")
  
  quant_residuals = UW_reg$residuals
  qqnorm(quant_residuals,pch=1,frame=T, main = "UW",
         make.plot = T, lwd=1)
  qqline(quant_residuals,col="red")
}
dev.off()

