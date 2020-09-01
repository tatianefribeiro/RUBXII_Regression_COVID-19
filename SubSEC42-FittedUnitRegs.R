###################################################################################
# PAPER: A new regression model for the COVID-19 mortality rate in the U.S. States
# SUBSECTION: 4.2. Fitted regressions
# GOAL: Fiting the RUBXII, KW, UW, BETA, and SIMPLEX regression models
# AUTHOR: Tatiane Fontana Ribeiro
# LAST UPDATE: August 27, 2020
###################################################################################
rm(list = ls())
library(tidyverse)
set.seed(2020) 

##Directory
setwd("Please, define here your directory")

#####################              AUXILIARY FUNCTIONS          ############################### 
source("unit_regressions_fit.R") # It fits a unit reg model, provides LOOCV measure, and plots.
source("RUBXII_Kw_UW_functions.R")    #link function and others.

#COVID-19 data set
data.set <- read_csv("COVID19_US_states.csv")

View(data.set)
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
    HDI,
    GINI,
    BEDS,
    PR,
    SR,
    AT,
    T60,
    T90,
    T120
  )

# Covariates matrix
X <- model.matrix(MR~.,data = data.set.select)  

# FITTED RUBXII REGRESSION MODEL
RUBXII_reg = UnitReg.fit(z,X,regression = "RUBXII")
round(RUBXII_reg$summary,4) 
RUBXII_reg$diagnosticMEASURES
#(res=round(LOOCV.unit_reg(z,X, regression = "RUBXII"),4))

# MAE/mean(MR)
mae = as.numeric(res)
(comp = round(mae/mean(z),4))

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

# Plots
plots_quantres(RUBXII_reg$residuals,RUBXII_reg$q.fv,n)

#*****************************************************************************************
#                             Kw Regression Model
#*****************************************************************************************
#Selected covariates
data.set.select <- data.set%>%   
  dplyr::select(
    MR,
    PD,
    GINI,
    SR,
    AT,
    T60,
    T90,
    T120
  )

# Covariates matrix
X <- model.matrix(MR~.,data = data.set.select)  

# FITTED RUBXII REGRESSION MODEL
Kw_reg = UnitReg.fit(z,X,regression = "Kw")
round(Kw_reg$summary,4)                       
Kw_reg$diagnosticMEASURES
#(res_Kw=round(LOOCV.unit_reg(z,X, regression = "Kw"),4))

# MAE/mean(MR)
mae_Kw = as.numeric(res_Kw)
(comp_Kw = round(mae_Kw/mean(z),4))

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

# Plots
plots_quantres(Kw_reg$residuals,Kw_reg$mu.fv,n)


#*****************************************************************************************
#                             UW Regression Model
#*****************************************************************************************
#Selected covariates
data.set.select <- data.set%>%   
  dplyr::select(
    MR,
    PD,
    GINI,
    AT,
    T120
  )

# Covariates matrix
X <- model.matrix(MR~.,data = data.set.select)  

# FITTED RUBXII REGRESSION MODEL
UW_reg = UnitReg.fit(z,X,regression = "UW")
round(UW_reg$summary,4)                       
UW_reg$diagnosticMEASURES
(res_UW=round(LOOCV.unit_reg(z,X, regression = "UW"),4))

# MAE/mean(MR)
mae_UW = as.numeric(res_UW)
(comp_UW = round(mae_UW/mean(z),4))

#RESET-type test
unit_reg <- UnitReg.fit(z,X,regression = "UW")
q_hat2 <- (unit_reg$q.fv)^2
q_hat3 <- (unit_reg$q.fv)^3
data_aug_mod <- data.frame(z,X[,-1],omega_hat2,omega_hat3)
X_aug <- model.matrix(z~.,data = data_aug_mod)  
X <- X_aug
mod_aug <- UnitReg.fit(z,X,regression = "UW")
ell_rest <- unit_reg$max_ell
ell_unr <- mod_aug$max_ell
(w_LR <- 2*(ell_unr-ell_rest))
(p_value_LR <- 1-pchisq(w_LR,2)) 

# Plots
plots_quantres(UW_reg$residuals,UW_reg$mu.fv,n)

#*****************************************************************************************
#                             BETA Regression Model
#*****************************************************************************************
#Selected covariates
data.set.select <- data.set%>%   
  dplyr::select(
    MR,
    PD,
    HDI,
    PR,
    AT,
    T60,
    T90,
    T120
  )

# Covariates matrix
X <- model.matrix(MR~.,data = data.set.select)  

# FITTED RUBXII REGRESSION MODEL
BETA_reg = UnitReg.fit(z,X,regression = "BETA")
round(BETA_reg$summary,4)                          #!!!!!!!! NA, não sei como resolver.
BETA_reg$diagnosticMEASURES
#(res_BETA=round(LOOCV.unit_reg(z,X, regression = "BETA"),4))

# MAE/mean(MR)
mae_BETA = as.numeric(res_BETA)
(comp_BETA = round(mae_BETA/mean(z),4))

#RESET-type test
unit_reg <- UnitReg.fit(z,X,regression = "BETA")
mu_hat2 <- (unit_reg$mu.fv)^2
mu_hat3 <- (unit_reg$mu.fv)^3
data_aug_mod <- data.frame(z,X[,-1],mu_hat2,mu_hat3)
X_aug <- model.matrix(z~.,data = data_aug_mod)  
X <- X_aug
mod_aug <- UnitReg.fit(z,X,regression = "BETA")
ell_rest <- unit_reg$max_ell
ell_unr <- mod_aug$max_ell
(w_LR <- 2*(ell_unr-ell_rest))
(p_value_LR <- 1-pchisq(w_LR,2)) 

# Plots
plots_quantres(BETA_reg$residuals,BETA_reg$mu.fv,n)


#*****************************************************************************************
#                             SIMPLEX Regression Model
#*****************************************************************************************
#Selected covariates
data.set.select <- data.set%>%   
  dplyr::select(
    MR,
    GINI,
    BEDS,
    MA,
    T60,
    T90,
    T120
  )

# Covariates matrix
X <- model.matrix(MR~.,data = data.set.select)  

# FITTED RUBXII REGRESSION MODEL
SIMPLEX_reg = UnitReg.fit(z,X,regression = "SIMPLEX")
round(SIMPLEX_reg$summary,4)                        
SIMPLEX_reg$diagnosticMEASURES
#(res_SIMPLEX=round(LOOCV.unit_reg(z,X, regression = "SIMPLEX"),4))

# MAE/mean(MR)
mae_SIMPLEX = as.numeric(res_SIMPLEX)
(comp_SIMPLEX = round(mae_SIMPLEX/mean(z),4))

#RESET-type test
unit_reg <- UnitReg.fit(z,X,regression = "SIMPLEX")
mu_hat2 <- (unit_reg$mu.fv)^2
mu_hat3 <- (unit_reg$mu.fv)^3
data_aug_mod <- data.frame(z,X[,-1],mu_hat2,mu_hat3)
X_aug <- model.matrix(z~.,data = data_aug_mod)  
X <- X_aug
mod_aug <- UnitReg.fit(z,X,regression = "SIMPLEX")
ell_rest <- unit_reg$max_ell
ell_unr <- mod_aug$max_ell
(w_LR <- 2*(ell_unr-ell_rest))
(p_value_LR <- 1-pchisq(w_LR,2)) 

# Plots
plots_quantres(SIMPLEX_reg$residuals,SIMPLEX_reg$mu.fv,n)

