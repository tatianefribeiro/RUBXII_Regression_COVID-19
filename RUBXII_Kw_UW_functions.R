###################################################################################
# PAPER: A new regression model for the COVID-19 mortality rate in the U.S. States
# TITTLE: Reflexive unit Burr XII (RUBXII), Kumaraswamy (KW), and unit Weibull (UW)
#         regression models: some useful functions.
# GOAL: Providing some useful functions for the RUBXII, KW, and UW regression models
#       such as: link function, function to generate occurrences by the inversion 
#       method, log-likelihood, cumulative distribution probability (cdf), and 
#       probability density function (pdf).
# AUTHOR: Tatiane Fontana Ribeiro
# LAST UPDATE: August 27, 2020
###################################################################################

# Logit link function
lfunc <- function(beta_vec,X){  
  qi_logit_link <- 1/(1+exp(-(X%*%beta_vec)))
  return(qi_logit_link)
}

###################################################################################
############  RUBXII regression model - USEFUL FUNCTIONS   ########################
###################################################################################
r_RUBXII <- function(n,q_i,c)   #It generates occurences of Z_i ~ RUBXII (q_i, c)
{
  u = runif(n)
  z_RUBXII = 1-exp(-((1-u)^(log(1+(log(1/(1-q_i)))^c)/log(1-tau))-1)^(1/c)) #qf RUBXII qf
  return(z_RUBXII)
}

l_RUBXIIreg <- function(par){  ## RUBXII log-likelihood
  beta_vec <- par[-length(par)]   
  c <- par[length(par)]            
  q_i <- lfunc(beta_vec,X)  
  
  ell =  sum(
    log(
      c*log(1/(1-tau))*(log(1/(1-z)))^(c-1)*
        (1+(log(1/(1-z)))^c)^(log(1-tau)/log(1+(log(1/(1-q_i)))^c)-1)/  
        ((1-z)*log(1+(log(1/(1-q_i)))^c))    
      )
  )
  
  return(ell)
}

cdf_RUBXII <- function(z,q,c)   ## RUBXII cdf
{
  1-(1+(log(1/(1-z)))^c)^(log(1-tau)/log(1+(log(1/(1-q)))^c))
}

l_RUBXIIreg_aug <- function(par){  # RUBXII log-likelihood (augmented model)  
  beta_vec <- par[-length(par)]  
  c <- par[length(par)]          
  q_i <- lfunc(beta_vec,X_aug) 
  
  ell =  sum(
    log(
      c*log(1/(1-tau))*(log(1/(1-z)))^(c-1)*
        (1+(log(1/(1-z)))^c)^(log(1-tau)/log(1+(log(1/(1-q_i)))^c)-1)/  
        ((1-z)*log(1+(log(1/(1-q_i)))^c))   
    )
  )
  return(ell)
}

#########    NULL MODEL LOG-LIK   ##############
l_RUBXII_noreg <- function(par){  
  c <- par[1]
  q <- par[2]
  
  ell <- sum(
    log(
      c*log(1/(1-tau))*(log(1/(1-z)))^(c-1)*
        (1+(log(1/(1-z)))^c)^(log(1-tau)/log(1+(log(1/(1-q)))^c)-1)/  
        ((1-z)*log(1+(log(1/(1-q)))^c))   
    )
  )
  return(ell)
}


################################################################################
############  Kw regression model - USEFUL FUNCTIONS   #########################
################################################################################
r_Kw <- function(n,w_i,dp)   #It generates occurences of Z_i ~ Kw (w_i, dp)
{
  u = runif(n)
  z_Kw = (1-(1-u)^(log(1-w_i^(1/dp))/log(.5)))^dp #  Kw qf
  return(z_Kw)
}

l_Kwreg <- function(par){ # Kw Log-likelihood 
  beta_vec <- par[-length(par)]  
  dp <- par[length(par)]          
  w_i <- lfunc(beta_vec,X)   
  
  ell =  sum(
    log(
      log(.5)/(dp*log(1-w_i^(1/dp)))*
        z^(1/dp-1)*(1-z^(1/dp))^
        (log(0.5)/log(1-w_i^(1/dp))-1)
    )
  )
  
  return(ell)
}

cdf_Kw <- function(z,w,dp)   # Kw cdf
{
  1-(1-z^(1/dp))^(log(0.5)/log(1-w^(1/dp)))
}

l_Kwreg_aug <- function(par){  # Kw log-likelihood (augmented model)
  beta_vec <- par[-length(par)]  
  dp <- par[length(par)]          
  w_i <- lfunc(beta_vec,X_aug)   
  
  ell =  sum(
    log(
      log(.5)/(dp*log(1-w_i^(1/dp)))*
        z^(1/dp-1)*(1-z^(1/dp))^
        (log(0.5)/log(1-w_i^(1/dp))-1)
    )
  )
  
  return(ell)
}

#########    NULL MODEL LOG-LIK   ##############3
l_Kw_noreg <- function(par){ 
  w <- par[1]
  dp <- par[2]
  
  ell <- sum(
    log(
      log(.5)/(dp*log(1-w^(1/dp)))*
        z^(1/dp-1)*(1-z^(1/dp))^
        (log(0.5)/log(1-w^(1/dp))-1) 
    )
  )
  return(ell)
}


################################################################################
############  UW quantile regression model - USEFUL FUNCTIONS   ################
################################################################################
pdf_UW <- function(z,q,b){  # UW pdf
  b/z*(log(tau)/log(q))*(log(z)/log(q))^(b-1)*
    tau^((log(z)/log(q))^b)
}

cdf_UW <- function(z,q,b)   # UW cdf
{
  tau^((log(z)/log(q))^b)
}

r_UW <- function(n,q,b)    #It generates occurences of Z_i ~ UW (q_i, \beta)
{
  u = runif(n)
  z_UW = exp(log(q)*((log(u)/log(tau))^(1/b)))
  return(z_UW)
}

l_UWreg <- function(par){   # UW log-likelihood
  beta_vec <- par[-length(par)]   
  b <- par[length(par)]           
  
  q_i <- lfunc(beta_vec,X)   
  ell =  sum(log(pdf_UW(z,q_i,b)))
  return(ell)
}



l_UWreg_aug <- function(par){    # UW log-likelihood (augmented model)
  beta_vec <- par[-length(par)]  
  b <- par[length(par)]        
  
  q_i <- lfunc(beta_vec,X_aug) 
  ell =sum(log(pdf_UW(z,q_i,b)))
  return(ell)
}

#########    NULL MODEL LOG-LIK   ##############3
l_UW_noreg <- function(par){  
  q <- par[1]
  b <- par[2]
  
  ell <- sum(log(pdf_UW(z,q,b)))
  return(ell)
}
