########################################################################################
# PAPER: A new regression model for the COVID-19 mortality rate in the U.S. States
# SUBSECTION: 4.2. Fitted regressions
# GOAL: Providing functions to fit a RUBXII, KW, UW, BETA, and SIMPLEX regression model
# AUTHOR: Tatiane Fontana Ribeiro
# LAST UPDATE: August 27, 2020
########################################################################################

#Packages
library(nortest)
library(gamlss)
library(extRemes)

####################################################################################
#                        IT FITS A UNIT REGRESSION MODEL
####################################################################################
tau = 0.5

#Main function
UnitReg.fit <- function(z, X, n = NA, regression = "RUBXII"){ 
  
  source("RUBXII_Kw_UW_functions.R")
  source("UW_reg.R")
  
  if (any(regression == c("RUBXII", "Kw", "UW"))){
    if (regression == "RUBXII"){
      loglik <- l_RUBXIIreg
      cdf <- cdf_RUBXII
      ell_null <- l_RUBXII_noreg
      guess_noreg <- c(1,median(z)) 
    }  
    if (regression == "Kw"){
      loglik <- l_Kwreg
      cdf <- cdf_Kw
      ell_null <- l_Kw_noreg
      guess_noreg <- c(median(z),1) 
    } 
    if (regression == "UW"){
      
      #Fit Beta REG MODEL
      mod_UW <- gamlss(z~X[,-1], 
                         family = "UW", trace= F)
      
      #Quantile residuals
      mu_fv <- as.numeric(mod_UW$mu.fv)
      sigma_fv <- as.numeric(mod_UW$sigma.fv)
      quant_res_UW <- mod_UW$residuals 
      
      #goodness-of-fit measures
      good_fit_meas_mat <- matrix(NA, nrow = 1,ncol = 3)
      colnames(good_fit_meas_mat) <- c("AIC","R-square","p-val(AD)")
      rownames(good_fit_meas_mat) <- c(regression)
      good_fit_meas_mat[1,] <- c( AIC(mod_UW,k=2),
                                  Rsq(mod_UW),
                                  as.numeric(shapiro.test(quant_res_UW)[2]))
      
      
      #summary UW
      beta.coef <- summary(mod_UW)
      mles_UW <- beta.coef[,1]
      se_UW <- beta.coef[,2]
      t_value <- beta.coef[,3]
      p_value_UW <- beta.coef[,4]
      
      K <- ncol(X)+1
  
      #SEs
      summ_UW = matrix(NA,nrow = K, ncol = 4)
      row.names(summ_UW) <- c(names(data.frame(X)),"sigma")
      colnames(summ_UW) <- c("Estimate", "Std. Error", "t value","Pr(>|t|)")
      summ_UW[,1] <- mles_UW
      summ_UW[,2] <- se_UW
      summ_UW[,3] <- t_value
      summ_UW[,4] <- p_value_UW   #Se RODAR por aqui FUNCIONA.
    
      #round(summ_UW, digits = 4)
      
      log.lik.UW = gen.likelihood(mod_UW)
      ell_UW = -log.lik.UW()
      
      results = list(summary = summ_UW, 
                     UW.coeff = mles_UW[1:(K-1)],
                     sigma.coeff = mles_UW[K],
                     mu.fv = mu_fv, 
                     sigma.fv = sigma_fv, 
                     residuals = quant_res_UW,
                     diagnosticMEASURES = good_fit_meas_mat,
                     max_ell = ell_UW,
                     eta_hat = X%*%mles_UW[1:(K-1)]
      )
      
      return(results)
    } 
  }
  else {
    stop(paste(regression, "Regression not available, available regressions are \"RUBXII\", ",
               "\"Kw\", and \"UW\""))
  }
  
  #############################################
  #It fits a RUBXII or Kw or UW regression model
  #############################################
  n <- length(z)
  ### Guess ####
  t =  log(z/(1-z))  #response variable for MQO as initial guess 
  beta_guess = as.numeric(lm(t~X[,-1])$coe)  #fitted linear regression
  c_guess = 1.0  #guess for shape parameter
  guess = c(beta_guess,c_guess)
  
  ####  RUBXII REG MOD MLEs    ######
  mles_vec <- vector() #To salve the estimates
  res <- optim(guess,loglik,
               method = "BFGS",
               control = list(maxit = 100, reltol = 1e-38,
                              fnscale = -1), hessian = TRUE)
  mles_vec <- res$par 
  
  # MLESs
  hat_qVEC = lfunc(c(mles_vec[1:(dim(X)[2])]),X)
  hat_c = mles_vec[ncol(X)+1]
  
  
  ### STANDARD ERRORS 
  H = res$hessian
  J = -H
  J_inv = solve(J)
  
  ### SUMMARY FOR RUBXII REG MODEL 
  mles_se <- sqrt(diag(J_inv))   
  t_stat_vec <- mles_vec/mles_se    
  K = ncol(X)+1  
  df <- K 
  df_mod <- n-K   
  p_value <- 2*(1-pt(abs(t_stat_vec),df_mod))   
  
  summ = matrix(NA,nrow = K, ncol = 4)
  if(regression == "RUBXII"){
    row.names(summ) <- c(names(data.frame(X)),"c")
  }
  if(regression == "Kw"){
    row.names(summ) <- c(names(data.frame(X)),"d_p")
  }
  if(regression == "UW"){
    row.names(summ) <- c(names(data.frame(X)),"beta")
  }
  colnames(summ) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  summ[,1] <- mles_vec[1:K]
  summ[,2] <- mles_se[1:K]
  summ[,3] <- t_stat_vec[1:K]
  summ[,4] <- p_value[1:K]
  
  #######################################################
  ## DIAGNOSTIC ANALYSIS - GOODNESS-OF-FIT MEASURES  ####
  #######################################################
  #QUANTILE RESIDUALS
  quant_res = qnorm(cdf(z,hat_qVEC,hat_c))
  ad_stat = as.numeric(shapiro.test(quant_res)[2])  
  
  ##Generalized R2 (R2_G)
  ell_hat_reg = loglik(mles_vec) 
  
  mles_vec_noreg <- vector()
  res_noreg <- optim(guess_noreg,ell_null,
                     method = "BFGS",
                     control = list(maxit = 100, reltol = 1e-38,
                                    fnscale = -1), hessian = TRUE)
  mles_vec_noreg <- res_noreg$par   
  hat_c_noreg <- res_noreg$par[1]
  hat_q_noreg <- res_noreg$par[2]   
  
  ell_hat_NOreg = ell_null(mles_vec_noreg)   ##maximized_log-lik
  
  R2_G = 1-exp(-(2/n)*(ell_hat_reg-ell_hat_NOreg))
  
  ###   GAIC: useful to compare nested models ######
  GAIC <- function(phi){-2*loglik(mles_vec)+(phi*df)}
  #where df is the number of estimate parameters
  good_fit_meas_mat <- matrix(NA, nrow = 1,ncol = 3)
  colnames(good_fit_meas_mat) <- c("AIC","R-square","p-val(SW)")
  rownames(good_fit_meas_mat) <- c(regression)
  good_fit_meas_mat[1,] <- c(GAIC(2), R2_G,
                             as.numeric(shapiro.test(quant_res)[2]))
  
  good_fit_meas_mat  #goodness-fit measures results
  
  
  #Useful quantities
  coeffbeta_mat = matrix(NA,nrow = 1,ncol = K-1)
  beta_tag <- paste0(rep("x_",(K-2)),2:(K-1))   
  colnames(coeffbeta_mat) <- c("(Intercept)",beta_tag)
  coeffbeta_mat[1,] <- mles_vec[1:(K-1)]
  
  coeff_c <- matrix(NA,nrow = 1,ncol = 1)
  colnames(coeff_c) <- c("(Intercept)")
  coeff_c[1,1] <- mles_vec[K]
  
  max_ell <- loglik(mles_vec)
  
  beta_hat <-  mles_vec[1:(K-1)]
  eta_hat <- X%*%beta_hat
  
  #Results as a list
  results = list(summary = summ, 
                 q.coeff = coeffbeta_mat,
                 c.coeff = coeff_c,
                 q.fv = hat_qVEC, 
                 c.fv = hat_c, 
                 residuals = quant_res,
                 diagnosticMEASURES = good_fit_meas_mat,
                 max_ell = max_ell,
                 eta_hat = eta_hat)
  
  return(results)
}

####################################################################################
#                                RESIDUALS PLOT FUNCTION  
####################################################################################
plots_quantres <- function(quant_residuals, fitted_values, sample_size){
  PPP <- par(mfrow = c(1,1))
  #PLOT 1: hist and estimated density
  # hist(quant_residuals, xlab="Quantile residuals", ylab="Density", 
  #      xlim = c(-3,3), ylim = c(0,0.6),breaks = 9,
  #      col="grey", freq = F, main = "Histogram")
  # 
  # x = seq(from=-4, to=4, length=200)
  # curve(dnorm(x, mean=mean(quant_residuals),
  #             sd=sqrt(var(quant_residuals))), add=TRUE, col=1, lwd = 1.5)
  
  # ##PLOT 2: Worm plot (wp)
  # wp(resid = quant_residuals, col = 1, ylim.all =2)
  # title(main = "Worm plot")
  # ##PLOT 3: qrs versus index
  # plot(1:n,quant_residuals, main = "Residual plot",      #We expected no trend.
  #      xlab = "Index", pch=20,
  #      ylab = "Quantile residuals", ylim = c(-4,4))  
  # abline(h=-3, lty=3)
  # abline(h=-2, lty=2)
  # abline(h=0)
  # abline(h=2,lty=2)
  # abline(h=3, lty=3)
  
  # PLOT 2: NORMAL QQ-PLOT
  qqnorm(quant_residuals,pch=1,frame=T, main = "QQ-plot",
         make.plot = T, lwd=1)
  qqline(quant_residuals,col="red")
  return(par(PPP))
}
