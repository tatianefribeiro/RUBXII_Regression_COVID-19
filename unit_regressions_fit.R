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
  
  if (any(regression == c("RUBXII", "BETA", "SIMPLEX", "Kw", "UW"))){
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
      loglik <- l_UWreg
      cdf <- cdf_UW
      ell_null <- l_UW_noreg
      guess_noreg <- c(median(z),1) 
    } 
    if (regression == "BETA"){
      #Fit Beta REG MODEL
      mod_beta <- gamlss(z~X[,-1], 
                         family = BE(
                           mu.link = "logit",sigma.link = "logit"
                         ), trace= F)
      
      #Quantile residuals
      mu_fv <- as.numeric(mod_beta$mu.fv)
      sigma_fv <- as.numeric(mod_beta$sigma.fv)
      quant_res_beta <- mod_beta$residuals 
      
      #goodness-of-fit measures
      good_fit_meas_mat <- matrix(NA, nrow = 1,ncol = 3)
      colnames(good_fit_meas_mat) <- c("AIC","R-square","p-val(AD)")
      rownames(good_fit_meas_mat) <- c(regression)
      good_fit_meas_mat[1,] <- c( AIC(mod_beta,k=2),
                                  Rsq(mod_beta),
                                  as.numeric(ad.test(quant_res_beta)[2]))
      
      
      #summary BETA
      K <- ncol(X)+1
      mles_beta <- as.numeric(c(mod_beta$mu.coefficients,
                                mod_beta$sigma.coefficients))
      se_beta <- as.numeric(sqrt(diag(vcov(mod_beta))))
      t_value <- mles_beta/se_beta
      df_mod <- n-K
      p_value_beta <- 2*(1-pt(abs(t_value),df_mod))
      
      #SEs
      summ_BETA = matrix(NA,nrow = K, ncol = 4)
      row.names(summ_BETA) <- c(names(data.frame(X)),"sigma")
      colnames(summ_BETA) <- c("Estimate", "Std. Error", "t value","Pr(>|t|)")
      summ_BETA[,1] <- mles_beta
      summ_BETA[,2] <- se_beta
      summ_BETA[,3] <- t_value
      summ_BETA[,4] <- p_value_beta
      #round(summ_BETA, digits = 4)
      
      log.lik.beta = gen.likelihood(mod_beta)
      ell_beta = -log.lik.beta()
      
      results = list(summary = summ_BETA, 
                      beta.coeff = mles_beta[1:(K-1)],
                      sigma.coeff = mles_beta[K],
                      mu.fv = mu_fv, 
                      sigma.fv = sigma_fv, 
                      residuals = quant_res_beta,
                      diagnosticMEASURES = good_fit_meas_mat,
                      max_ell = ell_beta,
                      eta_hat = X%*%mles_beta[1:(K-1)]
                      )

      return(results)
    } 
    if (regression == "SIMPLEX"){
      #Fit SIMPLEX REG MODEL
      mod_SIMPLEX <- gamlss(z~X[,-1], 
                            family = SIMPLEX(
                              mu.link = "logit"), trace= F)
      
      mu_fv <- as.numeric(mod_SIMPLEX$mu.fv)
      sigma2_fv <- as.numeric(mod_SIMPLEX$sigma.fv)
      quant_res_SIMPLEX <- mod_SIMPLEX$residuals 
      
      #goodness-of-fit measures
      good_fit_meas_mat <- matrix(NA, nrow = 1,ncol = 3)
      colnames(good_fit_meas_mat) <- c("AIC","R-square","p-val(AD)")
      rownames(good_fit_meas_mat) <- c(regression)
      good_fit_meas_mat[1,] <- c( AIC(mod_SIMPLEX,k=2),
                                  Rsq(mod_SIMPLEX),
                                  as.numeric(ad.test(quant_res_SIMPLEX)[2]))
      
      #summary SIMPLEX
      K <- ncol(X)+1
      mles_SIMPLEX <- as.numeric(c(mod_SIMPLEX$mu.coefficients,
                                   mod_SIMPLEX$sigma.coefficients))
      se_SIMPLEX <- as.numeric(sqrt(diag(vcov(mod_SIMPLEX))))
      t_value <- mles_SIMPLEX/se_SIMPLEX
      df_mod <- n-K
      p_value_SIMPLEX <- 2*(1-pt(abs(t_value),df_mod))
      
      #SEs
      summ_SIMPLEX = matrix(NA,nrow = K, ncol = 4)
      row.names(summ_SIMPLEX) <- c(names(data.frame(X)),"sigma2")
      colnames(summ_SIMPLEX) <- c("Estimate", "Std. Error", "t value","Pr(>|t|)")
      summ_SIMPLEX[,1] <- mles_SIMPLEX
      summ_SIMPLEX[,2] <- se_SIMPLEX
      summ_SIMPLEX[,3] <- t_value
      summ_SIMPLEX[,4] <- p_value_SIMPLEX
      round(summ_SIMPLEX, digits = 4)
      
      log.lik.SIMPLEX = gen.likelihood(mod_SIMPLEX)
      ell_SIMPLEX = -log.lik.SIMPLEX()
      
      results = list(summary = summ_SIMPLEX, 
                     beta.coeff = mles_SIMPLEX[1:(K-1)],
                     sigma.coeff = mles_SIMPLEX[K],
                     mu.fv = mu_fv, 
                     sigma2.fv = sigma2_fv, 
                     residuals = quant_res_SIMPLEX,
                     diagnosticMEASURES = good_fit_meas_mat,
                     max_ell = ell_SIMPLEX,
                     eta_hat = X%*%mles_SIMPLEX[1:(K-1)]
      )
      
      return(results)
    } 
  }
  else {
    stop(paste(regression, "Regression not available, available regressions are \"RUBXII\", ",
               "\"beta\",","\"simplex\",","\"Kw\", and \"UW\""))
  }
  
  #############################################
  #It fits a RUBXII, Kw or UW regression model
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
  ad_stat = as.numeric(ad.test(quant_res)[2])  
  
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
  colnames(good_fit_meas_mat) <- c("AIC","R-square","p-val(AD)")
  rownames(good_fit_meas_mat) <- c(regression)
  good_fit_meas_mat[1,] <- c(GAIC(2), R2_G,
                             as.numeric(ad.test(quant_res)[2]))
  
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
#                                 LOOCV FUNCTION
####################################################################################
LOOCV.unit_reg <- function(z_comp,X_comp, regression = "Kw"){
  
  if (any(regression == c("RUBXII", "BETA", "SIMPLEX", "Kw", "UW"))){
   
    data_comp = data.frame(cbind(z_comp,X_comp[,-1]))
    
    n_comp.CV <- length(data_comp$z_comp)
    X_comp.CV <- model.matrix(z_comp~.,data = data_comp) 
    z_comp.CV <- data_comp$z_comp
    
    n_train.CV <- n_comp.CV-1

    # LOOCV LOOP
    mae_vec.CV = vector()
    for(c in 1:n_comp.CV){
      
      # Salve the cth observation
      z_test.CV <- z_comp.CV[c]
      X_test.CV <- X_comp.CV[c,]
      
      # Extract the cth observation
      z_train.CV <- z_comp.CV[-c]
      X_train.CV <- X_comp.CV[-c,]
      
      if (any(regression == c("RUBXII", "Kw", "UW"))){
        #Fit the regression model
        mod.CV = UnitReg.fit(z_train.CV,X_train.CV,n_train.CV,regression = "Kw")
        
        ## Do prediction for the validation observation
        z_predict.CV = lfunc(c(mod.CV$q.coeff),X_test.CV)
      }

      if (regression == "BETA"){
        #Fit the model
        mod_beta_train <- gamlss(z_train.CV~X_train.CV[,-1], 
                                 family = BE(
                                   mu.link = "logit",sigma.link = "logit"
                                 ), trace= F)
        
        ## prediction
        beta_vec_hat.CV = as.numeric(mod_beta_train$mu.coefficients)
        z_predict.CV = lfunc(c(beta_vec_hat.CV),X_test.CV)
      }
      
      if (regression == "SIMPLEX"){
        #Fit the model
        mod_simplex_train <- gamlss(z_train.CV~X_train.CV[,-1], 
                                    family = BE(
                                      mu.link = "logit",sigma.link = "logit"
                                    ), trace= F)
        
        ## prediction
        simplex_vec_hat.CV = as.numeric(mod_simplex_train$mu.coefficients)
        z_predict.CV = lfunc(c(simplex_vec_hat.CV),X_test.CV)
      }
      
      # adequacy measures
      mae.CV = mean(abs(z_test.CV-z_predict.CV))
      
      # salve the adeq. measures
      mae_vec.CV[c] = mae.CV
    }
    
    #Results matrix
    res_CV = matrix(NA, nrow = 1,ncol = 1)
    colnames(res_CV) <- c("MAE")
    res_CV[1,1] <- c(mean(mae_vec.CV))
    
    return(res_CV)
  }
    else {
      stop(paste(regression, "Regression not available, available regressions are \"RUBXII\", ",
                 "\"beta\",","\"simplex\",","\"Kw\", and \"UW\""))
    }
}

####################################################################################
#                                RESIDUALS PLOT FUNCTION  
####################################################################################
plots_quantres <- function(quant_residuals, fitted_values, sample_size){
  PPP <- par(mfrow = c(2,2))
  #PLOT 1: hist and estimated density
  hist(quant_residuals, xlab="Quantile residuals", ylab="Density", 
       xlim = c(-3,3), ylim = c(0,0.6),breaks = 9,
       col="grey", freq = F, main = "Histogram")
  
  x = seq(from=-4, to=4, length=200)
  curve(dnorm(x, mean=mean(quant_residuals),
              sd=sqrt(var(quant_residuals))), add=TRUE, col=1, lwd = 1.5)
  ##PLOT 2: Worm plot (wp)
  wp(resid = quant_residuals, col = 1, ylim.all =2)
  title(main = "Worm plot")
  ##PLOT 3: qrs versus index
  plot(1:n,quant_residuals, main = "Residual plot",      #We expected no trend.
       xlab = "Index", pch=20,
       ylab = "Quantile residuals")  
  abline(h=-2, lty=2)
  abline(h=0)
  abline(h=2,lty=2)
  ##PLOT 4: NORMAL QQ-PLOT
  qqnorm(quant_residuals,pch=1,frame=T, main = "QQ-plot",
         make.plot = T, lwd=1)
  return(par(PPP))
}
