###################################################################################
# PAPER: A new regression model for the COVID-19 mortality rate in the U.S. States
# SUBSECTION: 3.2. Simulation study
# GOAL: Conducting a simulation study for the RUBXII regression model
# AUTHOR: Tatiane Fontana Ribeiro
# LAST UPDATE: August 28, 2020
###################################################################################
rm(list = ls())
##########################     FIXED QUANTITIES     ################################
set.seed(2020)
R = 10000 #number of repplications
tau = 0.5  #known
#####################    *******   SCENARIO  ******  ################################
##*** ATTENTION ****
#CHANGE here + last line, file.Rdata    #######!!!!!!!!!!!!!!!!!!!!
beta_true = c(-0.5, -2.8)  #regression parameters
c = 3.2 #shape parameter
####################################################################################
par_true = c(beta_true,c)
beta1 = beta_true[1]
beta2 = beta_true[2]

#############  AUXILIARY FUNCTIONS AND SOME IMPORTANT COMENTS   ####################
r_RUBXII <- function(n)   #It generates occurences of Y_i ~ RUBXII (q_i, c)
{
  u = runif(n)
  y_RUBXII = 1-exp(-((1-u)^(log(1+(log(1/(1-q_i)))^c)/log(1-tau))-1)^(1/c)) #RUBXII qf
  return(y_RUBXII)
}

lfunc <- function(beta,x2){  #logit link function --> lfunc == logistic function
  #because logit funct is the inverse from the logistic or sigmoidal funct
  b1 <- beta[1]
  b2 <- beta[2]
  qi_logit_link <- 1/(1+exp(-(b1+b2*x2))) #=logit inverse div by exp(eta_i)
  return(qi_logit_link)    #Note that here eta_i = b1+b2*x_{i2} because x2 is
  #a n-dim vector with n obs from covariate X2.
}



l_RUBXIIreg <- function(par){  #Log-likelihood function
  beta1 <- par[1]
  beta2 <-par[2]
  c <- par[3]
  q_i <- lfunc(c(beta1,beta2),x2)  #It call the link function
  
  ell =  sum(
    log(
      c*log(1/(1-tau))*(log(1/(1-y)))^(c-1)*
        (1+(log(1/(1-y)))^c)^(log(1-tau)/log(1+(log(1/(1-q_i)))^c)-1)/  
        ((1-y)*log(1+(log(1/(1-q_i)))^c))    
    )
  )
  return(ell)
}
################################################################################
k = 1    #It index results final matrix              
vn = c(30,90,160,300)   #sample sizes vector
res_mat = matrix(NA,length(vn),10) #It save the final results

###########################    ARRAY BOXPLOT   #################################
arr_boxplot = array(NA, c(R,3,length(vn))) # array (obj, c(i, j, k))
# i=nrow; j=ncol and k=dimension
colnames(arr_boxplot) <- c("hat_beta1","hat_beta2","hat_c")
################################################################################

for (n in vn) {##Loop for sample size
  time_start = Sys.time()  # initial time of the hesimulation  
  x2 = rnorm(n)# runif(n,-3,3)  #constant for each n
  q_i = lfunc(c(beta1,beta2),x2) #q vector constant for each 
  #logit link function
  yr <- r_RUBXII(n) #It use previous q_i (i.e., beta1 and beta2)
  #and c defined in the start.
  
  for (j in 1:length(yr)) {
    if(yr[j] == 1) {yr[j] <- 0.9999}
  }
  
  
  #response variable for MQO as initial guess
  y =  log(yr/(1-yr))
  
  ##Initial guess
  beta_guess = lm(y~x2)  #fitted linear regression (simple in this case)
  beta1_guess <- as.numeric(beta_guess$coefficients[1]) #guess for beta1
  beta2_guess <- as.numeric(beta_guess$coefficients[2]) #guess for beta2
  c_guess = 1.0  #guess for shape parameter
  guess = c(beta1_guess,beta2_guess,c_guess)
  
  #guess = c(1,1,1)
  
  ##MC simulation
  bug = 0  #to count the bug numbers from MC simulation
  i = 1  #index MLEs matrix
  estim = matrix(NA,R,3)  #matrix for to salve the MLEs
  while (i <= R) {  ##Loop MC simulation
    #note the "<=". It is correct because "i" starts in 1
    #whether "i" to start in zero, you must to use "<"
    y <- r_RUBXII(n)  #Y ~ RUBXII (q_i, c)
    #Here, you must to have attention to the definition
    #from response vector from the log-lik;
    #mandatorily you must to define as y because you called offset
    #in the log-lik function (l_RUBXIIreg function).
    for (j in 1:length(y)) {
      if(y[j] == 1) {y[j] <- 0.9999}
    }
    res <- try(optim(guess,l_RUBXIIreg,# grr_analytic,-->> to use after checking der. in matrix not.
                     method = "BFGS",
                     control = list(maxit = 100, reltol = 1e-38,
                                    fnscale = -1)),
               silent = T)
    if(class(res)=="try-error" || res$conv != 0)
    {
      bug <- bug+1  #counting the bugs
    }else{
      estim[i,] <- res$par     
      arr_boxplot[i,,k] <- res$par #ARRAY to salve all estimates for each n
      #Note that in the place of j, it is void because 
      #we want all columns.
      i <- i+1  #it is update the repplication
    }
  }##end loop MC
  
  #(1) Useful quantities
  mean_MLEs <- apply(estim, 2, mean)   #mean of the MLEs (2 --> by columns)
  #from help: if estim = matrix, then 1 indicates rows, 2 ind. col.
  RB <- (mean_MLEs-par_true)/par_true*100  #percentage relative bias
  RMSE <- sqrt(apply(estim[,],2,var)+(mean_MLEs-par_true)^2)  # root mean square error
  #Notice that RB and RMSE are 3-dimensional vectors
  
  #Results matrix
  colnames(res_mat) <- c("beta_1","beta_2","c","n","RB_b1","RB_b2", "RB_c",
                         "RMSE_b1","RMSE_b2","RMSE_c")
  res_mat[k,] <- c(par_true,n,RB,RMSE)
  print(k)  #to help identify a possible infinite loop
  k = k+1  
  time_fin = Sys.time()  #final time of the simulation
} ##end loop n

##Useful quantities for file.tex
time = time_fin-time_start  #Comput the simulation time
stargazer::stargazer(res_mat, digits = 4)  #It build table for latex
time   #It prints the total time from the simualation of each scenario.
res_mat  #It prints the results final matrix

##It saves for boxplots and other graphics which you wish.
save.image("/home/tatiane/Insync/tfr1@de.ufpe.br/Google Drive/master_thesis/3-appCOVID19/new_models_dataset2/reg_RUBXII_scen4.RData")

