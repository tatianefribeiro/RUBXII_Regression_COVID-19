###################################################################################
# PAPER: A new regression model for the COVID-19 mortality rate in the U.S. States
# SECTION: 3.The proposed regression
# GOAL: Generating the plots of repar. RUBXII distribution (with tau=0.5 fixe)
# AUTHOR: Tatiane Fontana Ribeiro
# LAST UPDATE: August 28, 2020
###################################################################################
rm(list = ls())

#Inicializations
fromx = 0
tox = 1
yliminf = 0
ylimsup = 4.5
tau = 0.5


#Probability density function
pdf_RUBXII <- function(x)
{
c*log(1/tau)/log(1+(log(1/q))^c)*
  1/(1-x)*(log(1/(1-x)))^(c-1)*
    (1+(log(1/(1-x)))^c)^(log(tau)/log(1+(log(1/q))^c)-1)
  
}

integrate(pdf_RUBXII,0,1)
###############################################################################
par(mfrow=c(1,2))
c = 0.8
q = 0.5
curve(pdf_RUBXII,from=fromx, to=tox, add = FALSE, lty=1, 
      type = "l",xlab = expression("z"),
      ylab = expression("f(z)"),ylim =c(yliminf,ylimsup), 
      col =1, lwd = 4.0)

c = 3.3
q = 0.5
curve(pdf_RUBXII,from=fromx, to=tox, add = TRUE, lty=2, 
      type = "l", ylab = expression("f(y)"),
      ylim =c(yliminf,ylimsup), col = 2, lwd = 4.0)

c = 4.6
q = 0.1
curve(pdf_RUBXII,from=fromx, to=tox, add = TRUE, lty=4, 
      type = "l", ylab = expression("f(y)"),
      ylim =c(yliminf,ylimsup), col = 3, lwd = 4.0)

c = 1.2
q = 0.8
curve(pdf_RUBXII,from=fromx, to=tox, add = TRUE, lty=5, 
      type = "l", ylab = expression("f(y)"),
      ylim =c(yliminf,ylimsup), col = 4, lwd = 4.0)

legend("topleft", c(expression(paste(c, " = 0.8, ",q,' = 0.5')),
                     c(expression(paste(c, " = 3.3, ",q,' = 0.5')),
                       c(expression(paste(c, " = 4.6, ",q,' = 0.1')),
                         c(expression(paste(c, " = 1.2, ",q,' = 0.8')))))),
       col = c(1,2,3,4),
       lty= c(1,2,4,5),
       lwd = c(4,4,4,4), bty="n", cex = 1.3)

ylimsup = 6
c = 1.1
q = 0.9
curve(pdf_RUBXII,from=fromx, to=tox, add = FALSE, lty=1, 
      type = "l", ylab = expression("f(z)"),xlab = expression("z"),
      ylim =c(yliminf,ylimsup), col = 1, lwd = 4.0)

c = 4.4
q = 0.6
curve(pdf_RUBXII,from=fromx, to=tox, add = TRUE, lty=2, 
      type = "l", ylab = expression("f(y)"),
      ylim =c(yliminf,ylimsup), col = 2, lwd = 4.0)

c = 5.0
q = 0.4
curve(pdf_RUBXII,from=fromx, to=tox, add = TRUE, lty=4, 
      type = "l", ylab = expression("f(y)"),
      ylim =c(yliminf,ylimsup), col = 3, lwd = 4.0)

c = 9.3
q = 0.2
curve(pdf_RUBXII,from=fromx, to=tox, add = TRUE, lty=5, 
      type = "l", ylab = expression("f(y)"),
      ylim =c(yliminf,ylimsup), col = 4, lwd = 4.0)

legend("topright", c(expression(paste(c, " = 1.1, ",q,' = 0.9')),
                     c(expression(paste(c, " = 4.4, ",q,' = 0.6')),
                       c(expression(paste(c, " = 5.0, ",q,' = 0.4')),
                         c(expression(paste(c, " = 0.3, ",q,' = 0.2')))))),
       col = c(1,2,3,4),
       lty= c(1,2,4,5),
       lwd = c(4,4,4,4), bty="n", cex = 1.3)



