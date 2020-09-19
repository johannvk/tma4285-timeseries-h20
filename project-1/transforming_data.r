library(forecast)

# Load data:
cv19 = read.csv2("covid19.csv", header=TRUE, sep=";")
colnames(cv19) <- c("Dato", "Kum.Ant", "Nye.Tilf")
cv19$Dato = as.Date(cv19$Dato, "%d.%m.%y")
# View(cv19)
# BoxCox-transform first:
best.lambda = BoxCox.lambda(cv19$Nye.Tilf)
Nye.Tilf.BXCX = BoxCox(cv19$Nye.Tilf, best.lambda)

simple_mod = auto.arima(cv19$Nye.Tilf, max.order=10,
                        stepwise=FALSE, approximation=FALSE,
                        lambda=NULL)
var(simple_mod$residuals)/simple_mod$sigma2
par(cex=1.1)
Acf(simple_mod$residuals, 
    main="ACF for Naiive ARMA(4, 1, 4) Residuals", xlab="Lag (days)",
    col.lab=1, type="correlation")

# Perform 7-lag then 1-lag differencing: Y = (1-B)(1-B^7)X
Nye.Tilf.BXCXdiff71 = diff(diff(Nye.Tilf.BXCX, lag=7), lag=1)
Nye.Tilf.BXCXdiff71_MZ = Nye.Tilf.BXCXdiff71 - mean(Nye.Tilf.BXCXdiff71)

# Plotting results of data transformation:
# plot(Nye.Tilf.BXCXdiff71_MZ)
ggtsdisplay(Nye.Tilf.BXCXdiff71_MZ)


"----------- ARMA Model Estimation --------------"
# By choosing max.order = 10, we get an Arma(3, 5)-model.
# There we have some coefficients with absolute value greater than 1.
# Is this ok? Should we limit ourselves to models with coefficients
# with absolute value less than 1?

arma_Nye.Tilf = auto.arima(Nye.Tilf.BXCXdiff71_MZ, d=0, ic="aicc", max.order=10, 
                           stepwise=FALSE, approximation=FALSE, 
                           stationary=TRUE, lambda=NULL, parallel=TRUE)
arma_Nye.Tilf

# Diagnostic: Check That the "rescaled residuals" from the 
# ARMA-model have variance 1, and resemble white noise in
# the ACF-plot:
sd(arma_Nye.Tilf$residuals)/sqrt(arma_Nye.Tilf$sigma2)
Acf(arma_Nye.Tilf$residuals)

#------- Estimate of uncertainty in the phi- and theta-parameters found ------
#--------------------- by the Arima-function: --------------------------------

# Retrieve the estimates from the main ARMA-model:
Boot.Arma = function(ARMA.Obj, num_trials=1000) {
  # Function estimating mean- and variance of parameters
  # in an ARMA(p, q)-model. Assumes zero mean.
  
  phis = ARMA.Obj$model$phi
  thetas = ARMA.Obj$model$theta
  sigma2 = ARMA.Obj$sigma2
  
  p = length(phis)
  q = length(thetas)
  Arima_Ord = c(p, 0, q)
  
  # Number of simulation trials:
  m = num_trials
  # Steps in each trial:
  n_step = max(500, 2*length(ARMA.Obj$nobs))
  
  # Initializing storage for Bootstrap-estimates:
  # Stored Column-wise for each parameter.
  phi_estimates = matrix(0.0, nrow=m, ncol=p)
  theta_estimates = matrix(0.0, nrow=m, ncol=q)
  sigma2_vec = c()
  
  for(i in 1:m){
    sim_ts = arima.sim(list(order=Arima_Ord, ar=phis , ma=thetas), 
                       sd=sqrt(sigma2), n=n_step)
    arma_mod = Arima(sim_ts, order=Arima_Ord, include.mean=T, lambda=NULL)
    
    phi_estimates[i,] = arma_mod$model$phi
    theta_estimates[i,] = arma_mod$model$theta
    sigma2_vec = c(sigma2_vec, arma_mod$sigma2)
  }
  
  # Store mean and var in columns. One row per parameter:
  phi_boot = matrix(0.0, nrow=p, ncol=2)
  colnames(phi_boot) = c("mean", "var")
  for(i in 1:p){
    phi_boot[i, 1] = mean(phi_estimates[, i])
    phi_boot[i, 2] = var(phi_estimates[, i])
  }
  
  theta_boot = matrix(0.0, nrow=q, ncol=2)
  colnames(theta_boot) = c("mean", "var")
  for(i in 1:q){
    theta_boot[i, 1] = mean(theta_estimates[, i])
    theta_boot[i, 2] = var(theta_estimates[, i])
  } 
  
  sigma2_boot = c(mean(sigma2_vec), var(sigma2_vec))
  names(sigma2_boot) = c("mean", "var")
  ret_value = list()
  ret_value[["phi"]] = phi_boot
  ret_value[["theta"]] = theta_boot
  ret_value[["sigma2"]] = sigma2_boot
  
  return (ret_value)
}

bootStraps = Boot.Arma(arma_Nye.Tilf, num_trials = 50)
bootStraps$theta
str(bootStraps)

