library(tidyverse)
library(readxl)
library(ggplot2)
library(magrittr)
library(stargazer)

# Time series Packages:
library(forecast)
library(sarima)
library(tseries)

# Make dataframe
df <- read.csv("covid19.csv", header = T, sep = ";")
colnames(df) = c("Dato", "Kumulativt.antall", "Nye.tilfeller")
df$Dato = as.Date(df$Dato, "%d.%m.%y")

# Set the period of the data
df$Nye.tilfeller = ts(df$Nye.tilfeller, frequency = 7, start = (8+4/7))

# Bootstrapping functions:
Boot.ARMA = function(ARMA.Obj, num_trials=1000) {
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

Boot.SARIMA = function(ARMA.Obj, num_trials=1000, lambda="auto") {
  # Function estimating mean- and variance of parameters
  # in an ARMA(p, q)-model. Assumes zero mean.
  model_specification = ARMA.Obj$arma
  p = model_specification[[1]]
  q = model_specification[[2]] 
  P = model_specification[[3]]
  Q = model_specification[[4]]
  s = model_specification[[5]]
  d = model_specification[[6]]
  D = model_specification[[7]]
  
  obj_sigma2 = ARMA.Obj$sigma2
  Arima_Mod = c(p, d, q)
  Seasonal_Mod = c(P, D, Q)
  
  # ARMA coefficients:
  if (p > 0){
    phi_indexes = paste(rep("ar", p), as.character(1:p), sep="")
    phis = ARMA.Obj$coef[phi_indexes]
  }else{
    phis = as(0.0, "BJFilter")
  }
  if (q > 0){
    theta_indexes = paste(rep("ma", q), as.character(1:q), sep="")
    thetas = ARMA.Obj$coef[theta_indexes]
  } else { 
    thetas = as(0.0, "BJFilter")
  }
  
  # Seasonal parts:
  if (P > 0){
    Phi_indexes = paste(rep("sar", P), as.character(1:P), sep="")
    Phis = ARMA.Obj$coef[Phi_indexes]
  } else {    
    Phis = as(0.0, "BJFilter") 
  }
  if (Q > 0){
    Theta_indexes = paste(rep("sma", Q), as.character(1:Q), sep="")
    Thetas = ARMA.Obj$coef[Theta_indexes]
  } else {    
    Thetas = as(0.0, "BJFilter") 
  }
  
  # Number of simulation trials:
  m = num_trials
  # Steps in each trial:
  n_step = max(500, ARMA.Obj$nobs)
  
  sarima_model = list(ar=phis, ma=thetas, sar=Phis, sma=Thetas,
                      sigma2=obj_sigma2, iorder=d, siorder=D, nseasons=s)
  
  sarima_simulation = prepareSimSarima(model=sarima_model, n=n_step,
                                       x=list(before=ARMA.Obj$x[1:20]))
  
  # Initializing storage for Bootstrap-estimates:
  # Stored Column-wise for each parameter.
  parameter_estimates = matrix(0.0, nrow=m, ncol=length(ARMA.Obj$coef))
  colnames(parameter_estimates) = attributes(ARMA.Obj$coef)$names
  
  sigma2_vec = c()
  
  for(i in 1:m){
    sim_ts = sarima_simulation()
    arma_mod = Arima(sim_ts, order=Arima_Mod, 
                     seasonal=list(order=Seasonal_Mod, period=s), 
                     include.mean=T, lambda=lambda, method="CSS")
    parameter_estimates[i,] = arma_mod$coef
    sigma2_vec = c(sigma2_vec, arma_mod$sigma2)
  }
  
  ret_value = list()
  
  sigma2_boot = c(mean(sigma2_vec), var(sigma2_vec))
  names(sigma2_boot) = c("mean", "var")
  
  ret_value[["parameter_estimates"]] = parameter_estimates
  ret_value[["sigma2"]] = sigma2_boot
  
  return (ret_value)
}

# Naiive ARMA model. No data transformation.
Naiive_ARMA = auto.arima(df$Nye.tilfeller, max.order=15, 
                         lambda=NULL, # No BoxCox
                         stepwise=F, # Look through all models
                         approximation=F, # no approximations
                         seasonal=F, ic="aicc")
Naiive_ARMA
naiive_coef_variance = Boot.ARMA(Naiive_ARMA, num_trials=1000)
naiive_coef_variance

# Generate SARIMA model by optimizing AICC:
SARIMA = auto.arima(df$Nye.tilfeller, 
                    D = 1, 
                    lambda = "auto", # Find best BoxCox
                    stepwise = F, # Look through all models
                    approximation = F, # No approximations
                    trace = F, # Print information
                    max.order = 10,
                    allowdrift = F)

# Bootstrapping to find a measure of uncertainty in 
# the estimate for each parameter in the SARIMA model:
SARIMA_params = Boot.SARIMA(SARIMA, num_trials=1000)
coef_est = SARIMA_params$parameter_estimates
coef_std.dev = sqrt(apply(coef_est, MARGIN=2, var))
coef_std.dev
coef_est$sigma2

# Predict into the future basef on the best SARIMA-model:
pred = forecast(SARIMA, h = 14) # h = days ahead to predict

autoplot(SARIMA)
par(mfrow = c(1,1), cex = 1.0)
plot(pred, include = 50, xlab = "Weeks", ylab = "New cases")

