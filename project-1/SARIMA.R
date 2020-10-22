library(tidyverse)
library(readxl)
library(ggplot2)
library(magrittr)
library(forecast)
library(stargazer)
library(tseries)

library(sarima)

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

# setwd("C:\\Users\\Gunna\\Documents\\Skole\\?r 4\\Tidsrekker\\Proj 1")

# Make dataframe
df <- read.csv("covid19.csv", header = T,sep = ";")
colnames(df) = c("Dato", "Kumulativt.antall", "Nye.tilfeller")
df$Dato = as.Date(df$Dato, "%d.%m.%y")
# Set the period of the data
df$Nye.tilfeller = ts(df$Nye.tilfeller, frequency = 7, start = (8+4/7))

Naiive_ARMA = Arima(df$Nye.tilfeller, order=c(4,1,4), lambda=NULL)
naiive_coef = Boot.SARIMA(Naiive_ARMA, num_trials=1000, lambda=NULL)
naiive_est = naiive_coef$parameter_estimates
naiive_coef.std.dev = sqrt(apply(naiive_est, MARGIN=2, var))
naiive_coef.std.dev
Naiive_ARMA


# SARIMA = auto.arima(df$Nye.tilfeller
#                     , D = 1
#                     , lambda = "auto" #take best BoxCox
#                     , stepwise = F # look through all models
#                     , approximation = F # no approximations
#                     , trace = T # give information
#                     , max.order = 10
#                     , allowdrift = F
#                     )

SARIMA = Arima(df$Nye.tilfeller, order = c(3,1,2), seasonal = c(0,1,2), 
               lambda = "auto")
#Look at model info
SARIMA

SARIMA_params = Boot.SARIMA(SARIMA, num_trials=1000)
coef_est = SARIMA_params$parameter_estimates
coef_std.dev = sqrt(apply(coef_est, MARGIN=2, var))
coef_std.dev

# Prediction-plot:
pred = forecast(SARIMA, h = 14) # h = days ahead to predict

autoplot(SARIMA)
par(mfrow = c(1,1), cex = 1.05)
plot(pred, include = 50, xlab = "Weeks", ylab = "New cases", )


#######################################################
        # Add new data to Sarima-model:
#######################################################

df2 <- read.csv("../project-2/covid_19_new.csv", header = T,sep = ";")
colnames(df2) = c("Dato", "Kumulativt.antall", "Nye.tilfeller")
df2$Dato = as.Date(df2$Dato, "%d.%m.%y")
View(df2)

# Adjust the last two data-points for new cases:
df2$Nye.tilfeller[length(df2$Nye.tilfeller)-1] = 162
df2$Nye.tilfeller[length(df2$Nye.tilfeller)] = 107

# Set the period of the data along with its starting point:
df2$Nye.tilfeller = ts(df2$Nye.tilfeller, frequency = 7, start = (8+4/7))

SARIMA_refit = Arima(df2$Nye.tilfeller, model=SARIMA)
pred2 = forecast(SARIMA_refit, h = 14) # h = days ahead to predict
pred3 = forecast(SARIMA_refit, h = 14)
par(mfrow = c(1,1), cex = 1.05)

# How to predict several times, with different random-number inputs?
plot(pred2, include = 100, xlab = "Weeks", ylab = "New cases", ylim=c(0, 250))
plot(pred3, include = 100, xlab = "Weeks", ylab = "New cases", ylim=c(0, 250))

# ggtsdisplay(df2$Nye.tilfeller)
