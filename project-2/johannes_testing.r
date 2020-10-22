library(ggplot2)
library(ggTimeSeries)

# Time series packages:
library(forecast)
library(sarima)
library(tseries)

simulate.Sarima = function(Sarima.Obj, n=NULL, before.data=NULL){
  # Function returning a function, which when called returns an n-point 
  # simulation of a SARIMA-process, as given in Sarima.Obj.
  # n: The Number of data points the simulated time series should contain. 
  # If n is not specified, the default is as many data points as was used
  # in fitting the Sarima.Obj.
  # before.data: Data specifying the start of the simulated time series. 
  # If before.data is specified, the simulated time-series will start with 
  # those data, but will only return "fresh" simulated data-points.

  model_specification = Sarima.Obj$arma
  p = model_specification[[1]]
  q = model_specification[[2]] 
  P = model_specification[[3]]
  Q = model_specification[[4]]
  s = model_specification[[5]]
  d = model_specification[[6]]
  D = model_specification[[7]]
  
  obj_sigma2 = sqrt(Sarima.Obj$sigma2)
  Arima_Mod = c(p, d, q)
  Seasonal_Mod = c(P, D, Q)
  
  # ARMA coefficients:
  if (p > 0){
    phi_indexes = paste(rep("ar", p), as.character(1:p), sep="")
    phis = Sarima.Obj$coef[phi_indexes]
  }else{
    phis = as(0.0, "BJFilter")
  }
  if (q > 0){
    theta_indexes = paste(rep("ma", q), as.character(1:q), sep="")
    thetas = Sarima.Obj$coef[theta_indexes]
  } else { 
    thetas = as(0.0, "BJFilter")
  }
  
  # Seasonal parts:
  if (P > 0){
    Phi_indexes = paste(rep("sar", P), as.character(1:P), sep="")
    Phis = Sarima.Obj$coef[Phi_indexes]
  } else {    
    Phis = as(0.0, "BJFilter") 
  }
  if (Q > 0){
    Theta_indexes = paste(rep("sma", Q), as.character(1:Q), sep="")
    Thetas = Sarima.Obj$coef[Theta_indexes]
  } else {    
    Thetas = as(0.0, "BJFilter") 
  }
  
  # Steps in each trial:
  if (is.null(n)){
    n_step = Sarima.Obj$nobs
  }
  else{
    n_step = n
  }
  sarima_model = list(ar=phis, ma=thetas, sar=Phis, sma=Thetas,
                      sigma2=obj_sigma2, iorder=d, siorder=D, nseasons=s)
  
  # Specify whether or not to include "begin.data":
  if (is.null(before.data)){
    sarima_simulation = prepareSimSarima(model=sarima_model, n=n_step)
  }
  else{
    sarima_simulation = prepareSimSarima(model=sarima_model, n=n_step,
                                         x=list(before=before.data))
  }
  
  return(sarima_simulation)
}

cv19 = read.csv2("covid_19_new.csv", header=TRUE, sep=";")
colnames(cv19) <- c("Dato", "Kum.Ant", "Nye.Tilf")
cv19$Dato = as.Date(cv19$Dato, format="%d.%m.%y")

# Adjust the last two data-points for new cases:
cv19$Nye.Tilf[length(cv19$Nye.Tilf)-1] = 162
cv19$Nye.Tilf[length(cv19$Nye.Tilf)] = 107

cv19$Nye.Tilf = ts(cv19$Nye.Tilf, frequency=7, start=(8 + 4/7))
# ggtsdisplay(diff(diff((cv19$Nye.Tilf), 7), 1))
# print(cv19$Nye.Tilf, calendar=TRUE)

plot(cv19$Nye.Tilf)
# head(cv19)

# No Box-Cox transformation:
auto_SARIMA = auto.arima(cv19$Nye.Tilf, 
                                # d=1, 
                                ic='aicc',
                                max.order=12, 
                                # stationary=T, 
                                lambda=NULL, 
                                stepwise=F, approximation=F)
# Result: ARIMA(1,1,5)(1,0,1)[7], AICc=2222.6
no_box_cox = Arima(cv19$Nye.Tilf, order=c(1,1,5), 
                   seasonal=list(order=c(1,0,1), period=7), 
                   lambda = 'auto')

# Allow for Box-Cox transformation of data:
auto_SARIMA_lambda = auto.arima(cv19$Nye.Tilf, 
                         # d=1, 
                         ic='aicc',
                         max.order=12, 
                         # stationary=T, 
                         lambda='auto', 
                         stepwise=F, approximation=F)
# lambda= 0.1849332
# Result: ARIMA(3,1,2)(2,0,0)[7]: AIC=754.83, AICc=755.47, BIC=782.51
# With ic='bic':
# ARIMA(0,1,1)(2,0,0)[7]: AIC=762.36, AICc=762.54, BIC=776.2

fixed_SARIMA = Arima(cv19$Nye.Tilf, order=c(3,1,2), 
               seasonal=list(order=c(2,0,0),period=7), 
               lambda = 'auto')
str(fixed_SARIMA)
#  Testing forecast-package simulation: 
# Need three ':' to use an un-exported function directly.
sim_ts = forecast:::simulate.Arima(fixed_SARIMA, seed=223, future=F, 
                                   nsim=length(fixed_SARIMA$x), 
                                   lambda=fixed_SARIMA$lambda)
ndiffs(sim_ts)
plot.ts(sim_ts)

# Test simulating a new time series based on the fixed_SARIMA model:
sim_func = simulate.Sarima(fixed_SARIMA)
ts.plot(sim_func())

# autoplot(fixed_SARIMA)
Acf(fixed_SARIMA$residuals^2)

plot(forecast(auto_SARIMA_lambda, h=20), ylim=c(0, 400))
