# Plotting packages:
library(ggplot2)
library(ggTimeSeries)
library(gridExtra)

# Time series packages:
library(forecast)
library(sarima)
library(tseries)
library(TSA)

Boot.SARIMA = function(ARMA.Obj, num_trials=1000, 
                       lambda="auto", max.sigma2=10.0) {
  # Function for bootstrapping the standard deviations of 
  # SARIMA model parameters. 

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
  
  # Number of simulation trials:
  m = num_trials
  # Steps in each trial:
  period = s
  n_step = min(8*period, ARMA.Obj$nobs)
  
  # Initializing storage for Bootstrap-estimates:
  # Stored Column-wise for each parameter.
  parameter_estimates = matrix(0.0, nrow=m, ncol=length(ARMA.Obj$coef))
  colnames(parameter_estimates) = attributes(ARMA.Obj$coef)$names
  
  sigma2_vec = c()
  i = 1
  while(i <= m){
    # Some simulated timeseries lead to difficulties in optimizing
    # the Maximum-Likelihood function for the specified SARIMA model.
    # Those cases are skipped, and not counted towards the full number
    # of simulated time series/fitted models.
    tryCatch({
      sim_ts = ts(forecast:::simulate.Arima(ARMA.Obj,
                                            future=T, 
                                            nsim=n_step, 
                                            lambda=ARMA.Obj$lambda),
                  frequency=7, start=(0))
      arma_mod = Arima(sim_ts, order=Arima_Mod, 
                       seasonal=list(order=Seasonal_Mod, period=s), 
                       lambda=lambda, method="CSS-ML")
      parameter_estimates[i,] = arma_mod$coef
      sigma2_vec = c(sigma2_vec, arma_mod$sigma2)
      i = i + 1
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  # Remove the parameter estimates with sigma2 > max.sigma2. 
  # These models are unrealistic, and misleading the 
  # bootstrap estimates.
  parameter_estimates = parameter_estimates[sigma2_vec < max.sigma2, ]
  sigma2_vec = sigma2_vec[sigma2_vec < max.sigma2]
  
  ret_value = list()
  
  ret_value[["parameter.estimates"]] = parameter_estimates
  ret_value[["param.std.dev"]] = sqrt(apply(parameter_estimates, 
                                            MARGIN=2, var))
  ret_value[["sigma2.estimates"]] = sigma2_vec
  ret_value[["sigma2.std.dev"]] = sqrt(var(sigma2_vec)) 
  return (ret_value)
}

# Import data set and correct the last two datums:
cv19 = read.csv2("covid_19_new.csv", header=TRUE, sep=";")
colnames(cv19) <- c("Dato", "Kum.Ant", "Nye.Tilf")
cv19$Dato = as.Date(cv19$Dato, format="%d.%m.%y")

# Adjust the last two data-points for new cases:
cv19$Nye.Tilf[length(cv19$Nye.Tilf)-1] = 162
cv19$Nye.Tilf[length(cv19$Nye.Tilf)] = 107

cv19$Nye.Tilf = ts(cv19$Nye.Tilf, frequency=7, start=(8 + 4/7))

# New Sarima model for the full CV19-data. Found through 
# New model for new data: 
# Found from auto.arima():
# auto_SARIMA_lambda = auto.arima(cv19$Nye.Tilf, 
#                                 ic='aicc',
#                                 max.order=12, 
#                                 lambda='auto', 
#                                 stepwise=F, approximation=F)
# Result:
# lambda= 0.1849332
# Result: ARIMA(3,1,2)(2,0,0)[7]: AIC=754.83, AICc=755.47, BIC=782.51
SARIMA_8 = Arima(cv19$Nye.Tilf, order=c(3,1,2), 
                 seasonal=list(order=c(2,0,0),period=7), 
                 lambda = 'auto')

set.seed(123)
SARIMA_8.boot = Boot.SARIMA(SARIMA_8, num_trials=2000, 
                              max.sigma2=5*SARIMA_8$sigma2)
print("Done Bootstrapping!")

# Results:
length(SARIMA_8.boot$sigma2.estimates)
SARIMA_8.boot$param.std.dev
SARIMA_8.boot$sigma2.std.dev
