# Plotting packages:
library(ggplot2)
library(ggTimeSeries)
library(gridExtra)

# Time series packages:
library(forecast)
library(sarima)
library(tseries)
library(TSA)

# GARCH-packages:
library(rugarch)
library(scales)
library(car)
library(itsmr)


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
  n_step = max(10*period, ARMA.Obj$nobs)
  
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

autoplot(forecast:::forecast.Arima(SARIMA_8, h=14, bootstrap = T))

set.seed(321)
SARIMA_8.boot = Boot.SARIMA(SARIMA_8, num_trials=2000, 
                            max.sigma2=5*SARIMA_8$sigma2)
print("Done Bootstrapping!")

# Results:
length(SARIMA_8.boot$sigma2.estimates)
SARIMA_8.boot$param.std.dev
SARIMA_8.boot$sigma2.std.dev

########################################################################
############ Comparing old and new SARIMA-model:          ##############
########################################################################

# Import old data set:
cv19_old = read.csv2("covid_19_old.csv", header=TRUE, sep=";")
colnames(cv19_old) <- c("Dato", "Kum.Ant", "Nye.Tilf")
cv19_old$Dato = as.Date(cv19_old$Dato, format="%d.%m.%y")

cv19_old$Nye.Tilf = ts(cv19_old$Nye.Tilf, frequency=7, start=(8 + 4/7))

# Model from exercise 3:
SARIMA_3 = Arima(cv19_old$Nye.Tilf, order = c(3,1,2), 
                 seasonal=list(order=c(0,1,2),period=7), 
                 lambda = "auto")

# Predict two weeks ahead from the new cv19-data, with the old model:
SARIMA_3_refit = Arima(cv19$Nye.Tilf, model=SARIMA_3)
pred_3 = forecast(SARIMA_3_refit, h=14)

# Predict two weeks ahead with the new model estimated on the new data:
pred_8 = forecast(SARIMA_8, h=14)

# Plot the two predictions: Kan garantert bli penere plot!
par(mfrow = c(1,2), cex = 1.005)
ylim = c(0, 1.05*max(pred_8$model$x, 
                     #                     pred_3$upper, pred_8$upper,  # Gav veldig stor y-lim. 
                     pred_8$mean, pred_3$mean))
plot(pred_3, 
     ylim=ylim, 
     include=50, 
     main="14-Day Pred. Exerc. 3:
SARIMA(3,1,2)(0,1,2)[7]",
     xlab="Weeks", ylab="Daily New CV19 Cases")
plot(pred_8, 
     ylim=ylim, 
     include=50, 
     main="14-Day Pred. Exerc. 8:
SARIMA(3,1,2)(2,0,0)[7]",
     xlab="Weeks", ylab="Daily New CV19 Cases")

pred_3_plot = autoplot(pred_3, 
                       ylim=c(0, 400), 
                       include=50, 
                       main="14-Day Pred. Exerc. 3:
SARIMA(3,1,2)(0,1,2)[7]",
                       xlab="Weeks", ylab="Daily New CV19 Cases")
pred_8_plot = autoplot(pred_8, 
                       ylim=c(0, 400), 
                       include=50, 
                       main="14-Day Pred. Exerc. 8:
SARIMA(3,1,2)(2,0,0)[7]",
                       xlab="Weeks", ylab="Daily New CV19 Cases")
# Not as nice plots with autoplot/ggplot-functionality. 
# Confidence intervals are cut off and not shown properly.
# grid.arrange(pred_3_plot, pred_8_plot, ncol=2)

#################################################
#   Simulate five realizations of each model    #
#################################################

# Example use of simulate.Arima:
# sim_ts = forecast:::simulate.Arima(fixed_SARIMA, seed=223, future=F, 
#                                    nsim=length(fixed_SARIMA$x), 
#                                    lambda=fixed_SARIMA$lambda)

# Containers:
model3.ts = list(); model8.ts = list()
model3.pred = list(); model8.pred = list()
pred_h = 14  # Predict two weeks ahead.

seed0 = 1234; nobs = 8*7 # Simulate 8 Weeks of data first.
num.realizations = 6; 
for (i in 1:num.realizations) {
  # future=T: Simulated points conditional on the previously observed data.
  # Ensures positive outcomes for the simulated time series.
  
  model3.ts[[i]] = ts(forecast:::simulate.Arima(SARIMA_3_refit, seed=seed0 + i,
                                                future=T, 
                                                nsim=nobs, 
                                                lambda=SARIMA_3_refit$lambda),
                      frequency=7, start=(0))
  model3.pred[[i]] = forecast(Arima(model3.ts[[i]], model=SARIMA_3), h=pred_h)
  
  model8.ts[[i]] = ts(forecast:::simulate.Arima(SARIMA_8, seed=seed0 - i,
                                                future=T,  
                                                nsim=nobs, 
                                                lambda=SARIMA_8$lambda),
                      frequency=7, start=(0))
  model8.pred[[i]] = forecast(Arima(model8.ts[[i]], model=SARIMA_8), h=pred_h)
}

# Autoplot-centered plotting of simulated paths and their two-week predictions:
pred8.plots = list()
pred3.plots = list()
for (i in 1:(num.realizations)) {
  pred8.plots[[i]] = autoplot(model8.pred[[i]],
                              xlab="Weeks", 
                              ylab="New CV19-cases",
                              main=paste("Simulation ", as.character(i), sep="")
  )
  pred3.plots[[i]] = autoplot(model3.pred[[i]],
                              xlab="Weeks", 
                              ylab="New CV19-cases",
                              main=paste("Simulation ", as.character(i), sep="")
  )
}
do.call("grid.arrange", c(pred8.plots, ncol=3)) 
do.call("grid.arrange", c(pred3.plots, ncol=3)) 


###############################################################################
########################        GARCH-fitting      ############################
###############################################################################

# Laste data.
covid.data = read.csv("covid_19_new.csv", header=TRUE, sep=";")


# Endrer kolonnenavn.
colnames(covid.data)<- c("dato", "kumulativt.antall", "nye.tilf")
covid.data$dato <- as.Date(covid.data$dato, format="%d.%m.%y")



# Bruker oppdaterte tall fra FHI om antall meldte smittede 12 og 13 oktober. 
covid.data[length(covid.data$dato)-1, 3] = 162
covid.data[length(covid.data$dato), 3] = 107



# Lager timeseries-objekt. 
covid.nye.ts = as.ts(covid.data$nye.tilf)

# Definerer frekvens, og starttidspunkt i tidsrekken v?r.
covid.data$nye.tilf = ts(covid.data$nye.tilf, frequency=7,start = (8+4/7)) 


# Lager en SARIMA-modell ved hjelp av auto.arima. Bruker 7-dagers perioden spesifisert som frekvens i dataene.
sarima_ny = forecast::auto.arima(covid.data$nye.tilf,
                                 lambda="auto",
                                 stepwise=F,
                                 approximation=F,
                                 trace=T,
                                 max.order=10)

# Finne residuals fra sarima.
residuals.new = resid(sarima_ny)

# Henter ut (p+q) fra sarima-objekt.
sarima.order = sarima_ny$arma[1]+sarima_ny$arma[2]

# Ljung-Box test p? residuals fra sarima.
Box.test(residuals.new, lag=14, type=c("Ljung-Box"), fitdf = sarima.order)

# Ljung-Box p? squared residuals fra sarima. (aka. McLeod-Li)
Box.test(residuals.new^2, lag=14, type="Ljung-Box", fitdf=sarima.order)



# ----------------------- GARCH ------------------------------------

# Lager GARCH-modell med normal noise.
garch.model = ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(0,1)),
                         mean.model = list(armaOrder=c(0,0)),
                         distribution.model = "norm",
                         fixed.pars=list(mu=0)
)
# Tilpasser garch-modell til residuals fra sarima.
garch.fit = ugarchfit(garch.model, residuals.new)

garch_normal_order_sum = 0 + 1

# Lager GARCH-modell med t-noise.
garch.model.t = ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(1,1)),
                           mean.model = list(armaOrder=c(0,0)),
                           distribution.model = "std",
                           fixed.pars=list(mu=0)
)
# Tilpasser garch med t-noise til residuals fra sarima.
garch.fit.t = ugarchfit(garch.model.t, residuals.new)

garch_t_order_sum = 1 + 1



# ---------------  Residuals normal-----------------------
# Henter residuals fra Normal GARCH.
garch_fitted_norm = sigma(garch.fit)
residuals_garch_norm = as.vector(residuals.new) - as.vector(garch_fitted_norm)

# Tester
Box.test(residuals_garch_norm, type="Ljung",lag=14, fitdf = garch_normal_order_sum)
Box.test(residuals_garch_norm^2, type="Ljung",lag=14, fitdf=garch_normal_order_sum)

# Lage normal qq plot
qqnorm(residuals_garch_norm, ylab="GARCH residuals", xlab="Theoretical Normal quantiles", main="")
qqline(residuals_garch_norm)


# ------------------ Residuals t -----------------------------
# Henter residuals fra Student t GARCH
garch_fitted_t = sigma(garch.fit.t)
residuals_garch_t = as.vector(residuals.new) - as.vector(garch_fitted_t)

# Test på t residuals
Box.test(residuals_garch_t, type="Ljung",lag=14, fitdf=garch_t_order_sum)
Box.test(residuals_garch_t^2, type="Ljung",lag=14, fitdf=garch_t_order_sum)

# Finne frihetsgrader i t-fordeling i garch-modell
shape_est = garch.fit.t@fit$coef["shape"]

# Lage qq plot vs teoretisk t-fordeling
qqPlot(residuals_garch_t, distribution="t", df=shape_est, envelope=FALSE, 
       col.lines="black", lwd=1, ylab="GARCH residuals", xlab="Theoretical t quantiles",
       grid=FALSE)


# Tester fra residualene direkte om noen av P-verdiene er mindre
# enn alpha = 0.05-nivået for Ljung-Box-test ved lags = 1:lag.max.

alpha = 0.05
garch.ljung.box = LjungBoxTest(residuals_garch_norm, 
                               k=0, lag.max=30, SquaredQ=T)
any_p_less_alpha = any(garch.ljung.box[ , "pvalue"] < alpha)
any_p_less_alpha  # FALSE

garch.t.ljung.box = LjungBoxTest(residuals_garch_t, 
                                 k=0, lag.max=30, SquaredQ=T)
any_p_less_alpha = any(garch.t.ljung.box[ , "pvalue"] < alpha)
any_p_less_alpha  # FALSE
