# Plotting packages:
library(ggplot2)
library(ggTimeSeries)
library(gridExtra)

# Time series packages:
library(forecast)
library(sarima)
library(tseries)
library(TSA)

# Import data set and correct the last two datums:
cv19 = read.csv2("covid_19_new.csv", header=TRUE, sep=";")
colnames(cv19) <- c("Dato", "Kum.Ant", "Nye.Tilf")
cv19$Dato = as.Date(cv19$Dato, format="%d.%m.%y")

# Adjust the last two data-points for new cases:
cv19$Nye.Tilf[length(cv19$Nye.Tilf)-1] = 162
cv19$Nye.Tilf[length(cv19$Nye.Tilf)] = 107

cv19$Nye.Tilf = ts(cv19$Nye.Tilf, frequency=7, start=(8 + 4/7))

# Import old data set:
cv19_old = read.csv2("../project-1/covid19.csv", header=TRUE, sep=";")
colnames(cv19_old) <- c("Dato", "Kum.Ant", "Nye.Tilf")
cv19_old$Dato = as.Date(cv19_old$Dato, format="%d.%m.%y")

cv19_old$Nye.Tilf = ts(cv19_old$Nye.Tilf, frequency=7, start=(8 + 4/7))

# Model from exercise 3:
SARIMA_3 = Arima(cv19_old$Nye.Tilf, order = c(3,1,2), 
                 seasonal=list(order=c(0,1,2),period=7), 
                 lambda = "auto")

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

# Predict two weeks ahead from the new cv19-data, with the old model:
SARIMA_3_refit = Arima(cv19$Nye.Tilf, model=SARIMA_3)
pred_3 = forecast(SARIMA_3_refit, h=14)

# Predict two weeks ahead with the new model estimated on the new data:
pred_8 = forecast(SARIMA_8, h=14)

# Plot the two predictions: Kan garantert bli penere plot!
par(mfrow = c(1,2), cex = 1.01)
ylim = c(0, 1.05*max(pred_8$model$x, 
#                     pred_3$upper, pred_8$upper,  # Gav veldig stor y-lim. 
                     pred_8$mean, pred_3$mean))

plot(pred_3, ylim=ylim, include=50, 
     main="Two-Week Pred. Exerc. 3 Model:
     SARIMA(3,1,2)(0,1,2)[7]",
     xlab="Weeks", ylab="Daily New CV19 Cases")
plot(pred_8, ylim=ylim, include=50, 
     main="Two-Week Pred. by Exerc. 8 Model:
     SARIMA(3,1,2)(2,0,0)[7]",
     xlab="Weeks", ylab="Daily New CV19 Cases")


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
