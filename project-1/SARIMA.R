library(tidyverse)
library(readxl)
library(ggplot2)
library(magrittr)
library(forecast)
library(stargazer)
library(tseries)


setwd("C:\\Users\\Gunna\\Documents\\Skole\\År 4\\Tidsrekker\\Proj 1")

# Make dataframe
df <- read.csv("csv.csv", header = T,sep = ";")
colnames(df) = c("Dato", "Kumulativt.antall", "Nye.tilfeller")
df$Dato = as.Date(df$Dato, "%d.%m.%y")
# Set the period of the data
df$Nye.tilfeller = ts(df$Nye.tilfeller, frequency = 7, start = (8+4/7))

# SARIMA = auto.arima(df$Nye.tilfeller
#                     , D = 1
#                     , lambda = "auto" #take best BoxCox
#                     , stepwise = F # look through all models
#                     , approximation = F # no approximations
#                     , trace = T # give information
#                     , max.order = 10
#                     , allowdrift = F
#                     )

SARIMA = Arima(df$Nye.tilfeller, order = c(3,1,2), seasonal = c(0,1,2), lambda = "auto")
#Look at model info
SARIMA

pred = forecast(SARIMA, h = 14) # h = days ahead to predict


autoplot(SARIMA)
par(mfrow = c(1,1), cex = 1.5)
plot(pred, xlab = "Weeks", ylab = "New cases", )

