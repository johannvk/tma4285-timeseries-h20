# Laste libraries.
library(forecast)
library(stats)
library(ggplot2)
library(ggTimeSeries)
library(tseries)
library(rugarch)
library(scales)
library(TSA)
library(car)
library(itsmr)

# Laste data.
# covid.data = read.csv("C:\\Users\\isskj\\Documents\\NTNU\\Tidsrekker\\Project 2\\covid_ex8.csv", header=TRUE, sep=";")
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
Box.test(residuals.new, lag=10, type=c("Ljung-Box"), fitdf = sarima.order)

# Ljung-Box p? squared residuals fra sarima. (aka. McLeod-Li)
Box.test(residuals.new^2, lag=10, type="Ljung-Box", fitdf=sarima.order)



# ----------------------- GARCH ------------------------------------

# Lager GARCH-modell med normal noise.
garch.model = ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(0,1)),
                         mean.model = list(armaOrder=c(0,0)),
                         distribution.model = "norm"
                         )
# Tilpasser garch-modell til residuals fra sarima.
garch.fit = ugarchfit(garch.model, residuals.new)
garch.fit


# Lager GARCH-modell med t-noise.

garch.model.t = ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(1,1)),
                           mean.model = list(armaOrder=c(0,0)),
                           distribution.model = "std"
                           )
# Tilpasser garch med t-noise til residuals fra sarima.
garch.fit.t = ugarchfit(garch.model.t, residuals.new)

# Henter garch residuals.
residuals_garch_norm = residuals(garch.fit)
residuals_garch_t = residuals(garch.fit.t)

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

