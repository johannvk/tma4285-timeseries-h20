library(forecast)

# Load data:
cv19 = read.csv2("covid19.csv", header=TRUE, sep=";")
colnames(cv19) <- c("Dato", "Kum.Ant", "Nye.Tilf")
cv19$Dato = as.Date(cv19$Dato, "%d.%m.%y")

# BoxCox-transform first:
best.lambda = BoxCox.lambda(cv19$Nye.Tilf)
Nye.Tilf.BXCX = BoxCox(cv19$Nye.Tilf, best.lambda)

# Perform 7-lag then 1-lag differencing:
Nye.Tilf.BXCXdiff71 = diff(diff(Nye.Tilf.BXCX, lag=7), lag=1)
Nye.Tilf.BXCXdiff71_MZ = Nye.Tilf.BXCXdiff71 - mean(Nye.Tilf.BXCXdiff71)

# Plotting results of data transformation:
# plot(Nye.Tilf.BXCXdiff71_MZ)
# ggtsdisplay(Nye.Tilf.BXCXdiff71_MZ)


"----------- ARMA Model Estimation --------------"
# By choosing max.order = 10, we get an Arma(3, 5)-model.
# There we have some coefficients with absolute value greater than 1.
# Is this ok? Should we limit ourselves to models with coefficients
# with absolute value less than 1?

arma_Nye.Tilf = auto.arima(Nye.Tilf.BXCXdiff71_MZ, d=0, ic="aicc", max.order=5, 
                           stepwise=FALSE, approximation=FALSE, 
                           stationary=TRUE, lambda=NULL)
arma_Nye.Tilf
lines(arma_Nye.Tilf$fitted)

#------- Estimate of uncertainty in the phi- and theta-parameters found ------
#--------------------- by the Arima-function: --------------------------------

phi0 = arma_Nye.Tilf$model$phi
theta0 = arma_Nye.Tilf$model$theta
sigma2 = arma_Nye.Tilf$sigma2

# Number of simulation trials:
m = 2000
# Steps in each trial:
n_step = 500
# Order:
Arima_Ord = c(1, 0, 1)
phi_vec = c()
theta_vec = c()
sigma2_vec = c()

for(i in 1:m){
  sim_ts = arima.sim(list(order = c(1,0,1), ar=phi0 , ma=theta0), 
                     sd=sqrt(sigma2), n=n_step)
  arma_mod = Arima(sim_ts, order=Arima_Ord, include.mean=FALSE, lambda=NULL)
  
  phi_vec = c(phi_vec, arma_mod$model$phi)
  theta_vec = c(theta_vec, arma_mod$model$theta)
  sigma2_vec = c(sigma2_vec, arma_mod$sigma2)
}
mean(phi_vec)
sqrt(var(phi_vec))

mean(theta_vec)
sqrt(var(theta_vec))

mean(sigma2_vec)
sqrt(var(sigma2_vec))
