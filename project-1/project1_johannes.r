library(forecast)

cv19 = read.csv2("covid19.csv", header=TRUE, sep=";")

colnames(cv19) <- c("Dato", "Kum.Ant", "Nye.Tilf")
cv19$Dato = as.Date(cv19$Dato, "%d.%m.%y")

# Looking at what transformations of our data could result
# in what appears to be a stationary process:
hist(cv19$Nye.Tilf)
cv19$Nye.Tilf.Box_Cox = BoxCox(cv19$Nye.Tilf, 0.33)
Nye.Tilf_bc_lm = lm()
ggtsdisplay(box_cox_transformed_Nye.Tlf)

hist(box_cox_transformed_Nye.Tlf)
plot(cv19$Dato, box_cox_transformed_Nye.Tlf)

# Beginning experimentation with ARIMA-model:
arma_mod = Arima(cv19$Nye.Tilf - mean(cv19$Nye.Tilf), order=c(3, 1, 3))
str(arma_mod)
accuracy(arma_mod)
plot(arma_mod)

plot(cv19$Dato, mean_cor)
ggtsdisplay(cv19$Nye.Tilfeller - mean(cv19$Nye.Tilfeller))
str(arma_mod)

typeof(cv19)
cv19ts = as.ts(cv19, deltat=1/365)
str(cv19ts)
summary(cv19ts)
plot(cv19$Dato, log(cv19$Nye.Tilfeller))
colnames(cv19)
