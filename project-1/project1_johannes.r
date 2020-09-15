library(forecast)

# Load data:
cv19 = read.csv2("covid19.csv", header=TRUE, sep=";")
colnames(cv19) <- c("Dato", "Kum.Ant", "Nye.Tilf")
cv19$Dato = as.Date(cv19$Dato, "%d.%m.%y")

# Difference the new cases before polynomial fitting:
Nye.Tilf.diff7 = diff(cv19$Nye.Tilf, differences=1, lag=7)
Nye.Tilf.diff7.lifted = Nye.Tilf.diff7 - min(Nye.Tilf.diff7) + 1.0

# Lite tillit til kun diff7BoxCox:
Nye.Tilf.diff7BXCX.lambda = BoxCox.lambda(Nye.Tilf.diff7.lifted)
Nye.Tilf.diff7BXCX = BoxCox(Nye.Tilf.diff7.lifted, Nye.Tilf.diff7BXCX.lambda)
# ggtsdisplay(Nye.Tilf.diff7)
# ggtsdisplay(Nye.Tilf.diff7BXCX)

# Diff 7 - 1 - BoxCox:
Nye.Tilf.diff71 = diff(Nye.Tilf.diff7, 1)
Nye.Tilf.diff71.lifted = Nye.Tilf.diff71 - min(Nye.Tilf.diff71) + 1.0

Nye.Tilf.diff71BXCX.lambda = BoxCox.lambda(Nye.Tilf.diff71.lifted)
Nye.Tilf.diff71BXCX = BoxCox(Nye.Tilf.diff.lifted, Nye.Tilf.diff71BXCX.lambda)

# ggtsdisplay(Nye.Tilf.diff71)
ggtsdisplay(Nye.Tilf.diff71BXCX - mean(Nye.Tilf.diff71BXCX))
# n_diff_required = ndiffs(Nye.Tilf.diff71BXCX - mean(Nye.Tilf.diff71BXCX))

# Dato.diff7 = cv19$Dato[8:length(cv19$Dato)]
# lm_p4_Nye.Tilf = lm(Nye.Tilf.diff7~poly(Dato.diff7, degree=4))

# Shift residuals so that min(Nye.Tilf$residuals) = 1.0 if
# there were negative values present before.
# Nye.Tilf_res_min_value = min(lm_p4_Nye.Tilf$residuals)
# Nye.Tilf.diff71.shifted_res = lm_p4_Nye.Tilf$residuals - min(c(0, Nye.Tilf_res_min_value)) + 1
Nye.Tilf_res_lambda = BoxCox.lambda(cv19$Nye.Tilf.shifted_res)
Nye.Tilf.BXCX_Res = BoxCox(lm_p4_Nye.Tilf$residuals, Nye.Tilf_res_lambda)

ggtsdisplay(Nye.Tilf.BXCX_Res)

hist(box_cox_transformed_Nye.Tlf)
plot(cv19$Dato, box_cox_transformed_Nye.Tlf)

# Beginning experimentation with ARIMA-model:
arma_mod = Arima(Nye.Tilf.BXCX_Res, order=c(3, 0, 3))
str(arma_mod)
# Extract AICC:
arma_mod$aicc

accuracy(arma_mod)
plot(arma_mod)

ggtsdisplay(cv19$Nye.Tilfeller - mean(cv19$Nye.Tilfeller))
str(arma_mod)

# typeof(cv19)
# cv19ts = as.ts(cv19, deltat=1/365)
# str(cv19ts)
# summary(cv19ts)
# plot(cv19$Dato, log(cv19$Nye.Tilfeller))
# colnames(cv19)

# Perform differences to get a stationary timeseries, Box-Cox afterwards:
Nye.Tilf.diff7 = diff(cv19$Nye.Tilf, differences=1, lag=7)
Nye.Tilf.diff71 = diff(Nye.Tilf.diff7, differences=1, lag=1)

# Lift the data so that min(X) = 1.0 for BoxCox-transform:
Nye.Tilf.diff71.lifted = Nye.Tilf.diff71 - min(Nye.Tilf.diff71) + 1.0
best.lambda = BoxCox.lambda(Nye.Tilf.diff71.lifted)
Nye.Tilf.diff71BXCX = BoxCox(Nye.Tilf.diff.lifted, best.lambda)
mean_Nye.Tilf.diff71BXCX = mean(Nye.Tilf.diff71BXCX)
Nye.Tilf.diff71BXCX_MZ = Nye.Tilf.diff71BXCX - mean_Nye.Tilf.diff71BXCX