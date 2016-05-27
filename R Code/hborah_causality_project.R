# Load the libraries
library(vars)
library(urca)
library(pcalg)
library(stats)
library(lmtest)
library(forecast)
library(tseries)


# Read the input data
csv_file_path <- "/Users/Himangshu/Google Drive/NC_State_Spring_2016/Samatova!/Projects/Causal/Input Data/data.csv"
tetrad_csv_path <- "/Users/Himangshu/Google Drive/NC_State_Spring_2016/Samatova!/Projects/Causal/Input Data/tetrad_data.csv"

# Read the data as normal dataframe and a time series separately.
data_causal <- read.csv(csv_file_path, sep = ",", header = TRUE)
ts_causal <- ts(data_causal, frequency = 52)

# Build a VAR model
# Treat all the models as endogeneous as we are not sure about any external influence factors.
rprice_slice <- data_causal[, 2]
mprice_slice <- data_causal[, 3]
qty_slice <- data_causal[, 1]

# First glance of the data.
plot.ts(rprice_slice)
plot.ts(mprice_slice)
plot.ts(qty_slice)


# Check the Lag value according to Schwartz Information Criterion
VARselect(data_causal, lag.max = 10, type="const")
# The output of the above code
# $selection
# AIC(n)  HQ(n)  SC(n) FPE(n)
#   7      2      1      7

# SC criteria suggests a lag of 1, i.e p = 1
# build the VAR model for the 3 variables. Type is "const" as suggested.
var_model <- VAR(data_causal, p=1, type="const")

# Extract the residuals from the VAR model
qty_res <- var_model$varresult$Move$residuals
rprice_res <- var_model$varresult$RPRICE$residuals
mprice_res <- var_model$varresult$MPRICE$residuals

# Visualize the residuals.
plot(qty_res)
plot(mprice_res)
plot(rprice_res)


# Find lag using nDiffs.
qty_lag <- ndiffs(qty_res)
mprice_lag <- ndiffs(mprice_res)
rprice_lag <- ndiffs(rprice_res)

# Check for stationarity using the Augmented Dickey-Fuller test
qty_df <- ur.df(qty_res, lags = qty_lag, type="none")
mprice_df <- ur.df(mprice_res, lags = mprice_lag, type="none")
rprice_df <- ur.df(rprice_res, lags = rprice_lag, type="none")

summary(qty_df)
summary(mprice_df)
summary(rprice_df)

# All 3 significantly small p value indicates that we can reject the NULL Hypothesis
# that the series are non-staionary. hence they are indeed stationary.

# Check whether the variables follow a Gaussian distribution
ks.test(qty_res, "pnorm", mean(qty_res), sd(qty_res))
ks.test(rprice_res, "pnorm", mean=mean(rprice_res), sd=sd(rprice_res))
ks.test(mprice_res, "pnorm", mean=mean(mprice_res), sd=sd(mprice_res))

# Low p-value suggests that we reject the hypothesis that the true ditribution of the errors
# data is not same as the standard normal distribution "pnorm", the residuals
# are not normally distributed!
#
# Write the residuals to a csv file to build causal graphs using Tetrad software
res_com <- cbind(qty_res, rprice_res, mprice_res)
write.csv(res_com, file = tetrad_csv_path, row.names = FALSE)
#
# R Code for PC and LiNGAM
# PC algorithm
suffStat=list(C=cor(res_com), n=nrow(res_com))
pc_fit <- pc(suffStat, indepTest=gaussCItest, alpha=0.1, labels=colnames(res_com), skel.method="original")
plot(pc_fit, main="PC Output")


# LiNGAM algorithm
lingam_fit <- LINGAM(res_com, verbose = TRUE)
show(lingam_fit)

