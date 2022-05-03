#########################################
#   Time Series II Project              #
#   26Oct2021                           #
#   Splash Energy                       #
#########################################


#Load libraries
library(tseries)
library(forecast)
library(haven)
library(fma)
library(expsmooth)
library(lmtest)
library(zoo)
library(seasonal)
library(ggplot2)
library(seasonalview)
library(aTSA)
library(imputeTS)
library(prophet)
library(lubridate)
library(tidyverse)




#Read in data
setwd('C:/Users/Richard Pincus/Documents/Classes - MSA/AA502/Time Series/TS II/HW/project')
energy = read.csv('UK.csv')
head(energy)
tail(energy)




#Check for any missing months in a year
energy %>% mutate(year = format(my(energy$Date), format='%Y')) %>%
  mutate(mon = format(my(energy$Date), format='%M')) %>% 
  group_by(year) %>% dplyr::summarise(n())

#No missing months





#Create time series object
energy.ts = ts(energy$Hydro_energy, start = 2006, frequency = 12)

#Plot data
autoplot(energy.ts) # + geom_text(label = index(energy.ts), nudge_y = 1, check_overlap = T) 


#Split into training and validation
train = head(energy.ts, -17)
val = subset(energy.ts, start=(length(energy.ts)-16), end=(length(energy.ts)-5))
test = tail(energy.ts, 5)


autoplot(train)






###
###
### Facebook's Prophet Model
###
###

#Build a prophet model to forecast on validation
#Since it is just curve fitting no stationarity is needed to model


#Holidays to consider in model - Energy Act Passed on 2013-06-13 to eliminated energy
#production from coal 30% by 2020, 50% by 2025 and 80% by 2050
holidays <- data.frame(
  holiday = 'EnergyAct2013',
  ds = as.Date(c('2013-06-13', '2014-06-13', '2015-06-13',
                 '2016-06-13')),
  lower_window = 0,
  upper_window = 0
)
#Structure data
prophet.data <- data.frame(ds = seq(as.Date('2006-01-01'), as.Date('2019-11-01'), by = 'm'), y = train)

#Build model
Prof <- prophet(holidays = holidays)
Prof <- add_country_holidays(Prof, "UK")
Prof <- add_seasonality(Prof, name='monthly', period=30.5, fourier.order=6)
Prof <- fit.prophet(Prof, prophet.data)


#Forecast
forecast.data <- make_future_dataframe(Prof, periods = 12, freq = 'month')
predict(Prof, forecast.data)$yhat


#Plot forecasts
plot(Prof, predict(Prof, forecast.data))
prophet.forecast = tail(predict(Prof, forecast.data)$yhat, 12)
prophet.forecast = ts(prophet.forecast, start=c(2019,12), frequency=12)

plot(train, main = "UK Energy Prophet  Model Forecasts", xlab = "Date", ylab = "Hydorelectric Energy (GWh)", xlim = c(2006, 2021))
lines(prophet.forecast, col = "blue")
abline(v = 2019.917, col = "red", lty = "dashed")


# Calculate prediction errors from forecast
Prophet.error <- val - tail(predict(Prof, forecast.data)$yhat, 12)


# Calculate prediction error statistics (MAE and MAPE)
Prophet.MAE <- mean(abs(Prophet.error))
Prophet.MAPE <- mean(abs(Prophet.error)/abs(val))*100

Prophet.MAE #87.60173
Prophet.MAPE #14.19768










###
###
### Neural Network Model
###
###

#Stationarity needed for this model -> take 1 seasonal difference from last homework

#take seasonal difference
train1 = train %>% diff(lag = 12)


#plot data
autoplot(train1)


#Test if seasonal difference needs to be taken again 
train1 %>% nsdiffs() #no


#Decompose again using STL and plot
sdiff.train.stl = stl(train1, s.window=7)
autoplot(sdiff.train.stl)

#No trend -> test for stationarity now


#ADF for stationarity up to 2 lags
adf.test(train1)
#all significant lags -> this means ts is stationary




#Check ACF
acf0 = Acf(train1, lag=30)$acf
#q = 2 but I cannot use MA terms in Neural Network. I will try one with the full
#seaosn and one model with just the AR terms


#Check PACF
pacf0 = Pacf(train1, lag=30)$acf
#p = 1, P = 12




#Build model 1 - p=1, P=1
set.seed(18)
NN.Model1 <- nnetar(diff(train, 12), p=1)


#Forecast differences
NN.Forecast1 <- forecast::forecast(NN.Model1, h = 12)
plot(NN.Forecast1)


#Get forecasts in original data, not differenced data
energy.Forecast1 <- rep(NA, 12)

for(i in 1:12){
  energy.Forecast1[i] <- train[length(train) - 12 + i] + forecast::forecast(NN.Model1, h = 24)$mean[i]
}

energy.Forecast1 <- ts(energy.Forecast1, start = c(2019, 12), frequency = 12)

plot(train, main = "UK Energy Neural Network Model Forecasts", xlab = "Date", ylab = "Hydorelectric Energy (GWh)", xlim = c(2006, 2021))
lines(energy.Forecast1, col = "blue")
abline(v = 2019.917, col = "red", lty = "dashed")

NN.Model1 #sigma^2 estimated as 0.007944

#I think that means this model is pretty good



#Check ACF
acf1 = Acf(NN.Model1$residuals, lag=30)$acf
#q = 2 but I cannot use MA terms in Neural Network. I will try one with the full
#seaosn and one model with just the AR terms


#Check PACF
pacf1 = Pacf(NN.Model1$residuals, lag=30)$acf
#p = 1, P = 12


#Plot residuals
autoplot(NN.Model1$residuals)



#Check stationarity
adf.test(NN.Model1$residuals)
#All significant lags indicating stationary

#Check for white noise
checkresiduals(NN.Model1)




#Build a model on these residuals
NN.Model1.2 <- nnetar(tail(NN.Model1$residuals,143), P=1)


#Forecast differences 
NN.Forecast1.2 <- forecast::forecast(NN.Model1.2, h = 12)
plot(NN.Forecast1.2)


#Get forecasts in original data, not differenced data
energy.Forecast1.2 <- rep(NA, 12)

for(i in 1:12){
  energy.Forecast1.2[i] <- train[length(train) - 12 + i] + forecast::forecast(NN.Model1.2, h = 24)$mean[i]
}

energy.Forecast1.2 <- ts(energy.Forecast1.2, start = c(2019, 12), frequency = 12)

plot(train, main = "UK Energy Neural Network Model Forecasts", xlab = "Date", ylab = "Hydorelectric Energy (GWh)", xlim = c(2006, 2021))
lines(energy.Forecast1.2, col = "blue")
abline(v = 2019.917, col = "red", lty = "dashed")
















###
###
### Prophet again
###
###



#Holidays to consider in model - Energy Act Passed on 2013-06-13 to eliminated energy
#production from coal 30% by 2020, 50% by 2025 and 80% by 2050
holidays <- data.frame(
  holiday = 'EnergyAct2013',
  ds = as.Date(c('2013-06-13', '2014-06-13', '2015-06-13',
                 '2016-06-13')),
  lower_window = 0,
  upper_window = 0
)
#Structure data
prophet.data <- data.frame(ds = seq(as.Date('2006-01-01'), as.Date('2019-11-01'), by = 'm'), y = train)

#Build model
Prof <- prophet(holidays = holidays)
Prof <- add_country_holidays(Prof, "UK")
Prof <- add_seasonality(Prof, name='monthly', period=30.5, fourier.order=6)
Prof <- fit.prophet(Prof, prophet.data)


#Forecast
forecast.data <- make_future_dataframe(Prof, periods = 12, freq = 'month')
predict(Prof, forecast.data)$yhat


#Plot forecasts
plot(Prof, predict(Prof, forecast.data))
prophet.forecast = tail(predict(Prof, forecast.data)$yhat, 12)
prophet.forecast = ts(prophet.forecast, start=c(2019,12), frequency=12)

plot(train, main = "UK Energy Prophet  Model Forecasts", xlab = "Date", ylab = "Hydorelectric Energy (GWh)", xlim = c(2006, 2021))
lines(prophet.forecast, col = "blue")
abline(v = 2019.917, col = "red", lty = "dashed")
















###
###
### Model Validation
###
###


# Calculate prediction errors from forecast on Validation - NN2
NN.error1.2 <- val - energy.Forecast1.2

# Calculate prediction error statistics (MAE and MAPE)
NN.MAE1.2 <- mean(abs(NN.error1.2))
NN.MAPE1.2 <- mean(abs(NN.error1.2)/abs(val))*100

NN.MAE1.2
NN.MAPE1.2 #22.20719 wawa










###########################################
## Kari's code
###########################################

#Cutoff train at 11 years before validation starts
short_train = tail(train, n = (107))
#short all
short_all = tail(energy.ts, n = length(short_train) + 18)
short_all
#########################################################
# Establish pre-existing models
#########################################################
# Seasonal Decomposition
decomp_stl <- stl(short_all, s.window = 7)
plot(decomp_stl) #Inspect plot

#From this, the shorter data window appears to have additive seasonality + 
#trend that would be fit well using a piece-wise method

##Additive hw
hw_add <- hw(short_train, seasonal = "additive", initial = "optimal", h = 12) #did additive based off of decomposition
hw_add_error <- val-hw_add$mean #Capture error from prediction
mae_hw_add <- mean(abs(hw_add_error)) #MAE
mape_hw_ad <- mean(abs(hw_add_error)/abs(val)) #MAPE

##Multiplicative hw
hw_mult <- hw(short_train, seasonal = "multiplicative", initial = "optimal", h = 12) #did additive based off of decomposition
hw_mult_error <- val-hw_mult$mean #Capture error from prediction
mae_hw_mult <- mean(abs(hw_mult_error)) #MAE
mape_hw_mult <- mean(abs(hw_mult_error)/abs(val)) #MAPE

##Seasonal ARIMA
short_train %>% nsdiffs() #Calls for 1 difference - also indicates seasonality is stochastic
short_train %>% diff(lag = 6) %>% ggtsdisplay() #Inspect plot of seasonally differenced data
seas_diff <- short_train %>% diff(lag = 6) #Create seasonally-difference time series object
aTSA::adf.test(seas_diff) #Check Dickey-Fuller, series now stationary
#Fit AR & MA terms
Acf(short_train, lag=24)$acf #Significant spike at 1, 2 - use 2 MA terms
Pacf(short_train, lag=24)$acf #No AR terms needed

#Run Model
seas_a <- short_train %>% 
  Arima(order=c(2,0,1), seasonal=c(1,1,0))
#Check residuals
checkresiduals(seas_a) #Looks good! Also Ljung-Box indicates white noise
#Forecast
seas_forc <- forecast::forecast(seas_a, h = 12)
seas_forc_error <- val-seas_forc$mean #Calculate forecast error
mae_seas_forc <- mean(abs(seas_forc_error)) #MAE
mape_seas_forc <- mean(abs(seas_forc_error)/abs(val)) #MAPE

##auto arima
auto <- auto.arima(short_train, method="ML", seasonal = TRUE) #Run auto ARIMA algorithm to check results
summary(auto) #Inspect summary
auto_forc <- forecast::forecast(auto, h = 12) #Forecast using auto ARIMA model
auto_forc_error <- val-auto_forc$mean #Calculate auto ARIMA forecast error
mae_auto_forc <- mean(abs(auto_forc_error)) #MAE
mape_auto_forc <- mean(abs(auto_forc_error)/abs(val)) #MAPE

### Prophet Model ###

#Create data to read into prophet model
prophet.data <- data.frame(ds = seq(as.Date('2011-01-01'), as.Date('2019-11-01'), by = 'm'), y = short_train)
holidays <- data.frame(
  holiday = 'use',
  ds = as.Date(c('')),
  ds = as.Date(c('2013-06-01', '2014-01-01', '2014-07-01', '2015-12-01', '2016-06-01')),
  lower_window = 0,
  upper_window = 1
)
Prof <- prophet(holidays = holidays) #Initialize prophet model
Prof <- add_country_holidays(Prof, "UK") #Consider UK holidays
Prof <- add_seasonality(Prof, name='monthly', period=30.5, fourier.order=6) #Model annual season using Fourier vars
#Prof <- add_seasonality(Prof, name = "weekly", period = 7, fourier.order = 3)
Prof <- fit.prophet(Prof, prophet.data) #Fit prophet model

# Forecast Prophet for validation
forecast.data <- make_future_dataframe(Prof, periods = 12, freq = 'month')
predict(Prof, forecast.data)$yhat
prophet.forecast = ts(tail(predict(Prof, forecast.data)$yhat, 12), start=c(2019, 11), frequency = 12)

# Calculate prediction error from Prophet
Prophet.error <- val - tail(predict(Prof, forecast.data)$yhat, 12)

# Calculate prediction error statistics (MAE and MAPE)
Prophet.MAE <- mean(abs(Prophet.error))
Prophet.MAPE <- mean(abs(Prophet.error)/abs(val))*100

### Neural Network ###

#Set seed to make analysis cohesive
set.seed(18)
#Calculate on one seasonal difference
NN.Model <- nnetar(diff(short_train, 12), p = 10, P = 4) #Take seasonal diff, and use p, P based on auto.arima

#Create forecast
NN.Forecast <- forecast::forecast(NN.Model, h = 12)

#Shift forecast back to original data - remove seasonal diff
Pass.Forecast <- rep(NA, 12) #Create empty forecast vector

for(i in 1:12){ #Replace with undifferenced forecast value
  Pass.Forecast[i] <- short_train[length(short_train) - 12 + i] + forecast::forecast(NN.Model, h = 12)$mean[i]
}

#Put forecast values into time series object
Pass.Forecast <- ts(Pass.Forecast, start = c(2019, 11), frequency = 12)

# Calculate prediction errors from forecast
NN.error <- val - Pass.Forecast

# Calculate prediction error statistics (MAE and MAPE)
NN.MAE <- mean(abs(NN.error))
NN.MAPE <- mean(abs(NN.error)/abs(val))*100

## sketchy neural net process
sketchy_data = data.frame(matrix(nrow = 0, ncol = 3))
colnames(sketchy_data) = c('p', 'P', 'mape')
# for (n in 1:10){
#   for (j in 1:10){
#     NN.sketchy <- nnetar(diff(short_train, 12), p = n, P = j)
#     sketchy.Forecast <- rep(NA, 12) #Create empty forecast vector
#     for(i in 1:12){ #Replace with undifferenced forecast value
#       sketchy.Forecast[i] <- short_train[length(short_train) - 12 + i] + forecast::forecast(NN.sketchy, h = 12)$mean[i]
#       sketchy.error = val - sketchy.Forecast
#       sketchy.mape = mean(abs(sketchy.error)/abs(val))*100
#     }
#     sketchy_data = rbind(sketchy_data, c(n, j, sketchy.mape))
#   }
# }

#########################################################
# Weather inputs
#########################################################
#Fit a model
rain = read.csv("C:\\Users\\kjahn\\Documents\\MSA\\Fall\\Analytics\\Time Series 2\\Final Project\\rain.csv")
xreg <- rain$rain
short_train %>%
  Arima(order=c(0,0,2), seasonal=c(1,1,1)) %>%
  residuals() %>% ggtsdisplay()
S.ARIMA <- Arima(short_train, order=c(0,0,2), seasonal=list(order=c(1,1,1), period=12), xreg=xreg)
summary(S.ARIMA)
checkresiduals(S.ARIMA)#forecast the next five observations
forecast::forecast(S.ARIMA, h = 12)autoplot(forecast::forecast(S.ARIMA, h = 12)) + autolayer(fitted(S.ARIMA), series="Fitted") +
  ylab("Hydro Energy") +
  geom_vline(xintercept = 2018.91,color="orange",linetype="dashed")# Calculate prediction errors from forecast
S.ARIMA.error <- val - forecast::forecast(S.ARIMA, h = 6)$mean# Calculate prediction error statistics (MAE and MAPE)
S.ARIMA.MAE <- mean(abs(S.ARIMA.error))
S.ARIMA.MAPE <- mean(abs(S.ARIMA.error)/abs(val))*100
S.ARIMA.MAE
S.ARIMA.MAPE





##Arima
arima4 = Arima(short_train, order=c(0,0,2), seasonal = c(1,1,1))
arima4.forecast = forecast::forecast(arima4, h = 12)

arima5 = Arima(short_train, order=c(1,0,2), seasonal = c(1,1,1))
arima5.forecast = forecast::forecast(arima5, h = 12)

arima6 = Arima(short_train, order=c(1,0,0), seasonal = c(0,1,1))
arima6.forecast = forecast::forecast(arima6, h = 12)

#MAPE
error.arima=val-arima6.forecast$mean

#Calculate MAE & MAPE
MAE.arima=mean(abs(error.arima))
MAPE.arima=mean(abs(error.arima)/abs(val))



#########################################################
# Averaged models
#########################################################
For.Avg <- (hw_add$mean + tail(predict(Prof, forecast.data)$yhat, 12))/2
Avg.error <- validation - For.Avg

Avg.MAE <- mean(abs(Avg.error))
Avg.MAPE <- mean(abs(Avg.error)/abs(validation))*100


For.Avg2 <- (hw_add$mean + Pass.Forecast)/2
Avg.error2 <- val - For.Avg2

Avg.MAE2 <- mean(abs(Avg.error2))
Avg.MAPE2 <- mean(abs(Avg.error2)/abs(val))*100





#Weighted averages for HW Additive and NN models
wc.model = lm(tail(train, 95) ~ offset(tail(hw_add$fitted, 95)) + I(NN.Model$fitted - hw_add$fitted) - 1)
summary(wc.model)
nn_coef = coef(wc.model)[1]
hw_add_coef = 1 - nn_coef


For.Avg3 <- (hw_add_coef*hw_add$mean + nn_coef*Pass.Forecast)/2
Avg.error3 <- val - For.Avg3

Avg.MAE3 <- mean(abs(Avg.error3))
Avg.MAPE3 <- mean(abs(Avg.error3)/abs(val))*100



For.Avg4 <- (hw_add$mean + prophet.forecast)/2
Avg.error4 <- val - For.Avg4

Avg.MAE4 <- mean(abs(Avg.error4))
Avg.MAPE4 <- mean(abs(Avg.error4)/abs(val))*100



autoplot(val) + autolayer(hw_add$mean) + 
  autolayer(hw_mult$mean) + 
  autolayer(prophet.forecast) + 
  autolayer(Pass.Forecast) + 
  autolayer(arima4.forecast$mean) +
  autolayer(arima5.forecast$mean) + 
  autolayer(arima6.forecast$mean) 


autoplot(val) + autolayer(For.Avg2)
autoplot(val) + autolayer(prophet.forecast) + 
  autolayer(hw_mult$mean) + 
  autolayer(hw_add$mean)








#Build model 1 - p=1, P=12
set.seed(18)
NN.Model1 <- nnetar(diff(train, 12), p = 1, P=2)


#Forecast differences
NN.Forecast1 <- forecast::forecast(NN.Model1, h = 12)
plot(NN.Forecast1)


#Get forecasts in original data, not differenced data
energy.Forecast1 <- rep(NA, 12)

for(i in 1:12){
  energy.Forecast1[i] <- train[length(train) - 12 + i] + forecast::forecast(NN.Model1, h = 24)$mean[i]
}

energy.Forecast1 <- ts(energy.Forecast1, start = c(2019, 12), frequency = 12)

plot(train, main = "UK Energy Neural Network Model Forecasts", xlab = "Date", ylab = "Hydorelectric Energy (GWh)", xlim = c(2006, 2021))
lines(energy.Forecast1, col = "blue")
abline(v = 2019.917, col = "red", lty = "dashed")

NN.Model1 #sigma^2 estimated as 0.007944

#I think that means this model is pretty good

# summary(NN.Model1)
# 
# autoplot(train1) + 
#   autolayer(NN.Model1$fitted)





#Build model 2 - full season 
set.seed(18)
NN.Model2 <- nnetar(diff(train, 12), p = 12)


#Forecast differences
NN.Forecast2 <- forecast::forecast(NN.Model2, h = 12)
plot(NN.Forecast2)


#Get forecasts in original data, not differenced data
energy.Forecast2 <- rep(NA, 12)

for(i in 1:12){
  energy.Forecast2[i] <- train[length(train) - 12 + i] + forecast::forecast(NN.Model2, h = 24)$mean[i]
}

energy.Forecast2 <- ts(energy.Forecast2, start = c(2019, 12), frequency = 12)

plot(train, main = "UK Energy Neural Network Model Forecasts", xlab = "Date", ylab = "Hydorelectric Energy (GWh)", xlim = c(2006, 2021))
lines(energy.Forecast2, col = "blue")
abline(v = 2019.917, col = "red", lty = "dashed")

NN.Model2 #sigma^2 estimated as 866.4

autoplot(val) + 
  autolayer(energy.Forecast1) + 
  autolayer(energy.Forecast2)




For.Avg4 <- (energy.Forecast1 + prophet.forecast + hw_mult$mean)/3
Avg.error4 <- val - For.Avg4

Avg.MAE4 <- mean(abs(Avg.error4))
Avg.MAPE4 <- mean(abs(Avg.error4)/abs(val))*100



autoplot(val) + 
  autolayer(For.Avg4)

            