#########################################
#   Time Series II HW2                  #
#   16Oct2021                           #
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
setwd('C:/Users/Richard Pincus/Documents/Classes - MSA/AA502/Time Series/TS II/HW/HW2')
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
autoplot(energy.ts)

#Split into training and validation
train = head(energy.ts, -17)
val = subset(energy.ts, start=(length(energy.ts)-16), end=(length(energy.ts)-5))
test = tail(energy.ts, 5)





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




#Build model 1 - p=1, P=12
set.seed(18)
NN.Model1 <- nnetar(diff(train, 12), p = 1, P = 12)


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

#Not so good as the first NN


# Calculate prediction errors from forecast on Validation - NN1
NN.error1 <- val - energy.Forecast1

# Calculate prediction error statistics (MAE and MAPE)
NN.MAE1 <- mean(abs(NN.error1))
NN.MAPE1 <- mean(abs(NN.error1)/abs(val))*100

NN.MAE1
NN.MAPE1


# Calculate prediction errors from forecast on Validation - NN2
NN.error2 <- val - energy.Forecast2

# Calculate prediction error statistics (MAE and MAPE)
NN.MAE2 <- mean(abs(NN.error2))
NN.MAPE2 <- mean(abs(NN.error2)/abs(val))*100

NN.MAE2
NN.MAPE2





#Build model 2 - full season 
set.seed(18)
NN.Model3 <- nnetar(diff(train, 12), p = 12, P = )


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

#Not so good as the first NN



