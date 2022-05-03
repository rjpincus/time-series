#########################################
#   Time Series Analysis HW2            #
#   02Sep2021                           #
#   Ozone data ESM                      #
#########################################


#Load libraries
library(tidyverse)
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
library(lubridate)

#Set directory
setwd('C:/Users/Richard Pincus/Documents/Classes - MSA/AA502/Time Series/data/HW2')

#Read in data
energy = read.csv('UK.csv')
head(energy)
str(energy)

#Check if any months are missing in a year
energy %>% mutate(year = format(my(energy$Date), format='%Y')) %>%
                     mutate(mon = format(my(energy$Date), format='%M')) %>% 
                     group_by(year) %>% dplyr::summarise(n())

#Looks like 2021 is missing some months but I think this is where the data stops - April 2021


#Plot data first


#Create Time Series object
energy.ts = ts(energy$Hydro_energy, start=2006, frequency=12)

#Plot
autoplot(energy.ts)

#Split into training, validation, and test datas
train = energy[1:(nrow(energy)-17),]
val   = energy[(nrow(energy)-16):(nrow(energy)-5),]
test  = energy[(nrow(energy)-4):nrow(energy),]

#Create Time Series object for train
train.ts = ts(train$Hydro_energy, start=2006, frequency = 12)

#Plot train TS
autoplot(train.ts)





#Decompose data to see seasonality and trend using STL
energy.stl = stl(train.ts, s.window = 7)

#Plot STL
plot(energy.stl)

#Create plot of Data with Trend and Seasonal decomps overlaid
# summary(energy.stl$time.series)
# autoplot(train.ts) + 
#   geom_line(aes(y=energy.stl$time.series[,2],color="maroon")) +
#   geom_line(aes(y=energy.stl$time.series[,1] + mean(energy.stl$time.series[,2]), color='red')) + 
#   xlab('Year') + 
#   ylab('Net Generation of Electricity (UNIT)') + 
#   scale_colour_manual(name = 'Lines', 
#                       values =c('black'='black','red'='red', 'maroon'='maroon'), 
#                       labels = c('Data','Seasonality','Trend'))

#Create plot of Data with Trend decomp overlaid
summary(energy.stl$time.series)
autoplot(train.ts) + 
  geom_line(aes(y=energy.stl$time.series[,2],color="red")) +
  xlab('Year') + 
  ylab('Net Generation of Electricity (UNIT)') + 
  scale_colour_manual(name = 'Lines', 
                      values =c('black'='black','red'='red'), 
                      labels = c('Data','Trend'))

#Plot of Actual and Seasonally adjusted 
seas_adj=train.ts-energy.stl$time.series[,1]

autoplot(train.ts) +
  geom_line(aes(y=seas_adj,color="red")) +
  xlab('Year') + 
  ylab('Net Generation of Electricity (UNIT)') + 
  scale_colour_manual(name = 'Lines', 
                      values =c('black'='black','red'='red'), 
                      labels = c('Data','Seasonally \nAdjusted Data'))











#Create Simple Exponential Smoothing Model
energy.simp.esm = ses(train, initial='simple',h=17)
summary(energy.simp.esm)
plot(energy.simp.esm)
#MAPE=169.842

#Create Additive Holt ESM - Additive
energy.add <- hw(train.ts, seasonal = "additive", h=12)
summary(energy.add)
plot(energy.add)
#MAPE=16.39188

#Create Multiplicative Holt ESM - Multiplicative
energy.mult <- hw(train.ts, seasonal = "multiplicative", h=12)
summary(energy.mult)
plot(energy.mult)
#MAPE=16.91997

#Create Damped ESM just to check
energy.damp.esm <- holt(train.ts, initial = "optimal", h = 17, damped = TRUE)
summary(energy.damp.esm)
plot(energy.damp.esm)
#MAPE=20.72544




#I will check both Holt ESM Milti and Add on the Validation test set now

#Create validation time series
val.ts = ts(val$Hydro_energy, start=c(2019,12), frequency = 12)




#Calcualte Error against the validation TS for Additive model
error.add=val.ts-energy.add$mean

#Calculate MAPE for Additive model
MAPE=mean(abs(error.add)/abs(val.ts))*100
MAPE



#Calcualte Error against the validation TS for Multiplicative model
error.mult=val.ts-energy.mult$mean

#Calculate MAPE for Additive model
MAPE=mean(abs(error.mult)/abs(val.ts))*1--
MAPE


#Combine Validation and Training to retrain model
train2 = energy[1:(nrow(energy)-5),]

#Create Time Series object for train
train2.ts = ts(train2$Hydro_energy, start=2006, frequency = 12)





#Create Additive Holt ESM - Additive
energy.add2 <- hw(train2.ts, seasonal = "additive", h=5)
summary(energy.add2)
plot(energy.add2)
#MAPE=16.50171




#Create test Time Series
test.ts = ts(test$Hydro_energy, start=c(2020,12), frequency = 12)

#Calculate error against test data for final model
energy.add2=test.ts-energy.add2$mean

#Calculate MAPE for Additive model
MAPE=mean(abs(energy.add2)/abs(test.ts))*100
MAPE
#MAPE=12.5351

