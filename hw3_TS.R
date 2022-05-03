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

#Split into train
train = energy[1:(nrow(energy)-17),]
val   = energy[(nrow(energy)-16):(nrow(energy)-5),]
test  = energy[(nrow(energy)-4):nrow(energy),]



#Set alpha level
a = 1-pchisq(log(nrow(train)),1)
a

#Create a time series object
energy.ts = ts(train$Hydro_energy, start=2006, frequency=12)
energy.ts

#Plot
autoplot(energy.ts)

#STL
decomp_stl = stl(energy.ts, s.window=7)
plot(decomp_stl)

#Seasonally Adjust Data
seas_adj=energy.ts-decomp_stl$time.series[,1]
autoplot(seas_adj)

#Check stationarity of the season adjusted time series using Dickey Fuller test Type 2 - up to lag 2
a
aTSA::adf.test(seas_adj)
# aTSA::adf.test(energy.ts)

#Up to lag 2 in Type 2 there are significant p-vals so thus we have stationary

#Check if up to lag 2 less then alpha
aTSA::adf.test(energy.ts)$type2[1:3,3]<a

#Plot Autocorrelation and Partial autocorrelation
# acf1=Acf(energy.ts, lag=10)$acf
acf1=Acf(seas_adj, lag=10)$acf
acf1
#The ACF looks less then ideal, we would love to take a seasonal difference to remove some 
#dependency structure still present in this data - regular TS

#The ACF for the seasonal adjusted data looks much better although there seems to still be a little
#bit of structure in there. I think it would be due to the changing magnitude of the seasons. 

# Pacf1=Pacf(energy.ts, lag=10)$acf
Pacf1=Pacf(seas_adj, lag=10)$acf
Pacf1

#The PACF looks less then ideal, we would love to take a seasonal difference to remove some 
#dependency structure still present in this data - regular TS

#The PACF looks like regular exponential decay for this seasonally adjusted TS.




# #Create differenced column in seasonal adjust dataframe
# # energy$HE_diff = energy$Hydro_energy - lag(energy$Hydro_energy)
# seas_adj
# stats::lag(seas_adj)
# seas_adj_lag1 = seas_adj-stats::lag(seas_adj)

#Create differenced dataframe for energy data
energy.ts.diff1 = energy.ts - stats::lag(energy.ts)

#Get time series object
# energy.ts1 = ts(energy$HE_diff, start=2006, frequency = 12)

# #Check stationarity of the lagged season adj. time series using Dickey Fuller test - up to lag 2
# a
# aTSA::adf.test(seas_adj_lag1)
#Check stationarity of the differenced time series using Dickey Fuller test - up to lag 2
a
aTSA::adf.test(energy.ts.diff1)

# #Check if up to lag 2 less then alpha
# aTSA::adf.test(seas_adj_lag1)$type1[1:3,3]<a

#Check if up to lag 2 less then alpha
aTSA::adf.test(energy.ts.diff1)$type2[1:3,3]<a

#Up to lag 2 in Type 1 there are all non-significant p-vals nonetheless for the sake of the assignment
#we shall proceed. 



# #Plot Autocorrelation and Partial autocorrelation 
# acf1=Acf(seas_adj_lag1, lag=10)$acf
# acf1
# #ACF seems to fall within the error lines after lag 1 so it is an acceptable structure
# 
# Pacf1=Pacf(seas_adj_lag1, lag=10)$acf
# Pacf1 
# #The Partial ACF mostly falls within the errors lines after lag 1 so this is also an acceptable structure


#Plot Autocorrelation and Partial autocorrelation
acf2=Acf(energy.ts.diff1, lag=10)$acf
acf2
#The ACF looks less then ideal, we would love to take a seasonal difference to remove some 
#dependency structure still present in this data

Pacf2=Pacf(energy.ts.diff1, lag=10)$acf
Pacf2
#The PACF looks less then ideal, we would love to take a seasonal difference to remove some 
#dependency structure still present in this data


#Diffrerence TS looks disgusting, never do this again



###Take a second difference and check again
energy.ts.diff2 = energy.ts.diff1 - stats::lag(energy.ts.diff1)

#Check stationarity of the differenced time series using Dickey Fuller test - up to lag 2
a
aTSA::adf.test(energy.ts.diff2)


#Check if up to lag 2 less then alpha
aTSA::adf.test(energy.ts.diff2)$type1[1:3,3]<a

#Up to lag 2 there are all significant tests thus we have stationarity


#Plot Autocorrelation and Partial autocorrelation
acf3=Acf(energy.ts.diff2, lag=10)$acf
acf3
#The ACF looks fire after lag 1. Exponential decrease after 1 lag.

Pacf3=Pacf(energy.ts.diff2, lag=10)$acf
Pacf3
#The PACF is porbably fine. There are some lags that have partial autocorrelations above the error lines
#but nothing too concerning. 



#Use Ljung-Box test to test for White Noise from the first time series no difference
index1=seq(1,length(Pacf1))

White.LB <- rep(NA, 10)
for(i in 1:10){
  White.LB[i] <- Box.test(energy.ts, lag=i, type="Ljung-Box", fitdf = 0)$p.value
}

white.dat=data.frame(cbind(White.LB,index1))
colnames(white.dat)=c("pvalues","Lag")

#Plot Ljung-Box Test p-vals
ggplot(white.dat,aes(x=factor(Lag),y=pvalues))+geom_col()+labs(title="Ljung-Box test p-values",x="Lags",y="p-values")+coord_cartesian(ylim = c(0, 0.025))

#All p-values are significant indicating our series is now white noise. Modeling complete





### for seaosnally adjusted data

#Use Ljung-Box test to test for White Noise from the first time series no difference
index1=seq(1,length(Pacf1))

White.LB <- rep(NA, 10)
for(i in 1:10){
  White.LB[i] <- Box.test(seas_adj, lag=i, type="Ljung-Box", fitdf = 0)$p.value
}

white.dat=data.frame(cbind(White.LB,index1))
colnames(white.dat)=c("pvalues","Lag")

#Plot Ljung-Box Test p-vals
ggplot(white.dat,aes(x=factor(Lag),y=pvalues))+geom_col()+labs(title="Ljung-Box test p-values",x="Lags",y="p-values")+coord_cartesian(ylim = c(0, 0.025))

#All p-values are significant indicating our series still has a lot of strucutre to model out before we 
#reach white noise

