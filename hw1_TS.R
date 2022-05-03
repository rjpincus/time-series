#########################################
#   Time Series Analysis HW1            #
#   31Aug2021                           #
#   Exploring Ozone data                #
#########################################


#Load packages
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


#Set directory
setwd('C:/Users/Richard Pincus/Documents/Classes - MSA/AA502/Time Series/data')


#Read in data
ozone = read.csv('Ozone_Raleigh2.csv')
head(ozone)
str(ozone)


#Get list of all possible years
ozone = ozone %>% mutate(Date.yr=format(as.POSIXct(Date, format = "%m/%d/%Y"), format='%Y')) %>%
  mutate(Date.mo=format(as.POSIXct(Date, format = "%m/%d/%Y"), format='%m'))
years = data.frame(table(ozone$Date.yr))
years = years$Var1


#Sum up all missing days in each year
year_sums = ozone %>% select(Date.yr) %>% group_by(Date.yr) %>% summarise(count=n())
year_sums = year_sums %>% mutate(missing_days=ifelse(Date.yr=='2020', 151-count, 
                                                     ifelse(Date.yr=='2016', 366-count, 365-count)))
year_sums %>% summarise(sum(missing_days))


#Need to consider if the data just stops in May 2020 or is it missing for the rest of 2020


#Get the mean max Ozone Concentration for Jan 2017
ozone %>% group_by(Date.yr, Date.mo) %>% 
  summarise(mean(Daily.Max.8.hour.Ozone.Concentration)) %>% 
  filter(Date.mo=='01' & Date.yr=='2017')


#Create time plot of monthly average data
ozone_mon = ozone %>% group_by(Date.yr, Date.mo) %>% 
  summarise(avg = mean(Daily.Max.8.hour.Ozone.Concentration)) 
ozone_mon

  #Get Time Series object
ozone_mon_ts <- ts(ozone_mon$avg, start = 2014, frequency =12)
  #Create time plot
autoplot(ozone_mon_ts)+labs(title="Time Series plot for Monthly Ozone Avg.", x="Date",y="Ozone Concentration")



#Decompose data to see seasonality and trend using STL
decomp_stl <- stl(ozone_mon_ts, s.window = 7)

#Plot STL
plot(decomp_stl)

