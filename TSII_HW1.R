#########################################
#   Time Series II HW1                  #
#   05Oct2021                           #
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
setwd('C:/Users/Richard Pincus/Documents/Classes - MSA/AA502/Time Series/TS II/HW/HW1')
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


#Decompose using the X13 method
decomp_x13=seas(train)
summary(decomp_x13)
#view(decomp_x13) #very nice app that allows us to inspect decomp
autoplot(decomp_x13)


#Decompose using STL
decomp_stl=stl(train, s.window=7)
autoplot(decomp_stl)





###
###
### ESM
###
###



#Use an ESM to model and forecast one season before accounting for seasonality

#Try both additive and multiplicative,one could argue for either based on decompositions

#Additive Holt Winters ESM
hw.a.energy <- hw(train, seasonal = "additive", h=12)
summary(hw.a.energy)
#MAPE = 16.36955  
#MAE = 98.49377 
#AIC = 2510.038 
#BIC = 2563.044 


#plot forecasts for add
autoplot(hw.a.energy)+
  autolayer(fitted(hw.a.energy),series="Fitted")+ylab("Energy")+ geom_vline(xintercept = 2020.4166667,color="orange",linetype="dashed")



#Multiplicative Holt Winters ESM
hw.m.energy <- hw(train, seasonal = "multiplicative", h=12)
summary(hw.m.energy)
#MAPE = 17.26937  
#MAE = 100.916 
#AIC = 2515.930 
#BIC = 2568.936 


#plot forecasts for mult
autoplot(hw.m.energy)+
  autolayer(fitted(hw.m.energy),series="Fitted")+ylab("Energy")+ geom_vline(xintercept = 2020.4166667,color="orange",linetype="dashed")

#The additive model has a slightly better MAPE by almost 1, but the error on the forecasts are
#much larger. I think the multiplicative will be more accurate in validation but both will be tested.


#Check automatic search procedure
ets.energy = ets(train)
summary(ets.energy)
#MAPE = 16.82263 

#Additive model still has best MAPE





###
###
### ARIMA
###
###


#Seeing both decompositions, seasonality is certainly present
#Use the Seasonal unit root test if a seasonal difference needs to be taken




#Test if seasonal difference needs to be taken
train %>% nsdiffs()

# Returns 1 difference needs to be taken



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

#Coninue modeling with ARIMA


#Check ACF
acf1 = Acf(train1, lag=30)$acf

#MA2 term

#Check PACF
pacf1 = Pacf(train1, lag=30)$acf

#no AR terms

#Check for white noise using Ljung Box
index1=seq(1,10)
White.LB <- rep(NA, 10)
for(i in 2:10){
  White.LB[i] <- Box.test(train1, lag=i, type="Ljung-Box", fitdf = 1)$p.value
}

white.dat=data.frame(cbind(White.LB,index1))
colnames(white.dat)=c("pvalues","Lag")

ggplot(white.dat,aes(x=factor(Lag),y=pvalues))+geom_col()+labs(title="Ljung-Box test p-values",x="Lags",y="p-values")+coord_cartesian(ylim = c(0, 0.025))

#All significant p-values. Hence we have some modeling to do.  










#Build ARIMA(0,0,2)
arima1 = Arima(train1, order=c(0,0,2))
summary(arima1)

#Plot residuals
ggtsdisplay(arima1$residuals)

#Check ACF and PACFs of residuals on model
acf2 = Acf(arima1$residuals, lag=24)
pacf2 = Pacf(arima1$residuals, lag=24)


#Check for white noise
index1=seq(1,10)
White.LB <- rep(NA, 10)
for(i in 4:10){
  White.LB[i] <- Box.test(arima1$residuals, lag=i, type="Ljung-Box", fitdf = 3)$p.value
}

white.dat=data.frame(cbind(White.LB,index1))
colnames(white.dat)=c("pvalues","Lag")

ggplot(white.dat,aes(x=factor(Lag),y=pvalues))+geom_col()+labs(title="Ljung-Box test p-values",x="Lags",y="p-values")+coord_cartesian(ylim = c(0, 0.025))

#Technically white noise but I think we can do better still






#Check autoarima
auto.mod = auto.arima(train)  
summary(auto.mod)
#Back to the drawing boards








#Try ARIMA(1,0,1)
arima2 = Arima(train1, order=c(1,0,1))
summary(arima2)

#Plot residuals
ggtsdisplay(arima2$residuals)

#Check ACF and PACFs of residuals on model
acf3 = Acf(arima2$residuals, lag=24)
pacf3 = Pacf(arima2$residuals, lag=24)


#Check for white noise
index1=seq(1,10)
White.LB <- rep(NA, 10)
for(i in 4:10){
  White.LB[i] <- Box.test(arima2$residuals, lag=i, type="Ljung-Box", fitdf = 3)$p.value
}

white.dat=data.frame(cbind(White.LB,index1))
colnames(white.dat)=c("pvalues","Lag")

ggplot(white.dat,aes(x=factor(Lag),y=pvalues))+geom_col()+labs(title="Ljung-Box test p-values",x="Lags",y="p-values")#+coord_cartesian(ylim = c(0, 0.025))



#Technically white noise but I think we can do better still








#Try some newly learned Seasonal ARIMA!!

#After reviewing the original ACF and PACF, I will try and ARIMA(0,0,0)(1,1,1)[12]
arima3 = Arima(train, order=c(0,0,0), seasonal = c(1,1,1))
summary(arima3)

#Plot residuals
ggtsdisplay(arima3$residuals)

#Check ACF and PACFs of residuals on model
acf4 = Acf(arima3$residuals, lag=24)
pacf4 = Pacf(arima3$residuals, lag=24)

#Check for white noise
index1=seq(1,10)
White.LB <- rep(NA, 10)
for(i in 4:10){
  White.LB[i] <- Box.test(arima3$residuals, lag=i, type="Ljung-Box", fitdf = 3)$p.value
}

white.dat=data.frame(cbind(White.LB,index1))
colnames(white.dat)=c("pvalues","Lag")

ggplot(white.dat,aes(x=factor(Lag),y=pvalues))+geom_col()+labs(title="Ljung-Box test p-values",x="Lags",y="p-values")+coord_cartesian(ylim = c(0, 0.025))

#No white noise here, try adding MA(2) term



#Try ARIMA(1,0,2)(1,1,1)[12]
arima4 = Arima(train, order=c(0,0,2), seasonal = c(1,1,1))
summary(arima4)

#Plot residuals
ggtsdisplay(arima4$residuals)
gghistogram(arima4$residuals) + ggtitle('Histogram of Residuals') + ylab('Count') + xlab('')

#Check ACF and PACFs of residuals on model
acf5 = Acf(arima4$residuals, lag=24)
pacf5 = Pacf(arima4$residuals, lag=24)

#Check for white noise
index1=seq(1,12)
White.LB <- rep(NA, 12)
for(i in 6:12){
  White.LB[i] <- Box.test(arima4$residuals, lag=i, type="Ljung-Box", fitdf = 5)$p.value
}

white.dat=data.frame(cbind(White.LB,index1))
colnames(white.dat)=c("pvalues","Lag")

ggplot(white.dat,aes(x=factor(Lag),y=pvalues))+geom_col()+labs(title="Ljung-Box test p-values",x="Lags",y="p-values")#+coord_cartesian(ylim = c(0, 0.025))

checkresiduals(arima4)

#White noise here, we may have a winner
#MAPE = 14.63323  





#Try ARIMA(1,0,2)(1,1,1)[12]
arima5 = Arima(train, order=c(1,0,2), seasonal = c(1,1,1))
summary(arima5)

#Plot residuals
ggtsdisplay(arima5$residuals)

#Check ACF and PACFs of residuals on model
acf6 = Acf(arima5$residuals, lag=24)
pacf6 = Pacf(arima5$residuals, lag=24)

#Check for white noise
index1=seq(1,12)
White.LB <- rep(NA, 12)
for(i in 6:12){
  White.LB[i] <- Box.test(arima5$residuals, lag=i, type="Ljung-Box", fitdf = 5)$p.value
}

white.dat=data.frame(cbind(White.LB,index1))
colnames(white.dat)=c("pvalues","Lag")

ggplot(white.dat,aes(x=factor(Lag),y=pvalues))+geom_col()+labs(title="Ljung-Box test p-values",x="Lags",y="p-values")#+coord_cartesian(ylim = c(0, 0.025))

#White noise here too but MAPE is a little higher. The prior model is my winner.
#MAPE = 14.65826   








#Final ARIMA model - ARIMA(0,0,2)(1,1,1)[12]
final.arima = arima4
backup.arima = arima5






#Final ESM mode 
final.esm = hw.m.energy
backup.esm = hw.a.energy








#Forecasts on validation data - ESM

#ESM error
error.esm=val-final.esm$mean

#Calculate MAE & MAPE
MAE.esm=mean(abs(error.esm))
MAPE.esm=mean(abs(error.esm)/abs(val))
#Better then additive


#ESM error - backup
error.esm1=val-backup.esm$mean

#Calculate MAE & MAPE
MAE.esm.b=mean(abs(error.esm1))
MAPE.esm.b=mean(abs(error.esm1)/abs(val))
#Worse then multicplicative




#Forecast using ARIMA model
arima.predicts = forecast::forecast(final.arima, h = 12)

#ARIMA error
error.arima=val-arima.predicts$mean

#Calculate MAE & MAPE
MAE.arima=mean(abs(error.arima))
MAPE.arima=mean(abs(error.arima)/abs(val))
#Not as good as multiplicative esm


#Forecast using ARIMA model - backup
arima.predicts1 = forecast::forecast(backup.arima, h = 12)

#ARIMA error
error.arima1=val-arima.predicts1$mean

#Calculate MAE & MAPE
MAE.arima.b=mean(abs(error.arima1))
MAPE.arima.b=mean(abs(error.arima1)/abs(val))
#worst


#Forecast using ARIMA model - auto
arima.predicts2 = forecast::forecast(auto.mod, h = 12)

#ARIMA error
error.arima2=val-arima.predicts2$mean

#Calculate MAE & MAPE
MAE.arima.c=mean(abs(error.arima2))
MAPE.arima.c=mean(abs(error.arima2)/abs(val))
#Bad do not use






#Vizualize ESM on full dataset
clrs <- c("blue", "red", "black")
autoplot(final.esm) + 
  autolayer(final.esm$mean, series='Forecasts') + 
  autolayer(fitted(final.esm), series='Model') + 
  autolayer(energy.ts, series='Raw') +
  guides(colour=guide_legend(title="Data series")) +
  scale_color_manual(values=clrs) + 
  ylim(c(110,1375)) 



#Vizualize ARIMA on full data set
autoplot(arima.predicts) + 
  autolayer(arima.predicts$mean, series='Forecasts') +
  autolayer(fitted(final.arima), series='Model') +
  autolayer(energy.ts, series='Raw') +
  guides(colour=guide_legend(title="Data series")) +
  scale_color_manual(values=clrs)+ 
  ylim(c(110,1375)) 






###
###Exporting data for visuals in report
###



#Export predictions from ESM and ARIMA
predicts.tab = rbind(val,final.esm$mean,arima.predicts$mean)
rownames(predicts.tab) = c('Validation', 'ESM', 'ARIMA')
colnames(predicts.tab) = c('Dec 2019', 'Jan 2020', 'Feb 2020', 'Mar 2020', 'Apr 2020', 'May 2020', 'Jun 2020', 'Jul 2020', 'Aug 2020', 'Sep 2020', 'Oct 2020', 'Nov 2020')

write.csv(predicts.tab, 'predictions.csv')


#Export raw series, predictions and intervals for ESM visualizations
preds1 = final.esm$mean
esm.lower = t(final.esm$lower)
colnames(esm.lower) = c('Dec 2019', 'Jan 2020', 'Feb 2020', 'Mar 2020', 'Apr 2020', 'May 2020', 'Jun 2020', 'Jul 2020', 'Aug 2020', 'Sep 2020', 'Oct 2020', 'Nov 2020')
esm.upper = t(final.esm$upper)
colnames(esm.lower) = c('Dec 2019', 'Jan 2020', 'Feb 2020', 'Mar 2020', 'Apr 2020', 'May 2020', 'Jun 2020', 'Jul 2020', 'Aug 2020', 'Sep 2020', 'Oct 2020', 'Nov 2020')

viz.tab1 = rbind(preds1, esm.lower, esm.upper, val)
rownames(viz.tab1) = c("Forecasts",'C.I. 80%', 'C.I. 95%', '80%', '95%', 'Validation')
viz.tab1

write.csv(viz.tab1, 'esm_viz.csv')


#Export raw series, predictions and intervals for ARIMA visualizations
preds2 = arima.predicts$mean
arima.lower = t(arima.predicts$lower)
colnames(arima.lower) = c('Dec 2019', 'Jan 2020', 'Feb 2020', 'Mar 2020', 'Apr 2020', 'May 2020', 'Jun 2020', 'Jul 2020', 'Aug 2020', 'Sep 2020', 'Oct 2020', 'Nov 2020')
arima.upper = t(arima.predicts$upper)
colnames(arima.upper) = c('Dec 2019', 'Jan 2020', 'Feb 2020', 'Mar 2020', 'Apr 2020', 'May 2020', 'Jun 2020', 'Jul 2020', 'Aug 2020', 'Sep 2020', 'Oct 2020', 'Nov 2020')

viz.tab2 = rbind(preds2, arima.lower, arima.upper, val)
rownames(viz.tab2) = c("Forecasts",'C.I. 80%', 'C.I. 95%', '80%', '95%', 'Validation')
viz.tab2

write.csv(viz.tab2, 'arima_viz.csv')



#Export raw series data for autoplot graph
energy.out = energy %>% mutate(Date = my(energy$Date))
write.csv(energy.out, 'energy_ts1.csv')


#Export residuals for diagnostic plts of ARIMA
write.csv(arima4$residuals, 'arima_resid.csv')
write.csv(Acf(arima4$residuals)$acf, 'acf_arima.csv')
out.acf = Acf(arima4$residuals)


