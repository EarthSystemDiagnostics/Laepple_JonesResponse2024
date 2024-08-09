
## Status/belongs to: Laepple et al., Matters arising, Nature 2024
## Output:Figure 3, Extended Data Figure 3 

basedrive=" " #Adapt to the path where the folder is downloaded e.g. "/Users/mustermensch/data/wais/"
setwd(paste(basedrive,"Laepple_JonesResponse2024",sep=""))

library(zoo)
library(tidyr)
library(ggplot2)

library(PaleoSpec) 
#PaleoSpec is an R package to assist in the spectral analysis of timeseries, 
#in particular timeseries of climate variables from observational, model, 
#and proxy paleoclimate data sources. PaleoSpec contains functions to analyse existing 
#timeseries and to generate timeseries with specific spectral properties.
#https://github.com/EarthSystemDiagnostics/paleospec

#Load various helper functions to handle time-series

source("./src/tools/ToolsTimeseries.R")

#from https://doi.org/10.15784/601603 WDC_seasonal_water_isotopes.xlsx 
data<-read.csv(file="./data/WDC_seasonal_water_isotopes.csv",skip=21,header=TRUE,sep=";",dec=",")

newTime <- (30:10900)

#Extract annual mean, summer and winter from the supplied dataset
avgAnnual<-binAvg(data[,5]*1000,data[,7],breaks=newTime+0.5)
avgSummer<-binAvg(data[,1]*1000,data[,2],breaks=newTime+0.5)
avgWinter<-binAvg(data[,3]*1000,data[,4],breaks=newTime+0.5)


dSummer<-pTs((avgSummer$avg-avgAnnual$avg),(newTime[-1])/1000) 
dWinter<-pTs((avgWinter$avg-avgAnnual$avg),(newTime[-1])/1000)

dSummer.1000<-rollmean(na.fill(dSummer),1000)
dWinter.1000<-rollmean(na.fill(dWinter),1000)
dSeasonalRange.1000 <- (dSummer.1000-dWinter.1000)/2

#https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05411-8/MediaObjects/41586_2022_5411_MOESM2_ESM.xlsx
#Read mean annual temperature from Cuffey via the supplied Source Data of Figure 2  of 
#Jones, T.R., Cuffey, K.M., Roberts, W.H.G., Markle, B.R., Steig, E.J., Stevens, C.M., 
#Valdes, P.J., Fudge, T.J., Sigl, M., Hughes, A.G., Morris, V., Vaughn, B.H., Garland, J., Vinther, B.M., Rozmiarek, K.S., 
#Brashear, C.A., White, J.W.C., 2023. Seasonal temperatures in West Antarctica during the Holocene. Nature 613, 292–297. https://doi.org/10.1038/s41586-022-05411-8

#Column 8= Age; Column 9 Mean Temp, 1,000-yr mean
#
temp<-read.csv("./data/Figure2PanelC.csv",sep=";",dec=",",skip=4,header=FALSE)
index <- !is.na(temp[,8])

STIME <- time(dSeasonalRange.1000)[1]
ETIME <- time(dSeasonalRange.1000)[length(time(dSeasonalRange.1000))]

cuffey.ts.1000<-window(pTs(temp[index,9],temp[index,8]),STIME,ETIME)

#https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05411-8/MediaObjects/41586_2022_5411_MOESM2_ESM.xlsx
#Read summer and winter reconstructed and modeled temperatures from Jones Figure 2
#Source Data of Figure 2  of 
#Jones, T.R., Cuffey, K.M., Roberts, W.H.G., Markle, B.R., Steig, E.J., Stevens, C.M., 
#Valdes, P.J., Fudge, T.J., Sigl, M., Hughes, A.G., Morris, V., Vaughn, B.H., Garland, J., Vinther, B.M., Rozmiarek, K.S., 
#Brashear, C.A., White, J.W.C., 2023. Seasonal temperatures in West Antarctica during the Holocene. Nature 613, 292–297. https://doi.org/10.1038/s41586-022-05411-8


figure.summer<-read.csv("./data/Figure2PanelA.csv",sep=";",dec=",",skip=4,header=FALSE)
figure.winter<-read.csv("./data/Figure2PanelWinter.csv",sep=";",dec=",",skip=4,header=FALSE)

#Convert to time-series objects
j23.summer <- window(pTs(figure.summer[,2],figure.summer[,1]),STIME,ETIME)
j23.winter <- window(pTs(figure.winter[,2],figure.winter[,1]),STIME,ETIME)

#Function to calculate anomalies relative to 1kyr BP, here 0.9-1.1kyrBP are choosen
n1k <- function(x) return(x-mean(window(x,0.9,1.1)))

### plot Figure 2, top panel
quartz(width=6,height=4.5)
plot(j23.summer,ylim=c(-5,2),xlim=c(11,0),xlab="Age (ka)",lwd=2,col="red",ylab="T Anom (K)",lty=3)
lines(n1k((dSeasonalRange.1000/6.96*5)+cuffey.ts.1000),col="red",lwd=3)

lines(figure.summer[,8],figure.summer[,9],col="grey",type="b",pch=19,cex=0.5) #ORBIT
#lines(figure.summer[,11],figure.summer[,12],col="violet",lty=2,cex=0.5) #MEBM
lines(figure.summer[,14],figure.summer[,15],col="darkgreen",type="b",pch=19,cex=0.5) #GLAC1D
lines(figure.summer[,17],figure.summer[,18],col="orange",type="b",pch=19,cex=0.5) #Ice6G
legend("bottomright",lwd=c(2,3),lty=c(3,1),col=c("red","red"),c("as published by J23","as implied assuming constant loss"),bty="n",)
legend("topleft",lwd=c(2,3),col=c("grey","darkgreen","orange"),pch=10,c("ORBIT","GLAC1D","ICE-6G"),bty="n")


### plot Figure 2, lower panel
quartz(width=6,height=4.5)
plot(j23.winter,ylim=c(-2,5),xlim=c(11,0),,xlab="Age (ka)",lwd=2,col="blue",ylab="T Anom (K)",lty=3)
lines(n1k(cuffey.ts.1000-(dSeasonalRange.1000/6.96*5)),col="blue",lwd=3)

lines(figure.winter[,8],figure.winter[,9],col="grey",type="b",pch=19,cex=0.5) #ORBIT
lines(figure.winter[,14],figure.winter[,15],col="darkgreen",type="b",pch=19,cex=0.5) #GLAC1D
lines(figure.winter[,17],figure.winter[,18],col="orange",type="b",pch=19,cex=0.5) #Ice6G
legend("bottomright",lwd=c(2,3),lty=c(3,1),col=c("blue","blue"),c("as published by J23","as implied assuming constant loss"),bty="n",)
legend("topright",lwd=c(2,3),col=c("grey","darkgreen","orange"),pch=10,c("ORBIT","GLAC1D","ICE-6G"),bty="n")


#Extended Data Figure 2:
## Read the 100yr average accumulation (3rd and 4th column)
##White, J., Bradley, E., Garland, J., Jones, T., Morris, V., Price, M., Vaughn, B., 2019. Stable Isotopes of Ice in the Transition and 
#Glacial Sections of the WAIS Divide Deep Ice Core. https://doi.org/10.15784/601274 from WDC_accumulation_annual.xls downloaded 03/2024
temp<-as.matrix(read.csv("./data/WDC_accumulation_50yr.csv",sep=";",dec=",",skip=7,header=FALSE))
index<-!is.na(temp[,4])
accum.ts<-pTs(rollmean(temp[index,5],10),rollmean(temp[index,4],10)/1000) #make a 10-point, thus effectivly 1000yr running average

accum.ts1000<-approx(c(time(accum.ts)),c(accum.ts),(time(dSeasonalRange.1000)))$y
cor.test(accum.ts1000,dSeasonalRange.1000) #Correlation = 0.90


#Figure plotting
quartz(width=6.5,height=4.5)
plot(dSeasonalRange.1000,lwd=2,,xlab="Age (ka)",ylab="dD amplitude", xlim=c(11,0),xaxs = 'i',yaxs = 'i',main="")
par(new=TRUE)
plot(accum.ts,col="blue",type="l",lwd=2,xlim=c(11,0),ylim=c(0.17-0.02,0.27+0.02),xaxs = 'i',yaxs = 'i',axes = FALSE, xlab ="", ylab ="")
axis(4,col="blue",col.lab="blue",col.axis="blue")
legend("topleft",col=c("blue","black"),lwd=2,c("Accumulation rate","seasonal amplitude in dD"))





