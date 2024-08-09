## Estimate loss of the seasonal amplitude in the WAIS record
## Status/belongs to: Laepple et al., Matters arising, Nature 2024
## Output:Figure 1, Extended Data Figure 1 and the loss estimates provided in the text

basedrive="" #Adapt to the path where the folder is downloaded e.g. "/Users/mustermensch/data/wais/"
setwd(paste(basedrive,"Laepple_JonesResponse2024",sep=""))

getwd()

library(zoo) 
library(PaleoSpec) 
#PaleoSpec is an R package to assist in the spectral analysis of timeseries, 
#in particular timeseries of climate variables from observational, model, 
#and proxy paleoclimate data sources. PaleoSpec contains functions to analyse existing 
#timeseries and to generate timeseries with specific spectral properties.
#https://github.com/EarthSystemDiagnostics/paleospec



#Load various helper functions to handle time-series
source("./src/tools/ToolsTimeseries.R")
#Load diffusion code for the diffusion correction
source("./src/tools/diffusion.R")


#from https://doi.org/10.15784/601603 WDC_seasonal_water_isotopes.xlsx 
data<-read.csv(file="./data/WDC_seasonal_water_isotopes.csv",skip=21,header=TRUE,sep=";",dec=",")

time.fill <- na.fill(data[,5]) #To allow a time-series object; Data get's not interpolated, only the time
wais.pTs <- pTs(cbind(data[,6],data[,7]),time.fill*1000,lat=rep(-79.468,2),lon=rep(-112.08,2),c("raw","corrected"))

#Daily ERA Interim T2M data extracted at the WAIS position lat=-79.468,lon=-112.08
#starting 1.1.1979
t.wais<-read.table(file="./data/t2mWAISeraInterim.txt",skip=2,header=FALSE)[,1]
t.wais.ts <- pTs(t.wais,seq(t.wais)/365+1979) ##Set the right time axis

#White, J., Bradley, E., Garland, J., Jones, T., Morris, V., Price, M., Vaughn, B., 2019. 
#Stable Isotopes of Ice in the Transition and Glacial Sections of the WAIS Divide Deep Ice Core. 
#https://doi.org/10.15784/601274 downloaded 03/2024
#WDC06A_age_dD_d18O.txt
#[“Age”] - yrs before 1950
#["dDcmave "]- 1/20th yr interpolated dD data
#["d18ocmave"] - 1/20th yr interpolated d18O data

data<-read.table("./data/WDC06A_age_dD_d18O.txt",skip=2)
dD.wais.ts<-pTs(rev(data[,2]),rev(1950-data[,1]))

### Figure 1
quartz(width=8,height=5.5)
plot(scale(rollmean(t.wais.ts,1),scale=FALSE),xlim=c(1979,2005),col="grey",lwd=0.5,xlab="Time (year)",ylab="Temperature anomaly (C)",ylim=c(-30,30))

lines(scale(rollmean(t.wais.ts,20),scale=FALSE),col="black",lwd=2)
lines(scale(window(dD.wais.ts,1979,2015)/6.96,scale=FALSE),col="red",lwd=2)

legend("topleft",lwd=2,col=c("black","red"),c("T reanalysis diffused","T isotope"),bty="n")


##### Caclulate Diffusion correction... 
#diffuse white noise with varying diffusion lengths modeled for the 
#upper part / section of interest of the WAIS core; replicate; take the expected = mean spectra


out<-TemporalDiffusionLength(rho.surface = 387,bdot=220,T=273.15-31.1,t.res=1/100,nt=(2005-1979+1)*100,dD=TRUE,P=794)
freq<-spectrum(ts(DiffuseRecord(rnorm(length(out[,2])),out[,2]),deltat=0.01),plot=FALSE)$freq
specDiffused<-replicate(1000,spectrum(ts(DiffuseRecord(rnorm(length(out[,2])),out[,2]*100),deltat=0.01),plot=FALSE)$spec)
transfer<-list(freq=freq,ratio=0.01/rowMeans(specDiffused))

##Calculate the spectra for Extended Data Figure 1
spec.dD.present<-SpecMTM(window(dD.wais.ts,1979,2005)/6.96)

#Correct by multiplying with the transfer function 
spec.dD.present.undiffused <- spec.dD.present
spec.dD.present.undiffused$spec = spec.dD.present$spec * approx(transfer$freq,transfer$ratio,spec.dD.present$freq)$y

spec.dD.past<-SpecMTM(na.fill(window(wais.pTs[,2],5000,5026))/6.96) ##Raw
spec.dD.past.diff<-SpecMTM(na.fill(window(wais.pTs[,1],5000,5026))/6.96) ##Diffused

spec.T<-SpecMTM(window(t.wais.ts,1979,2005))

quartz(width=8,height=5.5)
LPlot(AddConfInterval(spec.T),col="black",xlim=c(0.1,5),lwd=2,xlab="Frequency (1/year)",ylab="PSD Temperature",ylim=c(1e-5,1e3))
abline(v=c(1/1.2,1.2),lwd=2,lty=2)
LLines(AddConfInterval(spec.dD.present.undiffused),col="red")
LLines(AddConfInterval(spec.dD.past),col="brown")
LLines(AddConfInterval(spec.dD.past.diff),col="orange")


legend("bottomleft",bty="n",col=c("black","red","brown","orange"),lwd=2,c("Reanalysis 1979-2005","WAIS core 1979-2005 diffusion corrected","WAIS core 5ka BP, diffusion corrected","WAIS core 5ka BP, raw"))

#The ratio of the variance in the frequency band around the annual cycle (1/1.2-1.2) year-1: 
1/(sqrt(GetVarFromSpectra(spec.T,c(1/1.2,1.2))$var/GetVarFromSpectra(spec.dD.present.undiffused,c(1/1.2,1.2))$var))  #0.22 (1979-2005)

#Sensitivity on the frequency range
1/(sqrt(GetVarFromSpectra(spec.T,c(1/1.1,1.1))$var/GetVarFromSpectra(spec.dD.present.undiffused,c(1/1.1,1.1))$var)) #0.19
1/(sqrt(GetVarFromSpectra(spec.T,c(1/1.3,1.3))$var/GetVarFromSpectra(spec.dD.present.undiffused,c(1/1.3,1.3))$var)) #0.24

#Sensitiviy on the time-period; here looking at a time period 5kyr BP instead of the modern data
1/(sqrt(GetVarFromSpectra(spec.T,c(1/1.2,1.2))$var/GetVarFromSpectra(spec.dD.past,c(1/1.2,1.2))$var)) #0.21^2 (5ka BP diffusion corrected) 


###

#To test the sensitivity on the method, we also estimate the loss in the time domain
#Dataset from Jones, T.R., 2022. Seasonal temperatures in West Antarctica during the Holocene. https://doi.org/10.15784/601603

data<-read.csv(file="./data/WDC_seasonal_water_isotopes.csv",skip=21,header=TRUE,sep=";",dec=",")

newTime <- (30:1029) #look at the most recent millenia of the dataset
avgAnnual<-binAvg(data[,5]*1000,data[,7],breaks=newTime+0.5)
avgSummer<-binAvg(data[,1]*1000,data[,2],breaks=newTime+0.5)
avgWinter<-binAvg(data[,3]*1000,data[,4],breaks=newTime+0.5)

dSummer<-pTs((avgSummer$avg-avgAnnual$avg),newTime[-1]-0.5)  
dWinter<-pTs((avgWinter$avg-avgAnnual$avg),newTime[-1]-0.5)

(mean(c(dSummer),na.rm=TRUE)-mean(c(dWinter),na.rm=TRUE))/2 #Mean amplitude 15.9permil

t.wais.matrix <- c(t.wais.ts)[1:(37*365)]
dim(t.wais.matrix) <- c(365,37)
t.clim<-rowMeans(t.wais.matrix)

diff(range(t.clim))/2 #Mean amplitude of the ERA data at the WAIS site is 11.6K

