## Calculate the diffusion length that would be implied from the Jones data if the seasonality over the Holocene did not change
## Status/belongs to: Laepple et al., Matters arising, Nature 2024
## Output:Figure 3, Extended Data Figure 3 

basedrive=" " #Adapt to the path where the folder is downloaded e.g. "/Users/mustermensch/data/wais/"
setwd(paste(basedrive,"Laepple_JonesResponse2024",sep=""))

library(tidyr)
library(ggplot2)
library(dplyr)
library(zoo)
library(PaleoSpec) 
#PaleoSpec is an R package to assist in the spectral analysis of timeseries, 
#in particular timeseries of climate variables from observational, model, 
#and proxy paleoclimate data sources. PaleoSpec contains functions to analyse existing 
#timeseries and to generate timeseries with specific spectral properties.
#https://github.com/EarthSystemDiagnostics/paleospec

#Load various helper functions to handle time-series
source("./src/tools/ToolsTimeseries.R")



#### 


##' @title Return the variance in a frequency band
##' @param middleYear mid-point of the window
##' @param data pTs with two timeseries
##' @param f frequency range 
##' @param deltat width of the window
##' @return vector (variance of ts1, variance of ts2)
##' @author Thomas Laepple 
GetVariance <- function(middleYear,data,f=c(1/1.2,1.2),deltat=500)
{
    print(middleYear)
    ts <- window(data,middleYear-deltat/2,middleYear+deltat/2)
    s1<-SpecMTM(na.omit(na.fill(ts[,1])))
    s2<-SpecMTM(na.omit(na.fill(ts[,2])))        
    return(c(GetVarFromSpectra(s1,f=f)$var,GetVarFromSpectra(s2,f=f)$var))
}


#Data downloaded from Jones, T.R., 2022. Seasonal temperatures in West Antarctica during the Holocene. https://doi.org/10.15784/601603

data<-read.csv(file="./data/WDC_seasonal_water_isotopes.csv",skip=21,header=TRUE,sep=";",dec=",")
time.fill <- na.fill(data[,5]) #To allow a time-series object; Data get's not interpolated, only the time

wais.pTs <- pTs(cbind(data[,6],data[,7]),time.fill*1000,lat=rep(-79.468,2),lon=rep(-112.08,2),c("raw","corrected"))


##Plot Extended Data Figure 3
## Read the diffusion length provided by J23 in their Extended Data Figure 1 Source data

diffusion.jones23 <- read.csv("data/Jones23_Diffusion.csv",skip=10,sep=";",dec=",")
colnames(diffusion.jones23) <- c('age','sigma','sigmaplus1sd','sigmaminus1sd')

## and at the same midpoints, calculate the diffusion length from the variance around the annual cycle
midpoint <- diffusion.jones23$age*1000

# XXX Error in pTs(res, stats::time(x)[p1], lat = atb$lat, lon = a.... 'lat' and 'lon' need to be of the same length ----- could this also be related to the ts vs pts objects?
result<-sapply(midpoint,GetVariance,data=wais.pTs,deltat=140,f=c(1/1.2,1.2))

ratio.constant <- result[1,]/mean(result[2,])
ratio <- result[1,]/result[2,]

sigma.implied.constant<-sqrt(-log(ratio.constant))/(2*pi*1)
sigma.implied<-sqrt(-log(ratio))/(2*pi*1)

quartz(width=7,height=5)
plot(diffusion.jones23$age,diffusion.jones23$sigma,col="black",type="l",xlab="Age (ka)",ylab="diffusion length dD (yr)",xlim=c(11,0),ylim=c(0.16,0.34))
polygon(c(diffusion.jones23$age, rev(diffusion.jones23$age)), c(diffusion.jones23$sigmaminus1sd, rev(diffusion.jones23$sigmaplus1sd)),  col = ColTransparent('black', 0.3), border = NA)

lines(midpoint/1000,sigma.implied,col="red")

legend("bottomleft",col=c("black","red"),lwd=1,bty="n",c("J23","inferred from ratio of annual variance"))



## For Figure 3, the same on 500yr to get the same resolution as Kahle et al.,2018

#Read Figure data from Kahle et al., 2012, Fig. 11, right panel. 
#Kahle, E.C., Holme, C., Jones, T.R., Gkinis, V., Steig, E.J., 2018. 
#A Generalized Approach to Estimating Diffusion Length of Stable Water Isotopes 
#From Ice-Core Data. Journal of Geophysical Research: Earth Surface 123, 2377–2391. 
#https://doi.org/10.1029/2018JF004764
#As neither of the authors Kahle, Holme or Gkinis could supply the source data for 
#the figure, it was digitised from the pdf figure. The digitization error, quantified 
#by replication, is negligible compared to the method uncertainty.

kahle2018<-read.table("./data/Kahle2018_Fig11right_digitized.txt",skip=3,header=TRUE)

midpoint <- seq(from=400,to=10900,by=250)

result<-sapply(midpoint,GetVariance,data=wais.pTs,deltat=500)
ratio.constant <- result[1,]/mean(result[2,])
ratio <- result[1,]/result[2,]

sigma.implied.constant<-sqrt(-log(ratio.constant))/(2*pi*1)
sigma.implied<-sqrt(-log(ratio))/(2*pi*1)

## To get the diffusion lengths in the depth domain we need the layer thickness
#Data from Sigl, M., Buizert, C., Fudge, T.J., Winstrup, M., Cole-Dai, J., 
#McConnell, J.R., Ferris, D.G., Rhodes, R.H., Taylor, K.C., Welten, K.C., Woodruff, T.E., Adolphi, 
#F., Baggenstos, D., Brook, E.J., Caffee, M.W., Clow, G.D., Cheng, H., Cuffey, K.M., Dunbar, N.W., Edwards, 
#R.L., Edwards, L., Geng, L., Iverson, N., Koffman, B.G., Layman, L., Markle, B.R., Maselli, O.J., McGwire,
#K.C., Muscheler, R., Nishiizumi, K., Pasteris, D.R., Severinghaus, J.P., Sowers, T.A., Steig, E.J., 2019.
#WAIS Divide Deep ice core 0-68 ka WD2014 chronology. https://doi.org/10.1594/PANGAEA.902577
chron<- read.table("./data/WD2014_Chronology_BuizertSigl.tab",skip=54,sep="\t")
layerthickness <- diff(chron[,1])/diff(chron[,2])/1000 ##Layer thickness in meter w.e. per year
layerthickness.ts<-MakeEquidistant(chron[-1,2]*1000,layerthickness,time.target=midpoint)

## Merge datasets here that they can get one legend!
simulated <- tibble(age=(midpoint/1000)[],sigma_constant=(sigma.implied.constant*layerthickness.ts)[],sigma=(sigma.implied*layerthickness.ts)[])


###Table gymnastics to prepare the plotting

kahle2018.long <- kahle2018 %>%
  pivot_longer(cols = -age,
               names_to = c(".value", "set"),
               names_pattern = "(.*?)\\.(.*)")


MAXAGE <- 11

# Filter data for ages less than MAXAGE
filtered_df <- kahle2018.long %>%
  filter(age < MAXAGE) %>%
  pivot_wider(names_from = set, values_from = c(mean, lower, upper))

# Calculate standard deviation across the three series at each age point
filtered_df <- filtered_df %>%
  rowwise() %>%
  mutate(sd_across_series = sd(c(mean_Double, mean_Noise, mean_Jones), na.rm = TRUE))

# Calculate the pooled standard deviation
pooled_sd <- sqrt(sum((3 - 1) * (filtered_df$sd_across_series ^ 2)) / (3 * nrow(filtered_df) - 3))

# Calculate the mean of the three time series
filtered_df <- filtered_df %>%
  rowwise() %>%
  mutate(mean_of_series = mean(c(mean_Double, mean_Noise, mean_Jones), na.rm = TRUE))

meanseries <- filtered_df %>% select(age,mean_of_series) %>%  mutate(mean=mean_of_series,set = "mean_of_series")

## Plot Figure 3
ggplot(kahle2018.long, aes(x = age, y = mean, color = set)) +
  geom_line() +
  xlim(11, 0) +
  theme_minimal() +
  labs(x = "Age (ka)", y = expression(sigma[depth] ~ delta*D ~ (m))) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  theme(legend.position = c(0.70, 0.17)) +
  ggtitle("") +
  scale_color_manual(
    values = c(
      "Double" = "blue",
      "Noise" = "cyan",
      "Jones" = "red",
      "diffusion length for no change in Holocene seasonality" = "black"
    ),
    labels = c(
      "Double" = "Double-Gaussian Model (Kahle et al., 2018)",
      "Noise" = "Noise-Adding Technique (Kahle et al., 2018)",
      "Jones" = "First Order Model (Jones et al., 2017)",
      "diffusion length for no change in Holocene seasonality" = "diffusion length for no change in Holocene seasonality"
    ),
    breaks = c("diffusion length for no change in Holocene seasonality", "Double", "Noise", "Jones")
  ) +
  scale_fill_manual(
    values = c("Estimation Uncertainty" = "grey"),
    labels = c(
      "Estimation Uncertainty" = "estimation uncertainty"
    )
  ) +
  labs(colour = "", fill = "") +
  geom_line(data = simulated, aes(x = age, y = sigma_constant, color = "diffusion length for no change in Holocene seasonality"), lwd = 1.5) +
  geom_ribbon(data = meanseries, aes(x = age, ymin = mean_of_series - 2 * pooled_sd, ymax = mean_of_series + 2 * pooled_sd, fill = "Estimation Uncertainty"),
              alpha = 0.4,colour="NA")


