## Helper functions extracted from the pfields library to handle time-series and fields tlaepple@awi.de


addhistory<-function(x,newhist)
{
	newhist<-paste(date(),newhist)
	attr(x,"history")<-c(attr(x,"history"),newhist)
	return(x)
}




#Coordinate conversion from 1D<->2D
c1t2<-function(x,nLon)
{
      x<-x-1
	lat<-x%/%nLon+1
	lon<-x%%nLon+1
	return(list(lat=lat,lon=lon))
}

c2t1<-function(lat,lon,nLon)
{
	return(nLon*(lat-1)+(lon))
}

mergeattr <- function(data,source1,source2,newhistory='')
{
	result<-data
	temp1<-attributes(source1)
	temp2<-attributes(source2)
	attr(result,'lat')<-c(temp1$lat,temp2$lat)
	attr(result,'lon')<-c(temp1$lon,temp2$lon)
	attr(result,'name')<-c(temp1$name,temp2$name)
	attr(result,'history')<-c(temp1$history,paste(date(),newhistory))
	return(result)

}


copyattr <- function(data,source,newhistory='',cclass=TRUE)
{

	temp<-attributes(source)
	attr(data,'lat')<-temp$lat
	attr(data,'lon')<-temp$lon
	attr(data,'name')<-temp$name
	attr(data,'history')<-c(temp$history,paste(date(),newhistory))
	if (cclass) class(data)<-class(source)
	return(data)

}



is_pTs <- function(data) (sum(class(data) == 'pTs')>0)
is_pField <- function(data) (sum(class(data) == 'pField')>0)

summary.pTs <- function(x, ...)
{
 temp<-attributes(x)
 print('Proxy timeseries object')
 print(paste('Names: ',paste(temp$name,collapse=' / ')))
 print('History')
 print(temp$history)
 print("")
 cat("Time range: ",min(time(x))," - ",max(time(x)), "N:",length(time(x)),"\n")
 cat("Data range: ",min(x)," - ",max(x),"\n")
}

summary.pField <- function(x, ...)
{
 temp<-attributes(x)
 print('Proxy field object')
 print(paste('Names: ',paste(temp$name,collapse=' / ')))
 print('History')
 print(temp$history)
 print("")
 cat("Time range: ",min(time(x))," - ",max(time(x)), "N:",length(time(x)),"\n")
 cat("Data range: ",min(x)," - ",max(x),"\n")
 print("spatial extent ")
 cat('lat: ',min(temp$lat)," - ",max(temp$lat),"N:",length(temp$lat),"\n")
 cat('lon: ',min(temp$lon)," - ",max(temp$lon),"N:",length(temp$lon),"\n")


}

#apply a function on fields containing complete NA sets...
na.apply<-function(x,FUN,... )
  {
    index<-!is.na(colSums(x))
     x[,index]<-FUN(x[,index], ...)
   return(x)
  }



getlat <- function(data) return(attr(data,"lat"))
getlon <- function(data) return(attr(data,"lon"))
getname <- function(data) return(attr(data,"name"))
gethistory <- function(data) return(attr(data,"history"))



#apply FUN(field->scalar) for each timestep and gives back a timeseries
applyspace<-function(data,FUN)
{
     index<-!is.na(colSums(data))
   ts<-apply(data[,index],1,FUN)
   return(pTs(ts,time(data),name=getname(data)))
       }


#apply FUN(field->scalar) for each gridbox and gives back a single field
applytime<-function(data,FUN,newtime=NULL)
{
   if (is.null(newtime)) newtime<-mean(time(data))
   field<-apply(data,2,FUN)
   return(pField(field,newtime,getlat(data),getlon(data),name=getname(data)))
}


#return 2D Fields filled with lats and lons

latlonField <- function(data)
{
  lat<-getlat(data)
  lon<-getlon(data)

  nlat<-length(lat)
  nlon<-length(lon)

  lon2d<-rep(lon,nlat)
  lat2d<-rep(lat,each=nlon)


  return(list(lat2d=lat2d,lon2d=lon2d))
}



rollmean.pTs <- function(x, k, na.pad = TRUE, align = c("center", "left", "right"), ...)
{
  return(applyData(x,rollmean,k,na.pad, align, ...))
}

applyData<-function(x,fun,... )
  {
    x[]<-fun(as.vector(x),... )
    return(x)  }



## Converts a list of single pTs timeseries to one pTs object
## x = list containing the pTs objects; all need to have the same length
## returns a pTs object with all the timeseries of x, including the
#lat/lon and name information
list2pTs<-function(x)
{
    TOLERANCE = 0.01 #tolerance for different time steps

    N<-length(x) #Number of timeseries in the list
    newTime<-time(x[[1]])
    names<-lapply(x,getname)  #get all anmes
    lat<-lapply(x,getlat)   #get the latitudes
    lon<-lapply(x,getlon) #get the longitudes

    result<-pTs(NA,newTime,lat,lon,names)
    for (i in 1:length(x))
    {
       if  (sum((newTime-time(x[[i]]))^2) > TOLERANCE) stop("time steps
are different in the different timeseries")
       result[,i]<-x[[i]]
    }
return(result)
}


#scale each timeseries=point of the field

#Here only the history gets modified by removing the last 1 / 2 entries and
#added by Ops.pField and adding a new one...

scale.pField<-function(x, center = TRUE, scale = TRUE)
{
  
  x[]<-scale(unclass(x),center=center,scale=scale)
  hist<-attr(x,"history")
  hist<-hist[1:(length(hist)-sum(center,scale))]
  hist<-c(hist,paste(date(),"scale center=",center,"scale=",scale))
  attr(x,"history")<-hist
  return(x)
}


#
scale.pTs <- function(x, center = TRUE, scale = TRUE)
{
  hist<-gethistory(x)
  hist<-hist[1:(length(hist)-1)]
  result<-pTs(NextMethod(),time(x),getlat(x),getlon(x),getname(x),hist,date=FALSE)
  hist<-c(hist,paste(date(),"scale center=",center,"scale=",scale))
  return(addhistory(result,hist))
}



detrend <- function(x, ...) UseMethod("detrend")

detrend.default <- function(x, ...) print("No default detrend function")


detrend.pField<-function(x)
{
  nonmissing<-!is.na(x[,1])
  x[nonmissing,]<-lm(unclass(x)~seq(nrow(x)))$residuals
  return(addhistory(x,"detrend"))
}






detrend.pTs<- function(x)
{
  x[]<-lm(unclass(x)~seq(length(x)))$residuals
  return(addhistory(x,"detrend"))
  
}



rollmean.pField <- function(x,k,na.pad = FALSE, ...)
{
  x[]<-rollmean.default(x,k,na.pad=TRUE, ...)
  if (na.pad==FALSE) return(addhistory(na.omit(x),paste("rollmean",k)))
  else return(addhistory(x,paste("rollmean",k)))
  
}


rollmean.pTs <- function(x,k,na.pad = FALSE, ...)
{
  x[]<-rollmean.default(x,k,na.pad=TRUE, ...)
  if (na.pad==FALSE) return(addhistory(na.omit(x),paste("rollmean",k)))
  else return(addhistory(x,paste("rollmean",k)))
  
}



na.fill<-function(x,rule=1)
  #Fill a the missing values in a pTs timeseries by linear interpolation
  #Inputs
  #x: pTs timeseries
  #rule: interpolation rule used in approx
  #Outputs: interpolated pTs timeseries
{
  if (sum(is.na(x))==0) return(x)
  nonMissing<-!is.na(x)
  return(pTs(approx(c(time(x))[nonMissing],c(x)[nonMissing],c(time(x)),rule=rule)$y,time(x),getlat(x),getlon(x)))
}


#data[lat,lon,time]
pField <- function(data,time,lat=0,lon=0,name=" ",history=" ",date=TRUE)
{
  #constants
  TOL=1/400 #tolerance less than 1 day
  
  #check data
  if (length(data) <=1 ) 
    if (is.null(data[1])) data<-rep(NA,length(lat)*length(lon)*length(time))
  else data<-rep(data[1],length(lat)*length(lon)*length(time))
  
  if (length(data) != (length(lat)*length(lon)*length(time))) 
  {	
    stop("nLat*nLon*nTime != N_Elements(data), if you want to create an empty field, supply data=NA")	
  }
  
  #shape data in 2D array
  dim(data)<-c(length(lat)*length(lon),length(time))
  
  if (length(time) > 1) #real time series
  {
    if (abs(max(diff(time))-min(diff(time)))>TOL) stop("time steps are not equidistant")
    result<-ts(t(data),start=time[1],deltat=(time[2]-time[1]))
  }
  else result<-ts(t(data),start=time[1]) #or only one time step
  
  #put attributes and classes
  attr(result,'lat')<-lat
  attr(result,'lon')<-lon
  attr(result,'name')<-name
  if (date) attr(result,'history')<-paste(date(),history)
  else attr(result,'history')<-paste(history)
  
  
  attr(result,'oclass')<-class(result)
  attr(result,'nclass')<-c('pField','ts')
  class(result)<-attr(result,'nclass')
  
  invisible(result)
}





#data[x(t),i]

pTs<-function(data,time,lat=0,lon=0,name="",history="",date=TRUE)
{
  #constants
  TOL=1/400 #tolerance less than 1 day
  
  if (length(data) <=1 ) 
    if (is.null(data[1])) data<-matrix(NA,length(time),length(name))
  else if (length(name) == 1) data<-rep(data[1],length(time)) else
    data<-matrix(data[1],length(time),length(name))
  
  
  
  if (!is.null(ncol(data)) && (ncol(data)>1)) #multiple datasets
  {       #check data
    if (nrow(data) != length(time)) stop("nTime != N_Elements(data)")
    
  }
  #check data
  else if (length(data) != (length(time))) stop("nTime != N_Elements(data)")
  
  #shape data in 2D array
  
  if (length(time) > 1) #real time series
  {
    if (abs(max(diff(time))-min(diff(time)))>TOL) stop("time steps are not equidistant")
    result<-ts(data,start=time[1],deltat=(time[2]-time[1]))
  }
  else result<-ts(data,start=time[1]) #or only one time step
  
  #put attributes and classes
  attr(result,'lat')<-lat
  attr(result,'lon')<-lon
  attr(result,'name')<-name
  
  if (date) attr(result,'history')<-paste(date(),history)
  else attr(result,'history')<-paste(history)
  
  
  attr(result,'oclass')<-class(result)
  attr(result,'nclass')<-c('pTs','ts')
  class(result)<-attr(result,'nclass')
  
  invisible(result)
}
Ops.pField <- function (e1, e2) 
{
  LIMIT<-1
  TLIMIT<-1
  #Field and scalar, scalar and Field... simply take the normal operator... and change history
  #Field and Field: check if lat/lons are compatible... than transform in time series and back in pField
  #.Generic gives the the name of the operator
  
  nField <- 0
  if (is_pField(e1)) { nField<-nField+1;name1<-attr(e1,"name")}
  else name1<-e1[1]
  
  if (is_pField(e2)) {nField<-nField+1;name2<-attr(e2,"name")}
  else name2<-e2[1]
  
  newname<-paste(name1,.Generic,name2)
  
  if (nField == 2) 
  {
    temp1<-attributes(e1)
    temp2<-attributes(e2)
    time1<-time(e1)
    
    if ((length(temp1$lat) != length(temp1$lat)) |(length(temp1$lon) != length(temp1$lon))) stop("grids not compatible")
    if ((sum(abs(temp1$lat-temp2$lat))+sum(abs(temp1$lon-temp2$lon))) > LIMIT) stop("grids not compatible")
    
    if (abs(sum(c(time(e1))-c(time(e2))))>TLIMIT) 	warning("Operator applied on two timebases, new timebase = first timebase")
    
    #Trick, bring it on the same timebase... start as to be > 0 due to R-bug
    
    
    e1<-ts(e1,start=100)
    e2<-ts(e2,start=100)
    return(pField(t(NextMethod()),time1,temp1$lat,temp1$lon,newname,c(paste("A:",temp1$history),paste("B:",temp2$history),newname)))
  }
  else
  {
    e<-NextMethod(.Generic)
    if (!is.null(attr(e,'nclass'))) class(e)<-attr(e,'nclass')
    return(addhistory(e,newname))
  }
}




Ops.pTs <- function (e1, e2) 
{
  TLIMIT<-1
  #Field and scalar, scalar and Field... simply take the normal operator... and change history
  #Field and Field: check if lat/lons are compatible... than transform in time series and back in pField
  #.Generic gives the the name of the operator
  
  nField <- 0
  if (is_pTs(e1)) { nField<-nField+1;name1<-attr(e1,"name")}
  else name1<-e1[1]
  
  if (is_pTs(e2)) {nField<-nField+1;name2<-attr(e2,"name")}
  else name2<-e2[1]
  
  newname<-paste(name1,.Generic,name2)
  
  if (nField == 2) 
  {
    temp1<-attributes(e1)
    temp2<-attributes(e2)
    time1<-time(e1)
    #Trick, bring it on the same timebase... start as to be > 0 due to R-bug
    e1<-ts(e1,start=100)
    e2<-ts(e2,start=100)
    
    if (abs(sum(c(time(e1))-c(time(e2))))>TLIMIT) 	warning("Operator applied on two timebases, new timebase = first timebase")
    
    return(pTs(NextMethod(),time1,temp1$lat,temp1$lon,newname,c(paste("A:",temp1$history),paste("B:",temp2$history),newname)))
  }
  else
  {
    e<-NextMethod(.Generic)
    if (!is.null(attr(e,'nclass'))) class(e)<-attr(e,'nclass')
    return(addhistory(e,newname))
  }
}








#if p1 and p2 [x,y] or only p [x] are given: give back value
#if only p1 is given   [x,] :  pField: field at one time
#if only p2 is given   [,x]: pTs   : univariate time series

#Å¸ber die Dimensionen arbeiten...
#wenn [1:n,] dann length (dim=NULL) oder ncol (wenn dim != null) gleich ncol(x)
#damit kann dies von [1:n] unterschieden werden 

#FÅ‰lle zum testen: [a,b] [a:a1,b]  [a,b:b1] [a:a1,b:b1] 
# [a,] [a:a1,] [,b] [,b:b1] [a] [a:a1] [date] .... das ganze auf normalem Feld.. und Feld mit einem


"[.pField"<-function(x,p1,p2,...)
{
  result<-NextMethod("[")	
  temp<-attributes(x)
  
  if (!missing(p2)) l2<-length(p2)
  
  if (!missing(p1))
  {
    l1<-length(p1)
    if (missing(p2)) 
    {     
      #check in [a] oder [a,]
      if ((is.null(dim(result)) & (length(result)==ncol(x))) | !is.null(dim(result))) 
      {  
        #[a,]  give back field at times a:a1
        hist<-paste("[",p1[1],":",p1[l1],", ]",sep="")
        if (length(p1) > 1) result<-t(result)
        result<-pField(result,time(x)[p1],temp$lat,temp$lon,temp$name,temp$history,date=FALSE)
      }
      else
      {
        #[a]  directly give back value
        hist<-paste("[",p1[1],":",p1[l1],"]",sep="")						
      }
      
      
    }
    else 
    {
      #[a,b] directly give back value or area
      hist<-paste("[",p1[1],":",p1[l1],",",p2[1],":",p2[l2],"]",sep="")
    }
  }
  else if (!missing(p2)) 
  {
    #[,b]  give back time series at position b
    hist<-paste("[,",p2[1],":",p2[l2],"]",sep="")
    
    pos<-c1t2(p2,length(temp$lon))
    result<-pTs(result,time(x),temp$lat[pos$lat],temp$lon[pos$lon],temp$name,temp$history,date=FALSE)
    
  }
  
  if (!is.null(attr(result,"history"))) return(addhistory(result,hist))  else return(result)
  
  
}


#[a,] time area a out of timeseries
#[a] time area a out of timeseries
#[a,b] point
#[,b] time series with index b

"[.pTs"<-function(x,p1,p2,...)
{
  result<-NextMethod("[")	
  temp<-attributes(x)
  
  if (!missing(p2)) l2<-length(p2)
  
  if (!missing(p1))
  {
    l1<-length(p1)
    if (missing(p2)) 
    {     
      #check in [a] oder [a,]
      
      #[a,]  give back points at times a:a1
      hist<-paste("[",p1[1],":",p1[l1],",]",sep="")
      result<-pTs(result,time(x)[p1],temp$lat,temp$lon,temp$name,temp$history,date=FALSE)
      
      
      
      
    }
    else 
    {
      #[a,b] directly give back value or area
      hist<-paste("[",p1[1],":",p1[l1],",",p2[1],":",p2[l2],"]",sep="")
    }
  }
  else if (!missing(p2)) 
  {
    #[,b]  give back time series at position b
    hist<-paste("[,",p2[1],":",p2[l2],"]",sep="")
    
    pos<-c1t2(p2,length(temp$lat))
    result<-pTs(result,time(x),temp$lat[p2],temp$lon[p2],temp$name[p2],date=FALSE)
    
  }
  
  if (!is.null(attr(result,"history"))) return(addhistory(result,hist))  else return(result)
  
  
  
}


binAvg<-function(x,y,N=NULL,breaks=pretty(x,N),bFill=FALSE)
  
  #Averages y into bins according to the positon of a in the breaks
  #Either give N=number of breaks, or N+1 breaks
  #Breaks are defined as x>breaks[i], and x<=breaks[i+1]
  #see also stats.bin in library(fields) which returns a more comprehensive statistics
  #if fill=T, fill empty bins using linear interpolation from the neighbours to the center of the bin
  #Returns the breaks, centers, the averaged values and nobs, the number of observations averages
{
  NBIN <- length(breaks) - 1
  centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
  avg<-rep(NA,length(breaks)-1)
  nobs<-rep(NA,length(breaks)-1)
  for (i in 1:(length(breaks)-1)) {
    selection<-y[which((x>breaks[i])&(x<=breaks[i+1]))]
    avg[i]<-mean(na.omit(selection))
    nobs[i]<-sum(!is.na(selection))
  }
  
  
  if ((sum(is.na(avg))>0)&(bFill))
  {
    yInt<-approx(x,y,centers)$y
    missing<-is.na(avg)
    avg[missing]<-yInt[missing]
  }
  return(list(breaks=breaks,centers=centers,avg=avg,nobs=nobs))
}


