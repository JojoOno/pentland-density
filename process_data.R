



######################################################################################
#                                                                                    #
###############################  Spatial Modelling ################################### 
#                      
### Packages ###
require(sp)
require(sf)
require(rgdal)
require(lme4)
require(effects)
require(sjPlot)
library(mgcv)
require(geepack)
require(splines)
require(MuMIn)
require(raster)
require(ggplot2)
require(tidyverse)
require(ggspatial)
require(raster)
######################################################################################

load("data/GPS-COMBINED-FLOW.Rd") 

filter_lon_lat <- function(df, lon_min, lon_max, lat_min, lat_max) {
  # Subset the data frame to include only rows where lon and lat fall within the specified range
  filtered_df <- df[df$lonUTM >= lon_min & df$lonUTM <= lon_max & df$latUTM >= lat_min & df$latUTM <= lat_max, ]
  return(filtered_df)
}


remove <- c(65254, 65231,64303, "pv24-150-11", "pv24-151-11", "pv24-598-11", "pv24-580-11", 65446, 65243, 64312)

data <- gps_flow_comb %>%
  filter(hoflag==0) %>%
  filter(!(tag %in% remove)) %>%
  mutate(lonUTM=lonUTM*1000, latUTM=latUTM*1000) %>%
  filter_lon_lat(456533.58761129, 516953.619065722, 6481742.69728062, 6523097.60610944 )


######################################################################################

### Load Data ###

benthos <- read.asciigrid("Data/mapping/Benthic/D4_2018.asc/D4_2018.asc") 
proj4string(benthos) <- CRS("+init=epsg:4326")

# Don't reproject before finding points --> have to use original projection as converting to UTM creates an irregular 
# grid and therefore coerces to points which will not work with the 'over' function in sp
# benthos <- spTransform(benthos, CRS("+proj=utm +ellps=WGS84 +datum=WGS84 +zone=30 +north +units=km"))

## Crop the benthic data - It's fucking massive and takes an age to extract covariates for GPS data in its original format
coords <- cbind(data$Lon, data$Lat)
sp.points <- SpatialPoints(coords, proj4string=CRS("+init=epsg:4326"))
extent(sp.points)
benthos <- raster::crop(benthos, extent(extent(sp.points)))


######################################################################################
############################## Generate pseudo-absences ##############################
resolution <- 500 # grid-cell resolution in km

##### Establish land or sea   
xmin <- floor(456533.58761129)
xmax <- ceiling(516953.619065722)
ymin <- floor(6481742.69728062)
ymax <- ceiling(6523097.60610944)


rastermap <- raster(xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, resolution=resolution, crs="+proj=utm +ellps=WGS84 +datum=WGS84 +zone=30 +north +units=km")   # Using extended raster so usage is not underestimated.
res(rastermap) <- resolution # Resolution in km

rastermap[] <- 1:(rastermap@nrows * rastermap@ncols) #  Setting the values, #rows x #cols
plot(rastermap)
#setwd("C://users/jo26/Dropbox/PhD/Pentland Firth Data/MRSea.Data")
writeRaster(rastermap,"Data/GridIDs.asc", format='ascii', NAflag=-9999, overwrite=TRUE)
#summary(rastermap)

# Importing UK map (shape file)
UKmap <- readOGR("ScotLatLon", dsn="mapping/scotland/maximum-resolution") 
proj4string(UKmap) <- CRS("+init=epsg:4326")
UKmap <- spTransform(UKmap, CRS("+proj=utm +ellps=WGS84 +datum=WGS84 +zone=30 +north +units=m"))

plot(rastermap) # Check it with a plot
plot(UKmap, add=T) # Add on the coastline map

rastermap[] <- 1 # change all rastermap values to 1
landseaMap <- rasterize(UKmap, rastermap, mask=T, background=0) # Overlay & cookie cut the UK map to the rastermap
landseaMap[ is.na(landseaMap) ] <- 0    # Setting NAs to 0
plot(landseaMap)
writeRaster(landseaMap, "data/LandSeaMap.asc", format='ascii', NAflag=-9999, overwrite=TRUE)


##### Generate pseudo-absences at sea, accounting for "availability"

seals <- unique(data$tag)  
zeroN<-2 #  number of pseudoabsences per presence
i=1
for(i in 1:length(seals)) {
  
  d <- subset(data, data$tag==seals[i]) # subset by seal
  ss<-nrow(d)
  e <- extent(min(d$lonUTM), max(d$lonUTM), min(d$latUTM), max(d$latUTM)) # create extent that mimcs the limits of the seal's trips throughout the tagging period
  refmap <- crop(landseaMap, e)	 # creates a separate raster just for the extent of the seal's dispersal
  
  cellVal <- which(t(as.matrix(refmap) )== 0) # get All at sea cell
  samples <- sample(cellVal, zeroN*ss, replace=TRUE) # randomly sample from at sea possibilities
  coords <- xyFromCell(refmap, samples) # get coordinates for centre point of each sampled cell
  
  zx <- coords[,1]
  zy <- coords[,2]
  zt<-rep(d$datetime, zeroN)
  zs<-rep(d$tag, zeroN)
  
  data1<-data.frame("used"=rep(1,ss), "time"=as.POSIXct(d$datetime, "%Y/%m/%d %H:%M:%S", tz="GMT"),"id"=d$tag, "Lon"=d$lonUTM, "Lat"=d$latUTM)
  data0<-data.frame("used"=rep(0,zeroN*ss), "time"=as.POSIXct(zt, "%Y/%m/%d %H:%M:%S", tz="GMT"),"id"=zs, "Lon"=zx, "Lat"=zy)
  datax<-rbind(data1, data0)
  
  if(i==1){mod_data <- datax} else
  {mod_data <- rbind(mod_data, datax)}
  
  print(i)
}


#zeroN<-6 # Number of zeros corresponding to each telemetry point
#ss<-nrow(atsea) # Pooled sample size
#zx<-runif(zeroN*ss, min(atsea$Lon), max(atsea$Lon))
#zy<-runif(zeroN*ss, min(atsea$Lat), max(atsea$Lat))
#zt<-rep(atsea[,7], zeroN)
#zs<-rep(atsea[,4], zeroN)
#data1<-data.frame("used"=rep(1,ss), "time"=as.POSIXct(atsea$datetime, "%Y/%m/%d %H:%M:%S", tz="GMT"),"id"=atsea$Tag, "Lon"=atsea$Lon, "Lat"=atsea$Lat)
#data0<-data.frame("used"=rep(0,zeroN*ss), "time"=as.POSIXct(zt, "%Y/%m/%d %H:%M:%S", tz="GMT"),"id"=zs, "Lon"=zx, "Lat"=zy)
#data<-rbind(data1, data0)
#data<-cbind(data, xr=round(data$Lon), yr=round(data$Lat))
#st<-cbind(data$xr, data$yr)
#n<-nrow(data)


##### Generate depth and mean tidal flow values for each location, including pseudo-absences

coordinates(mod_data) <- ~Lon+Lat 
proj4string(mod_data) <- CRS("+proj=utm +ellps=WGS84 +datum=WGS84 +zone=30 +north +units=m")
data2 <- spTransform(mod_data, CRS("+init=epsg:4326"))

benthos <- raster(benthos)
data$seadepth <- extract(benthos,data2)
data <- data.frame(data)
data[is.na(data)]<-0
data$seadepth <- ifelse(data$seadepth > 0, 0, data$seadepth)


## define impact as after the turbines were installed

mod_data$impact <- ifelse(mod_data$time<"2016-10-31 00:00:00", 0, 1)  


#separate files as pre 2014 requires different licence
Tides16 <- read.csv("data/2016-2018HWLW.csv") 
Tides11 <- read.csv("data/2010-2012HWLW.csv")
Tides <- rbind(Tides11, Tides16)
Tides$Tide <- "NA"                                                  
Tides$Tide[1] <- ifelse(Tides$height[1] < Tides$height[2], "LW", "HW")  
for (i in 2:NROW(Tides)){                                           # loop for remaining rows
  Tides$Tide[i] <- ifelse(Tides$height[i] < Tides$height[i-1], "LW", "HW")  # if height lower than row above = "LW", else "HW"
}

nearest.time <- function(x, y, key = "datetime"){y[which.min(abs(difftime(x, y[,key]))),key]}


# subset tide data to HW only
HW2016.2018 <- Tides[Tides$Tide == "HW",]

# create posixct DateTime columns
HW2016.2018$datetime <- as.POSIXct(paste(HW2016.2018$date, HW2016.2018$time), format="%d/%m/%Y %H:%M", origin="1970-01-01", tz = "GMT")

# new column of time of nearest high water
mod_data$NearestHW <- as.POSIXct(mapply(function(x) nearest.time(x, HW2016.2018), mod_data$time), origin="1970-01-01", tz = "GMT")

# column of time around nearest high water
mod_data$TimeAroundHW <- difftime(mod_data$time, mod_data$NearestHW, units = "hours")
mod_data$TimeAroundHW <- as.numeric(mod_data$TimeAroundHW)



save(mod_data, file=paste("data/PresAbsDat.Rd"))
write.csv(mod_data, file="data/PresAbsDat.csv")



