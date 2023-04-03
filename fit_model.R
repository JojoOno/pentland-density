
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Author: Joe Onoufriou
## Date: 2023
## Purpose: Produced for Ewan Edwards at Xodus for density estimates within the MeyGen array for collision risk modelling
## Description: Code to model the change in distribution of north coast seals throughout the tidal cycle and predict the densities across 500 x 500 metre grid-cells throughout the range. 
## Notes: Can take up to 3 days to run the final model.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


################
### Packages ###
require(sp)
require(maptools)
require(rgdal)
require(raster)
require(splines)
require(MRSea)
require(randtests)
require(ggplot2)
require(tidyverse)
require(rgeos)
require(sf)
require(ggspatial)
require(sandwich)
require(geepack)
require(car)
require(wesanderson)
require(extrafont)
require(RColorBrewer)
require(mapr)
require(gridExtra)

#######################################

coast <- st_read("mapping/scotland/maximum-resolution/GSSHS_British_Isles poly.shp") %>%
  st_set_crs(4326) %>%
  st_transform(32630) %>%
  as(Class="Spatial")

pent <- st_read("mapping/scotland/maximum-resolution/ScotLatLon.shp") %>%
  st_set_crs(4326) %>%
  st_transform(32630) 

pent_lr <- st_read("mapping/scotland/maximum-resolution/GSSHS_British_Isles poly.shp") %>%
  st_set_crs(4326) %>%
  st_transform(32630) 

meygen <- st_read("mapping/wave-and-tidal-proposed-sites/TCE_Lease_Tide_20160919.shp") %>%
  st_set_crs(4326) %>%
  st_transform(32630) %>%
  filter(Name_Prop=="Inner Sound") %>%
  st_cast("POLYGON")


########################################
load("data/PresAbsDat.Rd")

mod_data <- mod_data %>%
  mutate(x.pos=Lon, y.pos=Lat) %>%
  mutate(response=used)  %>%
  mutate(TideState4=ifelse(between(TimeAroundHW, -1.5,1.5), "HW", TimeAroundHW)) %>%
  mutate(TideState4=ifelse(between(TimeAroundHW, 1.5,5.5), "EBB", TideState4)) %>%
  mutate(TideState4=ifelse(between(TimeAroundHW, 5.5,7)|between(TimeAroundHW, -7,-5.5), "LW", TideState4)) %>%
  mutate(TideState4=ifelse(between(TimeAroundHW, -5.5,-1.5), "FLOOD", TideState4)) 

########################################################
### subset for first run otherwise GEE takes forever ###

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#mod_dat <- sample_n(mod_data, 10000) # sample when first running just to check everything is in order
mod_dat <- mod_data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

########################################################

#######################################
### Make a grid of candidate knots ###

latnew<-seq(min(round(mod_data$y.pos))-500,max(round(mod_data$y.pos))+500,by=2000)
lonnew<-seq(min(round(mod_data$x.pos))-500,max(round(mod_data$x.pos))+500,by=2000)
knotgrid<-expand.grid(lonnew,latnew)
names(knotgrid)<-c("X","Y")
grd<-SpatialPoints(knotgrid, 
                   proj4string=CRS("+proj=utm +ellps=WGS84 +datum=WGS84 +zone=30 +north +units=m"))
knotgrid$land<-over(grd, coast)
knotgrid[,1][!is.na(knotgrid$land[1])]<-NA
knotgrid[,2][!is.na(knotgrid$land[1])]<-NA

knotgrid<-knotgrid[,1:2]
rownames(knotgrid)<-1:dim(knotgrid)[1]
knotgrid <- na.omit(knotgrid)

##############################
### Raster for predictions ###

raster <- raster("Data/GRIDIDs.asc")
proj4string(raster) <- "+proj=utm +ellps=WGS84 +datum=WGS84 +zone=30 +north +units=m"
r.pts <- rasterToPoints(raster, spatial=TRUE)


################## Modelling ###################
#####
# make prediction grid
#####

predictDataFLOOD <- data.frame("x.pos"=r.pts@coords[,1], "y.pos"=r.pts@coords[,2],
                               "TideState4"=rep("FLOOD"))

predictDataEBB <- data.frame("x.pos"=r.pts@coords[,1], "y.pos"=r.pts@coords[,2],
                             "TideState4"=rep("EBB"))

predictDataHW <- data.frame("x.pos"=r.pts@coords[,1], "y.pos"=r.pts@coords[,2],
                            "TideState4"=rep("HW"))

predictDataLW <- data.frame("x.pos"=r.pts@coords[,1], "y.pos"=r.pts@coords[,2],
                            "TideState4"=rep("LW"))

preddata <- rbind(predictDataFLOOD, predictDataEBB, predictDataHW, predictDataLW)

dists <- makeDists(cbind(preddata$x.pos, preddata$y.pos),
                   na.omit(knotgrid), knotmat = FALSE)$dataDist

################# set initial model without the spline terms and build SALSA ####################
initialModel <- glm(response~1, family="binomial", data=mod_dat)  

distMats <- makeDists(cbind(mod_dat$x.pos, mod_dat$y.pos), knotgrid)
r_seq <- getRadiiChoices(numberofradii=8, distMats$dataDist, basis='gaussian')

salsa2dlist <- list(fitnessMeasure = "AIC", knotgrid = knotgrid,
                    startKnots = c(6), minKnots = c(4), maxKnots = c(25), r_seq = r_seq,
                    gap = c(1), interactionTerm = c("as.factor(TideState4)"))

salsa2dOutput_k6 <- runSALSA2D(initialModel, salsa2dlist,
                               distMats$dataDist, distMats$knotDist)
#Warning messages:  1: glm.fit: fitted probabilities numerically 0 or 1 occurred ---- Does not matter at this stage as we are going to re-fit using GEE


bestModel <- salsa2dOutput_k6$bestModel
splineParams <- bestModel$splineParams


radiusIndices <- splineParams[[1]]$radiusIndices
aR <- splineParams[[1]]$knotPos 
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii

mod_dat <- mod_dat[order(mod_dat$id),] 
mod_dat$blockid <- mod_dat$id

MRSeaGam.imp <- make.gamMRSea(bestModel, panelid=mod_dat$id) ## re-fit using GEE - i.e. fit robust sandwich based estaimtes of variance - EWAN - not actually an issue here as we're not calculating CIs for your stuff given how difficult it is to interpret confidence intervals for density estimate across the range for each grid-cell
summary(MRSeaGam.imp)
anova(MRSeaGam.imp)

chosenknots<-knotgrid[bestModel$splineParams[[1]]$knotPos,]

ggplot() +
  ggspatial::annotation_spatial(data=pent) +
  geom_point(data=data.frame(knotgrid), aes(x=X, y=Y)) +
  geom_point(aes(x=X, y=Y, colour="red", size=2), data=data.frame(chosenknots), alpha=4/5, show.legend = FALSE)+
  theme_bw() + xlab('Easting (Km)') + ylab('Northing (Km)')


######
#create prediction grid with land removed 
######

pred_grd<-st_as_sf(preddata, coords=c("x.pos", "y.pos"), crs=32630) %>%
  mutate(x.pos = unlist(map(.$geometry,1)),y.pos = unlist(map(.$geometry,2)))

pred_grd$land <- sapply(st_intersects(pred_grd,pent_lr), function(z) if (length(z)==0) NA_integer_ else z[1])

preddata_sea <- filter(pred_grd, is.na(land))
# Note that you can skip the intermediate step of converting the filtered data back to a data frame and keep it as an sf object if you prefer
dists <- makeDists(cbind(preddata_sea$x.pos, preddata_sea$y.pos),
                   na.omit(knotgrid), knotmat = FALSE)$dataDist

predslink <- predict.gamMRSea(newdata=preddata_sea, g2k=dists, object=MRSeaGam.imp, type = "link")

predsexp <- exp(predslink)


result_grid <- preddata_sea[1:6821,] %>%
  select(-TideState4) %>%
  mutate(flood=predsexp[preddata_sea$TideState4=="FLOOD"], ebb= predsexp[preddata_sea$TideState4=="EBB"],
          HW=predsexp[preddata_sea$TideState4=="HW"], LW=predsexp[preddata_sea$TideState4=="LW"])%>%
    mutate(predsHW.perc = HW/sum(HW)*100, predsLW.perc = LW/sum(LW)*100,
         predsFLOOD.perc = flood/sum(flood)*100, predsEBB.perc=ebb/sum(ebb)*100)



dims<-getPlotdimensions(x.pos=preddata$x.pos, preddata$y.pos,
                        segmentWidth=500, segmentLength=500)

par(mfrow=c(2,2))


quilt.plot(preddata$x.pos[preddata$TideState4=="FLOOD"],
           preddata$y.pos[preddata$TideState4=="FLOOD"],
           predsexp[preddata$TideState4=="FLOOD"], asp=1, nrow=dims[1], ncol=dims[2], zlim=c(0, 3),
           main="FLOOD")

quilt.plot(predictDataHW$x.pos,
           predictDataHW$y.pos,
           predsHW, asp=1, nrow=dims[1], ncol=dims[2], zlim=c(0, 3),
           main="HW")

quilt.plot(predictDataLW$x.pos,
           predictDataLW$y.pos,
           predsHW, asp=1, nrow=dims[1], ncol=dims[2], zlim=c(0, 3),
           main="LW")

quilt.plot(predictDataEBB$x.pos,
           predictDataEBB$y.pos,
           predsEBB, asp=1, nrow=dims[1], ncol=dims[2], zlim=c(0, 3),
           main="EBB")


################## Plots ##################


font_import()
loadfonts(device = "win")
pal <- wes_palette("Zissou1", 11, type = "continuous")
pal2 <- palette(rev(brewer.pal(n = 11, name = "Spectral")))

lease  <- st_read("Data/mapping/wave-and-tidal-proposed-sites/TCE_Lease_Tide_20160919.shp") %>%
  st_set_crs(4326) %>%
  st_transform(32630)



pflood <- ggplot() +
  theme_classic() +
  geom_tile(data=result_grid, aes(x.pos, y.pos, fill = ifelse(predsFLOOD.perc>0.5, 0.5, predsFLOOD.perc))) +
  scale_fill_gradientn(colours=pal2, guide = guide_colourbar(title.position = "top", barwidth = 20, draw.ulim = FALSE, draw.llim = FALSE))+
  #scale_fill_viridis_c(guide = guide_colourbar(title.position = "top", barwidth = 20, draw.ulim = FALSE, draw.llim = FALSE))+
  #geom_point(data=predictDATA.grid.grid[marker==2,], aes(x=x.pos, y=y.pos), shape=3)+
  #geom_point(data=predictDATA.grid.grid[marker==1,], aes(x=x.pos, y=y.pos), shape="-", size=5) +
  annotation_spatial(data=pent_lr) +
  #geom_sf(data=lease, fill="darkorchid4", alpha=0.5)+
  theme(panel.grid.major = element_line(colour = "transparent"))+
  labs(fill = "Estimated Percentage")+
  theme(legend.direction = "horizontal", legend.position = "top", legend.box = "vertical",
        text=element_text(size=16, family="serif"))+
  xlab("Longitude") + ylab("Latitude")+
  ggtitle("Flooding Tide")

pflood

pebb <- ggplot() +
  theme_classic() +
  geom_tile(data=result_grid, aes(x.pos, y.pos, fill = ifelse(predsEBB.perc>0.5, 0.5, predsEBB.perc))) +
  scale_fill_gradientn(colours=pal2,  guide = guide_colourbar(title.position = "top", barwidth = 20, draw.ulim = FALSE, draw.llim = FALSE))+
  #scale_fill_viridis_c(limits = c(0, 0.5), guide = guide_colourbar(title.position = "top", barwidth = 20, draw.ulim = FALSE, draw.llim = FALSE))+
  
  #geom_point(data=predictDATA.grid[marker==2,], aes(x=x.pos, y=y.pos), shape=3)+
  #geom_point(data=predictDATA.grid[marker==1,], aes(x=x.pos, y=y.pos), shape="-", size=5) +
  annotation_spatial(data=pent_lr) +
  # geom_sf(data=lease, fill="darkorchid4", alpha=0.5)+
  theme(panel.grid.major = element_line(colour = "transparent"))+
  labs(fill = "Estimated Percentage")+
  theme(legend.direction = "horizontal", legend.position = "top", legend.box = "vertical",
        text=element_text(size=16, family="serif"))+
  xlab("Longitude") + ylab("Latitude")+
  ggtitle("Ebbing Tide")


plw <- ggplot() +
  theme_classic() +
  geom_tile(data=result_grid, aes(x.pos, y.pos, fill = ifelse(predsLW.perc>0.5, 0.5, predsLW.perc))) +
   scale_fill_gradientn(colours=pal2,  guide = guide_colourbar(title.position = "top", barwidth = 20, draw.ulim = FALSE, draw.llim = FALSE))+
  #scale_fill_viridis_c(limits = c(0, 0.5), guide = guide_colourbar(title.position = "top", barwidth = 20, draw.ulim = FALSE, draw.llim = FALSE))+
  
  #geom_point(data=predictDATA.grid[marker==2,], aes(x=x.pos, y=y.pos), shape=3)+
  #geom_point(data=predictDATA.grid[marker==1,], aes(x=x.pos, y=y.pos), shape="-", size=5) +
  annotation_spatial(data=pent_lr) +
  # geom_sf(data=lease, fill="darkorchid4", alpha=0.5)+
  theme(panel.grid.major = element_line(colour = "transparent"))+
  labs(fill = "Estimated Percentage")+
  theme(legend.direction = "horizontal", legend.position = "top", legend.box = "vertical",
        text=element_text(size=16, family="serif"))+
  xlab("Longitude") + ylab("Latitude")+
  ggtitle("Low Tide")

phw <- ggplot() +
  theme_classic() +
  geom_tile(data=result_grid, aes(x.pos, y.pos, fill = ifelse(predsHW.perc>0.5, 0.5, predsHW.perc))) +
  scale_fill_gradientn(colours=pal2, guide = guide_colourbar(title.position = "top", barwidth = 20, draw.ulim = FALSE, draw.llim = FALSE))+
  #scale_fill_viridis_c(limits = c(0, 0.5), guide = guide_colourbar(title.position = "top", barwidth = 20, draw.ulim = FALSE, draw.llim = FALSE))+
  
  #geom_point(data=predictDATA.grid[marker==2,], aes(x=x.pos, y=y.pos), shape=3)+
  #geom_point(data=predictDATA.grid[marker==1,], aes(x=x.pos, y=y.pos), shape="-", size=5) +
  annotation_spatial(data=pent_lr) +
  # geom_sf(data=lease, fill="darkorchid4", alpha=0.5)+
  theme(panel.grid.major = element_line(colour = "transparent"))+
  labs(fill = "Estimated Percentage")+
  theme(legend.direction = "horizontal", legend.position = "top", legend.box = "vertical",
        text=element_text(size=16, family="serif"))+
  xlab("Longitude") + ylab("Latitude")+
  ggtitle("High Tide")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(phw)
grid.arrange(legend,
             phw + theme(legend.position="none"),
             pebb+ theme(legend.position="none"),
             plw + theme(legend.position="none"),
             pflood + theme(legend.position="none"), 
             
             layout_matrix = rbind(c(1,1), 
                                   c(2,3), 
                                   c(4,5)),
             heights=c(0.1,0.7,0.7))



chosenknots_tide<-knotgrid[salsa2dOutput_k6$bestModel$splineParams[[1]]$knotPos,]


ggplot(data.frame(knotgrid)) +
  geom_point(aes(x=X, y=Y)) +
  geom_point(aes(x=X, y=Y, size=2), colour=wes_palette("Zissou1")[1], data=data.frame(chosenknots_tide), alpha=4/5, show.legend = FALSE) + 
  theme_bw() + xlab('Longitude') + ylab('Latitude') +
  annotation_spatial(data=pent_lr)


##########
## Calculate total percentage throughout MeyGen lease site
###########

# Read in the "predictions" data frame and convert it to an sf object
predictions_sf <- st_as_sf(result_grid[,c(1,2,9:12)], coords = c("x.pos", "y.pos"), crs = 32630)

# Initialize an empty list to store the sums of densities for each column
sum_densities_list <- list()
avg_densities_list <- list()
tide_names <- c("High Water", "Low Water", "Flood", "Ebb")

# Iterate over each tidal state
for (i in 3:6) {
  # Perform a spatial join to select only the densities that are within "meygen"
  densities_within_meygen <- predictions_sf %>% 
    st_join(meygen, join = st_within) %>% 
    filter(!is.na(.[[7]])) %>% # remove rows with NA values in the column of interest
    st_drop_geometry() %>% # drop the geometry column after the join
    select(i) %>% 
    pull()
  
  # Calculate the sum of the densities for this column and store it in the list
  sum_densities_list[[tide_names[i-2]]] <- sum(densities_within_meygen)
  avg_densities_list[[tide_names[i-2]]] <- mean(densities_within_meygen)
  }

print(sum_densities_list) # give you the total percentage of the north coast population (not including Orkney) that are within the MeyGen site during different states of tide

#$`High Water`
#[1] 0.8113636

#$`Low Water`
#[1] 0.4068655

#$Flood
#[1] 0.8908913

#$Ebb
#[1] 0.8395395

print(avg_densities_list) # give you the mean percentage of the north coast population (not including Orkney) that are within the MeyGen site during different states of tide

#$`High Water`
#[1] 0.05795454

#$`Low Water`
#[1] 0.02906182

#$Flood
#[1] 0.0636351

#$Ebb
#[1] 0.0599671

####
write.csv(result_grid, file="results/prediction-grid.csv")
