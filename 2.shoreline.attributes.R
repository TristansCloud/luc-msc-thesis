# === general ==================================================================

## 2.shoreline.attributes.R
## calculate the shoreline development and the difference in shoreline elevation
## and the lake surface elevation
#
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.0.2 (2020-06-22)"
#
#

# === table of contents ==================================================================

# 1. load libraries and read in data
#
# 2. calculate shoreline development
#
# 3. shoreline - lake surface elevation
#    - transform the lakes into ESRI:102002 to use meters 
#    - calculate the buffer
#    - save the difference between the buffer and the lake polygon
#    - transform back to EPSG:4269
#    - find lake elevation
#    - find difference in lake and shore elevation and write to csv

# ===  1. load libraries and read in data ==================================================================
library(raster)
library(sf)
library(tidyverse)

lakes <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/Stuart 2017 and Kosciuch 2020 sample lakes/Stuart 2017 and Kosciuch 2020 sample lakes")
DEM <- raster("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/BC TRIM 25m DEM/CDEM Vancouver Island EPSG4269.tif")

# ===  2. calculate shoreline development ==================================================================
# shoreline development is the lake perimeter / perimeter of a perfect circle
# with the same surface area as the lake. In the FWA lake polygons, AREA_SQM is
# the surface area in meters squared, and FEAT_LEN is the perimeter in meters.
# To convert from area to circumference, I used Circumference = 2*sqrt(pi*Area)
shore_dev <- lakes$FEAT_LEN/(sqrt(lakes$AREA_SQM*pi)*2)
shore_dev <- data.frame(lakes$WTRBDPLD,shore_dev); colnames(shore_dev)[1] <- "WTRBDPLD"
write.csv(shore_dev,"C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/shoreline_development.csv")

# ===  3. shoreline - lake surface elevation ==================================================================
# Following a method modified from Heathcote et al 2015:
# buffer (km) = 0.075km + 0.5*(Area(km^2)/pi)^0.5 assuming km and km^2

# I first transform the lakes to ESRI:102002, Canada_Lambert_Conformal_Conic
# so I can correctly use meters as the unit of distance
lcc_lakes <- st_transform(lakes, "ESRI:102002")

# create a vector with the lake surface area in km^2
AREA_KM2 <- lcc_lakes$AREA_HA/100
# Find the buffer distance in meters
BUFFER_M <- (0.075+0.5*sqrt(AREA_KM2/pi))*1000
# create the buffered polygons
buff_lakes <- st_buffer(lcc_lakes,BUFFER_M)
# create an empty simple feature
buffer_polygons <- st_sf(st_sfc(crs = "ESRI:102002"))
# remove the lake surface from the poylgon
for(i in 1:dim(buff_lakes)[1]){
  buffer_polygons <- rbind(buffer_polygons,st_difference(buff_lakes[i,],lcc_lakes[i,]))
}
# transform back to ESPG:4269
buffer_polygons <- st_transform(buffer_polygons,4269)

# find the average elevation within the buffer
shoreline_elevation <- data.frame(buffer_polygons$WTRBDPLD,raster::extract(DEM,buffer_polygons,fun=mean))
colnames(shoreline_elevation) <- c("WTRBDPLD","mean_shoreline_elevation(m)")

# lake elevation
lake_centers <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/Stuart 2017 and Kosciuch 2020 sample lakes/moved_lake_centroids_for_v_net")
lake_elevation <- raster::extract(DEM,lake_centers)
# subtract lake elevation from shoreline elevation (I checked that the lake elevations are in the same order as the shoreline elevations)
shoreline_elevation$`diff_shore_lake_elevation(m)` <- shoreline_elevation$`mean_shoreline_elevation(m)`-lake_elevation
# add lake elevations to the data.frame
shoreline_elevation$lake_elevation <- lake_elevation

write.csv(shoreline_elevation,"C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/lake_and_shoreline_elevations.csv")

