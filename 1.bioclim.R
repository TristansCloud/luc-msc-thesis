# === general ==================================================================

## 1.bioclim.R
## get max and min temps and precipitation by month for each lake using PRISM
## climatology data from https://data.pacificclimate.org/portal/bc_prism/map/.
#
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.0.2 (2020-06-22)"
#
#
# === table of contents ==================================================================

# 1. load libraries and import data
#
# 2. Min temp (Celsius) per lake
#
# 3. Max temp (Celsius) per lake
#
# 4. Lake precipitation (mm) per month
#
# 5. Drainage basin precipitation (m^3) per month

# ===  1. load libraries and import data ==================================================================

library(sf)
library(raster)
library(ncdf4)
library(tidyverse)

lakes <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/Stuart 2017 and Kosciuch 2020 sample lakes/Stuart 2017 and Kosciuch 2020 sample lakes")
basins <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/drainage basins/merged_drainage_basins"); basins <- st_make_valid(basins)

# not sure if the filename should end with ".nc.nc" but it does on my machine
precipition <- brick("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/PRISM bioclim/precipitation monthly/pr_mClimMean_PRISM_historical_19810101-20101231.nc.nc")
maxtemp <- brick("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/PRISM bioclim/max temp monthly/tasmax_mClimMean_PRISM_historical_19810101-20101231.nc.nc")
mintemp <- brick("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/PRISM bioclim/min temp monthly/tasmin_mClimMean_PRISM_historical_19810101-20101231.nc.nc")

# ===  2. Min temp per lake ==================================================================

lake_mintemp <- raster::extract(mintemp,lakes,fun=mean,small=TRUE)

#lakes$WTRBDPLD is a unique identifier of each lake polygon
lake_mintemp <- cbind(lakes$WTRBDPLD, lake_mintemp)
colnames(lake_mintemp) <- c("WTRBDPLD",paste(month.name,"min temp"))

write.csv(lake_mintemp,"C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/lake_mintemp.csv")

# ===  3. Max temp per lake ==================================================================

lake_maxtemp <- raster::extract(maxtemp,lakes,fun=mean,small=TRUE)

#lakes$WTRBDPLD is a unique identifier of each lake polygon
lake_maxtemp <- cbind(lakes$WTRBDPLD, lake_maxtemp)
colnames(lake_maxtemp) <- c("WTRBDPLD",paste(month.name,"max temp"))

write.csv(lake_maxtemp,"C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/lake_maxtemp.csv")

# ===  4. Lake precipitation (mm) per month ==================================================================

# average lake precipitation in mm
lake_precipitation_mean <- raster::extract(precipition,lakes,fun=mean,small=TRUE)
# add unique lake identifier WTRBDPLD
lake_precipitation_mean <- cbind(lakes$WTRBDPLD,lake_precipitation_mean)
# update column names
colnames(lake_precipitation_mean) <- c("WTRBDPLD",paste(month.name,"lake precipitation (mm)"))

write.csv(lake_precipitation_mean, "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/lake_precipitation_mean.csv")

# ===  5. Drainage basin precipitation (m^3) per month ==================================================================

# average basin precipitation in mm
wtrshd_precipitation_mean <- raster::extract(precipition,basins,fun=mean,small=TRUE)
# mm precipitation * 0.001 meters/mm * basin area (m^2) = volume of precipitation (m^3)
wtrshd_precipitation_vol <- wtrshd_precipitation_mean*0.001*st_area(basins)           

# basins$value is WTRBDPLD, I should change this in the shapefile
wtrshd_precipitation_vol <- cbind(basins$value,wtrshd_precipitation_vol) 
colnames(wtrshd_precipitation_vol) <- c("WTRBDPLD",paste(month.name,"basin precipitation (m^3)"))
# need this group_by as the basins has more polygons than lakes. This happened because some polygons self intersected near the edges and split into 2+ polgyons.
# I sum the volume of precipitation from each polygon of the lakes watershed.
wtrshd_precipitation_vol <- data.frame(wtrshd_precipitation_vol) %>% group_by(WTRBDPLD) %>% summarise(across(everything(),list(sum)))
colnames(wtrshd_precipitation_vol) <- c("WTRBDPLD",paste(month.name,"basin precipitation (m^3)")) # colnames get changed so I run it twice, very efficient :)


write.csv(wtrshd_precipitation_vol, "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/wtrshd_precipitation_vol.csv")


