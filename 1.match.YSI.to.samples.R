# === general ==================================================================

## 1.match.YSI.to.samples.R

# Join the 2021 YSI samples to 2021 notebook data

# to do: add misty, theimer, and beaver lake YSI samples using YSI GPS most likely

# === table of contents ==================================================================

# 1. load libraries and import data
#
# 2. match YSI to projects
#
# 3. join thesis YSI to FWA lake shapefile

# ===  1. load libraries and import data ==================================================================
library(tidyverse)
library(lubridate)
library(sf)

YSI <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2021/Summer 2021 notebook data/All projects YSI measurements/formatted YSI all projects.csv")
cline <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2021/Summer 2021 notebook data/All projects YSI measurements/YSI measurement times for clines.csv", fileEncoding = 'UTF-8-BOM')
thesis <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2021/Summer 2021 notebook data/All projects YSI measurements/2021 kosciuch thesis water and ysi.csv")
GPS <- st_read("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2021/Summer 2021 GIS Data/Summer 2021 GPS points/Working copy of 2021 GPS points.gpx") %>% 
  filter(nchar(name) <= 3) %>%  # drop points that are not named 001, 002, etc.
  mutate(name = as.integer(name))
lakes <- st_read("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/GIS data/Stuart 2017 and Kosciuch 2020 sample lakes/Stuart 2017 and Kosciuch 2020 sample lakes")


# ===  2. match YSI to projects ==================================================================

## cline
ysi_cline <- data.frame()

for(i in unique(cline$date.m.d.yyyy)){
  temp_YSI <- YSI[YSI$DATE == i,]
  temp_cline <- cline[cline$date.m.d.yyyy == i,]
  
  x <- temp_cline$YSI.time
  y <- temp_YSI$TIME
  
  x1 <- hms::as_hms(format(strptime(x, "%I:%M:%S %p"), "%H:%M:%S"))
  y1 <- hms::as_hms(format(strptime(y, "%I:%M:%S %p"), "%H:%M:%S"))
  
  for(time in 1:length(x1)){
    diff = min(abs(difftime(x1[time],y1)))
    # if the difference is less than 12 seconds it's a likely match
    if(diff < 12){
      temp_output <- cbind(temp_cline[time,],temp_YSI[diff == abs(difftime(x1[time],y1)),])
      ysi_cline <- rbind(ysi_cline,temp_output)
    }
  }
}

# check if the number of samples is the same after matching to YSI
dim(ysi_cline)[1] == dim(cline)[1]
# FALSE because one of the Marius samples had an ambiguous YSI time. The ysi points have GPS so they should be usable still
cline$YSI.time[!cline$YSI.time %in% ysi_cline$YSI.time]
# "1:48:45 PM"

# write.csv(ysi_cline,"C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2021/Summer 2021 notebook data/All projects YSI measurements/analysis ready cline ysi data.csv")

## thesis

ysi_thesis <- data.frame()

for(i in mdy(unique(thesis$date.m.d.yyyy))){
  temp_YSI <- YSI[mdy(YSI$DATE) == i,]
  temp_thesis <- thesis[mdy(thesis$date.m.d.yyyy) == i,]
  
  x <- temp_thesis$YSI.time
  y <- temp_YSI$TIME
  
  x1 <- hms::as_hms(format(strptime(x, "%H:%M:%S"), "%H:%M:%S"))
  y1 <- hms::as_hms(format(strptime(y, "%I:%M:%S %p"), "%H:%M:%S"))
  
  for(time in 1:length(x1)){
    diff = min(abs(difftime(x1[time],y1)))
    # if the difference is less than 12 seconds it's a likely match
    if(diff < 12){
      temp_output <- cbind(temp_thesis[time,],temp_YSI[diff == abs(difftime(x1[time],y1)),])
      ysi_thesis <- rbind(ysi_thesis,temp_output)
    }
  }
}

# has 3 samples with NA data, apparently the YSI did not record info for those samples.


# ===  3. join thesis YSI to FWA lake shapefile ==================================================================

GPS <- st_transform(GPS,st_crs(lakes))

GPS_lake_i <- st_nearest_feature(GPS, lakes)

GPS_lake <- cbind(GPS,lakes[GPS_lake_i,])

ysi_thesis <- left_join(ysi_thesis,GPS_lake[,c("name","WTRBDPLD")],by = c("GPS"="name"))
thesis <- left_join(thesis,GPS_lake[,c("name","WTRBDPLD")],by = c("GPS"="name"))

# write.csv(ysi_thesis[,!colnames(ysi_thesis) %in% "geometry"],"C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2021/Summer 2021 notebook data/All projects YSI measurements/analysis ready thesis ysi data.csv",row.names = FALSE)
# write.csv(thesis[,!colnames(thesis) %in% "geometry"],"C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2021/Summer 2021 notebook data/All projects YSI measurements/analysis ready thesis water sample data.csv",row.names = FALSE)


