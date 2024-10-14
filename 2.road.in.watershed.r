# === general ==================================================================

## 2.road.in.watershed.r
## get the length of roads in each lake's watershed
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
# 2. find the intersection of the watershed and roads, then calculate road length
#    - roads are separated into paved and unpaved
#
# 3. join road length tables
#
# 4. testing that overlapping watershed polygons will function correctly
#    - they do function correctly

# ===  1. load libraries and import data ==================================================================

library(sf)
library(tidyverse)

basins <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/drainage basins/merged_drainage_basins")
roads <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/BCGW_7113060B_1634341996709_9936/DRA_DGTL_ROAD_ATLAS_MPAR_SP")
basins <- st_transform(basins,st_crs(roads))

# ===  2. find the intersection of the watershed and roads, then calculate road length ==================================================================

# separate paved and unpaved roads
paved <- roads %>% filter(RDSURFACE == "paved")
unpaved <- roads %>% filter(RDSURFACE != "paved")

## paved
ints_paved <- st_intersection(paved,basins)
paved_length <- data.frame(tapply(st_length(ints_paved), ints_paved$value,sum))
# add rownames (WTRBDPLD) as a column, then remove the rownames
paved_length <- setNames(cbind(rownames(paved_length), paved_length, row.names = NULL), 
                           c("WTRBDPLD", "paved_road_length(m)"))

## unpaved
ints_unpaved <- st_intersection(unpaved,basins)
unpaved_length <- data.frame(tapply(st_length(ints_unpaved), ints_unpaved$value,sum))
# add rownames (WTRBDPLD) as a column, then remove the rownames
unpaved_length <- setNames(cbind(rownames(unpaved_length), unpaved_length, row.names = NULL), 
         c("WTRBDPLD", "unpaved_road_length(m)"))

# ===  3. join road length tables ==================================================================

## join tables
road_length <- full_join(unpaved_length,paved_length)
# change paved length NA to 0
road_length$`paved_road_length(m)`[is.na(road_length$`paved_road_length(m)`)] <- 0
# write.csv(road_length, "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/drainage basins/road_length.csv")

# ===  4. testing that overlapping watershed polygons will function correctly ==================================================================

### Proof that overlapping polygons don't cause a problem with finding road lengths
# two watersheds with overlapping boundaries
plot(basins[basins$value %in% c(710003484, 710003485),1])
ints2 <- st_intersection(roads,basins[basins$value %in% c(710003484, 710003485),])
tapply(st_length(ints2), ints2$value,sum)
# 710003484 710003485 
# 769576   1417115

# checking if this method gave the same answer as the individual polygons
ints3 <- st_intersection(roads,basins[basins$value %in% c(710003485),])
ints4 <- st_intersection(roads,basins[basins$value %in% c(710003484),])
tapply(st_length(ints3), ints3$value,sum) # 1417115 
tapply(st_length(ints4), ints4$value,sum) # 769576 
# gives the same answer, I don't have to use a for loop!
