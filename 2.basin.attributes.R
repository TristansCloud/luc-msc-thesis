# === general ==================================================================

## 2.basin.attributes.R
## surface area (m^2) of drainage basins
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
# 2. Drainage basin surface area


# ===  1. load libraries and import data ==================================================================
library(sf)
library(tidyverse)

basins <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/drainage basins/merged_drainage_basins"); basins <- st_make_valid(basins)

# ===  2. Drainage basin surface area ==================================================================

basin_area <- st_area(basins); basin_area <- cbind(basins$value,basin_area)
basin_area <- data.frame(basin_area) %>% 
  group_by(V1) %>% summarise(basin_area = sum(basin_area))
colnames(basin_area) <- c("WTRBDPLD", "basin area (m^2)")

write.csv(basin_area, "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/basin_area.csv")

