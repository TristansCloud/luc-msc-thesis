# === general ==================================================================

## 1.lake.order.R 
## Find the lake order using the maximum stream order
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.1.2 (2021-11-01)"
#

# === table of contents ==================================================================

# 1. import data and setup libraries
#
# 2. intersect FWA lake and stream polygons to find the maximum stream order within the lake.

# ===  1. load libraries and import data ==================================================================

library(sf)
library(tidyverse)

lakes <- st_read("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/GIS data/Stuart 2017 and Kosciuch 2020 sample lakes/Stuart 2017 and Kosciuch 2020 sample lakes")
streams <- st_read("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/GIS data/NAD83 CRS FWA files/FWA_VanIsld_streams/FWA_STREAM_NETWORKS_SP/FWSTRMNTWR_line.shp")

# ===  2. intersect FWA lake and stream polygons to find the maximum stream order within the lake. ==================================================================

ints_streams <- st_intersection(streams,lakes) # takes ~10 minutes to run
lake_order <- data.frame(tapply(ints_streams$STRMRDR, ints_streams$WTRBDPLD,max))

# add rownames (WTRBDPLD) as a column, then remove the rownames
lake_order <- setNames(cbind(rownames(lake_order), lake_order, row.names = NULL), 
                           c("WTRBDPLD", "lake_order"))

# write.csv(lake_order,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/lake_order.csv",row.names = FALSE)



