# === general ==================================================================

## 1.watershed.R 
## this script find the distance as the fish swims between lakes. Currently doesn't work.
#
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.0.2 (2020-06-22)"
#
#

# === table of contents ==================================================================

# 1. import data
#
# i may need to do v.net operation=connect first (http://manpages.ubuntu.com/manpages/bionic/man1/v.net.path.1grass.html)
#
# 2. v.overlay: calculates intersections
#
# 3. v.net.distance
#
#
#
# I think a good way to do this is to find the intersection of the FWA lake polygon and the stream network, then do the path analysis between all intersections.
# I would have to do some for loops to find which distance is the shortest.
# I really want to test if the non connected streams will be a problem
#

# === 1. import data ==================================================================

library(tidyverse)
library(rgrass7)
library(sf)

lakes <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/Stuart 2017 and Kosciuch 2020 sample lakes/Stuart 2017 and Kosciuch 2020 sample lakes.shp")

streams <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/NAD83 CRS FWA files/FWA_VanIsld_streams/FWA_STREAM_NETWORKS_SP/FWSTRMNTWR_line.shp")

# === 2. find intersection points ==================================================================

# filter intersect to keep only points
ls_intersect <- st_intersection(lakes, streams) %>% 
  filter(st_is(.,"POINT"))

# === 3. setup GRASS ==================================================================

# import stream network
execGRASS("v.import",input = "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/NAD83 CRS FWA files/FWA_VanIsld_streams/FWA_STREAM_NETWORKS_SP/FWSTRMNTWR_line.shp",
          output = "FWA_streams")

# set the GRASS mapset to the extent of the stream network
execGRASS("g.region",vector = "FWA_streams")

# execGRASS("v.import",input = "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/2020 Sampling/2020 sample lakes/2020 sample lakes.shp",
#           output = "sampled_lakes")

# === 2. v.net.distance ==================================================================

# import nodes to GRASS
writeVECT()

# v.net.distance must use points as the "from" input, but can have boundaries as the "to" input. Lake polygon = boundary?
execGRASS("v.net.distance",vector = "FWA_streams")

###### Developing ############

# two lakes that are linked by streams: WTRBDPLD:710017417 (Little mud lake) & WTRBDPLD:710017268 (Roberts lake)

# Two that are not linked: WTRBDPLD:710021978 (Durrance lake in VICT) & WTRBDPLD:710017268 (Roberts lake)

lilmud_points <- ls_intersect %>% filter(WTRBDPLD == 710017417)
writeVECT(lilmud_points,"lilmud_points")

Roberts_polygon <- lakes[lakes$WTRBDPLD == 710017268,] %>% filter(!is.na(WTRBDPLD))
writeVECT(Roberts_polygon, "Roberts_polygon",v.in.ogr_flags=c("overwrite"))

# list vector layers
execGRASS("g.list",type="vector")

# doesn't work, need to do v.net first!
# execGRASS("v.net.distance",input = "FWA_streams", from_layer = "lilmud_points", 
#           to_layer = "Roberts_polygon", to_type = "boundary", output = "lilmud_roberts")

writeVECT(ls_intersect, "ls_intersect")

execGRASS("v.net", input = "FWA_streams", points = "ls_intersect", operation = "connect", output = "ls_network", threshold = 0.1)

# this never finished running but still created an output. The output was corrupted.
execGRASS("v.net.allpairs",input = "ls_network", output = "ls_net_allpairs")

ls_net_allpairs <- readVECT("ls_net_allpairs")


