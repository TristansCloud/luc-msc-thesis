# === general ==================================================================

## 1.setup.GRASS.R 
## this script creates a GRASS mapset which can be used for the rest of the analyses
#
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.0.2 (2020-06-22)"
#
#

# ==============================================================================

library(tidyverse)
library(rgrass7)

initGRASS(gisBase = "C:/OSGEO4~1/apps/grass/grass78", home = tempdir(), mapset = "PERMANENT")

# print mapset name
execGRASS("g.mapset", flags = c("p","verbose"))

# project mapset to EPSG:4269
execGRASS("g.proj", flags = "c", epsg = 4269)

# update the mapset's projection
execGRASS("g.region", flags = "d")

# print the mapset's CRS to confirm it has been changed
execGRASS("g.region", flags = "p")

use_sf()
