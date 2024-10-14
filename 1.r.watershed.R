# === general ==================================================================

## 1.watershed.R 
## this script executes the GRASS command r.watershed on the carved DEM. The output of this script is used by 2.r.basin.R
#
# input script: 0.carve.DEM.R
#
#
# There are some lines of code commented out. This code has filepaths for 
# running on my linux server, and should be ignored :)
#
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.0.2 (2020-06-22)"
#
#

# === table of contents ==================================================================

# 1. load libraries and setup GRASS
#
# 2. import CDEM
# 
# 3. r. watershed
#
# 4. save r.watershed output

# ===  1. load libraries and setup GRASS ==================================================================

library(rgrass7)

initGRASS(gisBase = "C:/OSGEO4~1/apps/grass/grass78", home = tempdir(), mapset = "PERMANENT") # for OSGEO4W install of grass v7.8
# initGRASS(gisBase = "/usr/lib/grass78", home = tempdir(), mapset = "PERMANENT") # for linux install of grass v7.8

# print mapset name
execGRASS("g.mapset", flags = c("p","verbose"))

# project mapset to EPSG:4269
execGRASS("g.proj", flags = "c", epsg = 4269)

# update the mapset's projection
execGRASS("g.region", flags = "d")

# print the mapset's CRS to confirm it has been changed
execGRASS("g.region", flags = "p")


# ===  2. import data ==================================================================

# import digital elevation model (DEM)
execGRASS("r.import",input = "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/BC TRIM 25m DEM/CDEM_EPSG4269_3m_deep_streams_from_sum.tif", # /home/tk/GIS_data/CDEM_EPSG4269_3m_deep_streams_from_sum.tif
          output = "van_island_DEM")
# execGRASS("r.import",input = "/home/tk/GIS_data/CDEM_EPSG4269_3m_deep_streams_from_sum.tif", output = "van_island_DEM")

# set the GRASS mapset to the extent of the DEM
execGRASS("g.region",raster = "van_island_DEM")


# ===  3. r.watershed ==================================================================

execGRASS("r.watershed", elevation = "van_island_DEM",accumulation="carved_flow_accumulation", drainage = "carved_flow_direction", flags = "overwrite")

# show information about the map that was created
execGRASS("r.info",map = "carved_flow_accumulation")
execGRASS("r.info",map = "carved_flow_direction")



# ===  4. save r.watershed output ==================================================================

execGRASS("r.out.gdal", input = "carved_flow_accumulation", output = "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/BC TRIM 25m DEM/accum_3m_deep_streams_from_sum.tif",
          format = "GTiff")
# execGRASS("r.out.gdal", input = "carved_flow_accumulation", output = "/home/tk/GIS_data/accum_3m_deep_streams_from_sum.tif", format = "GTiff")

execGRASS("r.out.gdal", input = "carved_flow_direction", output = "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/BC TRIM 25m DEM/flow_direction_3m_deep_streams_from_sum.tif",
          format = "GTiff")
# execGRASS("r.out.gdal", input = "carved_flow_direction", output = "/home/tk/GIS_data/flow_direction_3m_deep_streams_from_sum.tif", format = "GTiff")
