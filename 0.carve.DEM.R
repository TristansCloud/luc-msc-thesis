# === general ==================================================================

## this script carves the CDEM by 3 meters using the FWA stream network
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
# 2. create a raster from the stream network with streams having cell value of -3
# 
# 3. subtract stream raster from the DEM



# ===  1. load libraries and setup GRASS ==================================================================
library(raster)
library(rgdal)

streams <- readOGR("/home/tk/GIS_data/FWSTRMNTWR_line","FWSTRMNTWR_line")
DEM <- raster("/home/tk/GIS_data/CDEM Vancouver Island EPSG4269.tif")


# ===  2. create a raster from the stream network with streams having cell value of -3 ==================================================================

stream.raster <- raster()
extent(stream.raster) <- extent(DEM)
res(stream.raster) <- res(DEM)
stream.rastert <- rasterize(streams, stream.raster, -3) # runs slow, but less slow than r.carve!
writeRaster(stream.rastert, filename = "/home/tk/GIS_data/stream_raster_3m_deep.tif", format="GTiff")


# ===  3. subtract stream raster from the DEM ==================================================================

carved_DEM <- sum(DEM, stream.rastert, na.rm = TRUE) # runs very quickly, but sets the min and max cell value to +- infinity.. The raster looks fine otherwise, all cells look normal to me and the subtraction worked
writeRaster(carved_DEM, filename = "/home/tk/GIS_data/CDEM_EPSG4269_3m_deep_streams_from_sum.tif", format="GTiff")
