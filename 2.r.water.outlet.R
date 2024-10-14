# === general ==================================================================

## 2.r.stream.basin.R 
## this script executes the GRASS command r.stream.basin with the Canadian 25m resolution CDEM of
## Vancouver Island. It first finds the outlet coordinates of the lake, then uses these outlets
## as the drainage point for each lake. 
# authors: T. Kosciuch
#
# input script: 1.r.watershed.R
# 
# R version
R.version.string
# "R version 4.0.2 (2020-06-22)"
#

# === table of contents ==================================================================

# 1. Load libraries and setup GRASS
#
# 2. Import data: BC Freshwater atlas (FWA) lake shapefiles, flow accumulation raster, and Canadian digital elevation model (CDEM)
#
# 3. r.stream.basin. This step is a for loop with multiple substeps:
# 
#  3.1 Find the lake drainage outlet's coordinates for each lake using the maximum value of the flow accumulation raster contained within the lake's polygon
#
#  3.2 Use the drainage outlet and CDEM as inputs for r.water.outlet, then turn the raster drainage basin into a polygon
#


# === 1. Load libraries and setup GRASS ==================================================================

library(rgrass7)
library(raster)
library(sf)
# library(rgdal)

initGRASS(gisBase = "C:/OSGEO4~1/apps/grass/grass78", home = tempdir(), mapset = "PERMANENT", override=FALSE)

# print mapset name
execGRASS("g.mapset", flags = c("p","verbose"))

# project mapset to EPSG:4269
execGRASS("g.proj", flags = "c", epsg = 4269)

# update the mapset's projection
execGRASS("g.region", flags = "d")

# print the mapset's CRS to confirm it has been changed
execGRASS("g.region", flags = "p")

# === 2.Import data ==================================================================

lakes <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/Stuart 2017 and Kosciuch 2020 sample lakes/Stuart 2017 and Kosciuch 2020 sample lakes")
# lakes <- st_read("/home/tk/GIS_data/Stuart 2017 and Kosciuch 2020 sample lakes")


accum <- raster("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/BC TRIM 25m DEM/flow_accumulation.tif") # flow_accumulation.tif used the 3m carved DEM
# accum <- raster("/home/tk/GIS_data/accum_3m_deep_streams_from_sum.tif")

# adjust the projection to that of the accumulation raster
lakes <- st_transform(lakes, projection(accum))

# === 3. r.basin ==================================================================

## 3.1 lake drainage outlet's coordinates

# these are the columns of the outlet dataframe. WTRBDKGRPC, WTRBDPLD, and WTRBDK together should be enough to id any polygon individually or other polygons part of the same lake
x <- vector(mode = "numeric", length = dim(lakes)[1]) ; y <- vector(mode = "numeric", length = dim(lakes)[1]) ;
WTRBDPLD <- vector(mode = "numeric", length = dim(lakes)[1]) ; WTRBDK <- vector(mode = "numeric", length = dim(lakes)[1]) # WTRBDKGRPC <- vector(mode = "double", length = dim(lakes)[1]) ; # this one wasn't working

outlets <- data.frame(x,y,WTRBDPLD, WTRBDK)

for(i in 1:dim(lakes)[1]){ 
  
  tempmask <- mask(crop(accum,lakes[i,]),lakes[i,])
  
  tempcord.i <- which.max(tempmask)
  
  tempcord.i <- xyFromCell(tempmask,tempcord.i)
  
  outlets[i,c("x","y")] <- tempcord.i
  
  outlets[i,c("WTRBDPLD","WTRBDK")] <- with(lakes[i,], c(WTRBDPLD,WTRBDK))
  
  # plot(tempmask)+points(outlets[i,])
  
}

write.csv(outlets, "/home/tk/GIS_data/outlets_3m_deep_streams_from_sum.csv")
outlets <- read.csv("/home/tk/GIS_data/outlets_3m_deep_streams_from_sum.csv")

##  3.2 Use the drainage outlet and CDEM as inputs for r.basin. I merged 
##  the individual polygons into one layer in QGIS and fixed issues where
##  polygons self intersected using st_make_valid() from the sf package.

#import accumulation raster - uncomment if you want to overwrite the file
# execGRASS("r.import",input = "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/BC TRIM 25m DEM/flow_direction_3m_deep_streams_from_sum.tif", 
#           output = "carved_flow_direction", flags = "overwrite")
# execGRASS("r.import",input = "/home/tk/GIS_data/flow_direction_3m_deep_streams_from_sum.tif", output = "carved_flow_direction", flags = "overwrite")


lakes.WTRBDPLD <- vector(length = dim(lakes)[1], mode = "character")

for(i in 1:dim(lakes)[1]){ 
  
  coords <- c(outlets[i,"x"], outlets[i,"y"])
  
  basin.name <- as.character(lakes[i,]$WTRBDPLD)
  
  # create drainage basin raster.
  execGRASS("r.water.outlet", input = "carved_flow_direction", output = basin.name, coordinates = coords, flags = "overwrite")
  
  basin.vector.name <- paste("vector", basin.name, sep = "_")
  
  # Turn the raster drainage basin into a polygon and add it to a vector layer
  execGRASS("r.to.vect", input = basin.name, output = basin.vector.name, type = "area",  flags = c("overwrite","s"))
  
  # update the vectors attribute table with the lake's WTRBDPLD value
  execGRASS("v.db.update", map = basin.vector.name, column = "value", value = basin.name)
  
  # write the polygon to a shapefile
  basin.filepath <- paste("/home/tk/GIS_data/thesis_drainage_basins",basin.vector.name, sep = "/")
  
  execGRASS("v.out.ogr", input = basin.vector.name, output = basin.filepath, format = "ESRI_Shapefile")
  
  lakes.WTRBDPLD[i] <- basin.name
}


