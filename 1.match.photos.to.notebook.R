# === general ==================================================================

## 1.match.photos.to.notebook.R

# join the photo IDs to the notebook trap level data (depth, substrate)

# === table of contents ==================================================================

# 1. load libraries and import data
#
# 2. match depth and distance from shore
#
# 3. one hot encode substrate categories

# ===  1. load libraries and import data ==================================================================
library(tidyverse)

notebook <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2020/Summer 2020 datasheets/2020 tkosciuch fieldwork notebook data/2020 tkosciuch fieldwork notebook.csv", fileEncoding = 'UTF-8-BOM') %>% 
  distinct()
photos <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/size_corrected_lateral_lengths_mm.csv")

# ===  1. load libraries and import data ==================================================================

photos$lake_id <- substring(photos$id,1,nchar(photos$id)-6)

for(i in photos)


