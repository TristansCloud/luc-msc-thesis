

## This work is now in 4.2BPLs.PCA.elastic.net.R







# === general ==================================================================

## 1.lake.order.R 
## Find the lake order using the maximum stream order
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.1.2 (2021-11-01)"

# === table of contents ==================================================================

# This script uses old data, there are now two datasets available that have both
# kosciuch and stuart data:
#   lake means with parasite data: "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_kosciuch_stuart_lake_mean_parasitenotsc.csv"
#   individuals with parasite data: "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_kosciuch_stuart_morpho.csv"

# 1. import and format landmarks and setup libraries
#
# 2. PCA and visualize results

# ===  1. import and format landmarks and setup libraries ==================================================================
library(dplyr)
library(ggfortify)
library(geomorph)
library(abind)

# For trait lengths when I eventually add parasite data in:
# traits <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_kosciuch_stuart_lake_mean_parasitenotsc.csv")
# Until I add parasite data, use this dataset:
# traits <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_mean_size_corrected_trait_lengths.csv")


# read in landmark data and match landmarks between yoel and kosciuch
match_landmarks <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/match_yoel_kosciuch_landmarks.csv") %>% 
  filter(!is.na(kosciuch_landmark))
yoel_landmarks <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Yoels 2017 project/OG files/05.geomo_nom.csv") %>% 
  select(c("fishID.univ", "filename", "Csize", contains("Coord"))) %>% 
  filter(substr(fishID.univ,4,4) == "L") %>% # keep only lake fish
  mutate(lake = substr(fishID.univ,1,3)) %>% 
  select(!contains("gpaCoord.lm18"))
# this is where Yoel got landmark locations: https://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2009.00665.x
datapath <- "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/OJJ file backups"
## Lateral image data import. lateral_missing are samples that were found to be missing from the lateral dataset so need to be added in.
## there is more info in the README in the datapath.
lateral <- read.csv(file.path(datapath,"Mar1lateral_without_missing_fish.csv"), fileEncoding = 'UTF-8-BOM')
lateral_missing <- read.csv(file.path(datapath,"Dec10lateralmissingfish.csv"), fileEncoding = 'UTF-8-BOM')
tk_landmarks <- bind_rows(lateral,lateral_missing) %>% 
  select(contains(c("id","scale.15.mm",match_landmarks$kosciuch_landmark))) %>% 
  mutate(lake = substr(id,1,11)) %>% 
  filter(!lake %in% c("MH00385SALM", "PY00098SALM"))


# take the mean values of each lake
yoel_lake_lm <- yoel_landmarks %>% 
  select(-c(fishID.univ, filename)) %>% 
  group_by(lake) %>% 
  summarise_all(mean)
  
tk_lake_lm <- tk_landmarks %>% 
  select(-c(id)) %>% 
  group_by(lake) %>% 
  summarise_all(mean)

# prep to combine landmarks
temp_yoel_lm <- c(paste(match_landmarks$yoel_landmark,"x", sep = "."),
                  paste(match_landmarks$yoel_landmark,"y", sep = "."))
lm_order <- order(temp_yoel_lm)
temp_yoel_lm <- temp_yoel_lm[lm_order]
temp_tk_lm <- c(paste(match_landmarks$kosciuch_landmark,"x", sep = "."),
                paste(match_landmarks$kosciuch_landmark,"y", sep = "."))
temp_tk_lm <- temp_tk_lm[lm_order]

yoel_lake_lm <- yoel_lake_lm[,c("lake", "Csize", temp_yoel_lm)]
tk_lake_lm <- tk_lake_lm[,c("lake", "scale.15.mm", temp_tk_lm)]


# create geomorph matrices of coordinates and combine coordinates
  
yoel_coords <- yoel_lake_lm %>%
  select(-c("Csize","lake"))
yoel_coords <- arrayspecs(as.matrix(yoel_coords), p = 18, k = 2)
tk_coords <- tk_lake_lm %>% 
  select(-c("lake", "scale.15.mm"))
tk_coords <- arrayspecs(as.matrix(tk_coords), p = 18, k = 2) %>%
  rotate.coords(., type = "flipY") # my coordinates are dorsal side down so I flip them
plotAllSpecimens(yoel_coords, label = T)
plotAllSpecimens(gpagen(tk_coords)$coords, label = T)

coords <- gpagen(abind(yoel_coords, tk_coords))

plotAllSpecimens(coords$coords, label = T)

# ===  2. PCA and visualize results ==================================================================

### Body shape PCA
coord_pca <- gm.prcomp(coords$coords)
pca_plot <- plot(coord_pca)
t <- picknplot.shape(pca_plot) # , method = "points", mag = 2, links = coord_pca$shapes


### PCA on trait lengths

PCA <- prcomp(traits[,-which(names(traits) %in% c("id","lake_id","WTRBDPLD","dataset","watershed"))], scale. = TRUE, center = TRUE)
summary(PCA)
# plot first two principal components and outline by lake, saved as a png
png(file="C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/PCA.png",
    width=600, height=350)
autoplot(PCA, data = traits)
dev.off()

PCA$rotation

# write loadings and proportion of variance to txt file to be imported to a table in word
write.table(round(summary(PCA)$importance,3),"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/PCA_summary.txt")
write.table(round(PCA$rotation,3),"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/PCA_loadings.txt")
