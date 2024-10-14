# === general ==================================================================

## 2.icp.metals.formatting.R 
## Create analysis read ICP data
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.1.2 (2021-11-01)"
#

# === table of contents ==================================================================

# 1. load libraries and import data
#
# 2. delete rows with no data and drop unnecessary columns
#
# 3. merge lake ID info with ICP data
#
# 4. fix issue tags
#
# 5. write data

# ===  1. load libraries and import data ==================================================================

library(dplyr)
library(tidyr)
library(sf)
library(ggfortify)

ICP <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/ICP MS metals/ICP_metals_csv_pre_formatting.csv", fileEncoding = 'UTF-8-BOM')
water_sample_metadata <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2021/Summer 2021 notebook data/All projects YSI measurements/analysis ready thesis water sample data.csv") %>% 
  filter(lake.surface.or.bottom == "surface")

# ===  2. delete rows with no data and drop unnecessary columns ==================================================================

ICP <- ICP[,c("Parameter.","Client.Sample.ID.","Results..")]

# change to 0 if below detection limit
ICP$Results..[grepl("<",ICP$Results..)] <- 0 

# make the colnames nicer.
colnames(ICP) <- c("Parameter","Client.Sample.ID.","Results")

# Drop rows that are labels, not results
ICP <- ICP[!ICP$Parameter == "Dissolved Metals Filtration Location",]

# change results to numeric
ICP$Results <- as.numeric(ICP$Results)

# ===  3. merge lake ID info with ICP data ==================================================================

ICP_w_GPS <- ICP[grepl("GPS",ICP$Client.Sample.ID.),]
ICP_n_GPS <- ICP[!grepl("GPS",ICP$Client.Sample.ID.),]

# check the GPS is always at the end of the client ID. It is. Further, it is always the last 3 characters
unique(ICP_w_GPS$Client.Sample.ID.)

ICP_w_GPS$GPS <- as.integer(substring(ICP_w_GPS$Client.Sample.ID., nchar(ICP_w_GPS$Client.Sample.ID)-2))

# match water sample GPS to ICP GPS
matched_ICP_w_GPS <- left_join(ICP_w_GPS,water_sample_metadata[,c("GPS","WTRBDPLD","site")], by = ("GPS"))

# sites that don't have a GPS in the ICP sitename
unique(ICP_n_GPS$Client.Sample.ID.)
# "BEAVER LAKE" "MISTY LAKE"  "THEIMER" 

# Initialize the column
ICP_n_GPS$WTRBDPLD <- NA
ICP_n_GPS[ICP_n_GPS$Client.Sample.ID. == "THEIMER",]$WTRBDPLD <- 710014253
ICP_n_GPS[ICP_n_GPS$Client.Sample.ID. == "MISTY LAKE",]$WTRBDPLD <- 710014259
ICP_n_GPS[ICP_n_GPS$Client.Sample.ID. == "BEAVER LAKE",]$WTRBDPLD <- 710002615

 
# THEIMER == 710014253
# MISTY == 710014259
# BEAVER == two options for beaver lake, 710021962 in victoria, 710002615 in the
#    northern part of the island. Since I didn't catch any fish at victoria beaver
#    lake I didn't take a water sample there. 

# Bind ICP results together
ICP_out <- bind_rows(matched_ICP_w_GPS[,c("Parameter","Results","WTRBDPLD","site","Client.Sample.ID.")],ICP_n_GPS[,c("Parameter","Results","WTRBDPLD","Client.Sample.ID.")])

ICP_wide <- pivot_wider(ICP_out,names_from = "Parameter",values_from = "Results")

# ===  4. fix issue tags  ==================================================================

# There are three samples that don't have the right GPS and are ambigous in which lake they belong to.
# There are an additional 2 samples (00403SALM GPS 194, 00403SALM GPS 209) that don't have the correct GPS
# for the lake name listed but by elimination can be assigned to the correct lake.

# the three ambigous samples. Matches with ICP_wide$Client.Sample.ID.
issue_sites <- c("ORMOND LAKE GPS 164", "SWAN/SOREN LAKE GPS 196","CECIL GPS 165")

## The issue sites can't be properly assigned to a lake. The following code 
## drops them from the samples. The commented out code looks to see if the 
## issue sites cluster close to any other lakes, particularly the closest sample
ICP_out <- ICP_out %>% filter(!Client.Sample.ID. %in% issue_sites)
ICP_wide <- ICP_wide %>% filter(!Client.Sample.ID. %in% issue_sites)

# # based on the lake name (Ormond, swan/soren, and cecil) is there a lake with a known water sample immediately adjacent?
# # closest_lake matches with ICP_wide$site
# closest_lake <- c("00316SALM","moore","roberts")
# 
# # using PCA plots to see if the points are close to where they should be
# metal <- ICP_wide[,4:dim(ICP_wide)[2]]
# # drop columns with only 0 values
# metal <- metal[,colSums(metal) > 0]
# 
# metal_pca <- prcomp(metal, scale. = TRUE, center = TRUE)
# summary(metal_pca)
# metal_pca$rotation
# 
# id <- rep("normal", dim(metal)[1])
# id[ICP_wide$site == "00316SALM"] <- "set1"
# id[ICP_wide$site == "moore"] <- "set2"
# id[ICP_wide$site == "roberts"] <- "set3"
# id[ICP_wide$Client.Sample.ID. == "ORMOND LAKE GPS 164"] <- "set1"
# id[ICP_wide$Client.Sample.ID. == "SWAN/SOREN LAKE GPS 196"] <- "set2"
# id[ICP_wide$Client.Sample.ID. == "CECIL GPS 165"] <- "set3"
# 
# 
# t <- cbind(metal,id)
# autoplot(metal_pca, data = t, colour = 'id', legend = TRUE)
# 

# ===  5. write data ==================================================================

# write.csv(ICP_out,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/ICP_long_format.csv",row.names = FALSE)
# 
# write.csv(ICP_wide,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/ICP_wide_format.csv",row.names = FALSE)



