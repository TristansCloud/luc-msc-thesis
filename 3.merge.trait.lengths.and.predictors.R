# === general ==================================================================

## 3.merge.trait.lengths.and.predictors.R 
## Join FWA lake polygon information to the list of lakes I actually sampled.
## Merge 
#
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.1.2 (2021-11-01)"
#

# === table of contents ==================================================================

# 1. load libraries and import data
#
# 2. merge enviro predictors into one dataset
#
# 3. merge stuart 2017 morpho and tristan's collection morpho
#    - also account for missing features
#    - write individual specimens data and lake means as separate files
#
# 4. merge notebook and FWA lake data

# ===  1. load libraries and import data ==================================================================
library(tidyr)
library(dplyr)
library(sf)
library(geomorph)
library(readxl)
library(stringr)

# there are not more than 37 Kosciuch lakes with fish and 16 Stuart ones
notebook <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2020/Summer 2020 datasheets/2020 tkosciuch fieldwork notebook data/2020 tkosciuch fieldwork notebook.csv", fileEncoding = 'UTF-8-BOM') %>% 
  distinct()
kosciuch_lake_joining_info <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/notebook_lakes.csv") %>% 
  select(lake_name,tag_prefix,WTRBDKGRPC,WTRBDPLD) %>% 
  mutate(tag_prefix = toupper(tag_prefix),
         lake_name = toupper(lake_name)) %>% 
  mutate(lake_name = case_when(lake_name == "THETIS LAKE EAST (LOWER)" ~ "THETIS EAST",
                               lake_name == "THETIS LAKE WEST (UPPER)" ~ "THETIS WEST",
                               .default = lake_name)) %>% 
  distinct() %>% 
  mutate(watershed = paste(tag_prefix, WTRBDKGRPC, sep = ""))

kosciuch_notebook_specimens_WTRBDPLD <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/kosciuch_notebook_specimens.csv")


# Add these to create a dataset of non allometric scaled traits
lateral <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_lateral.csv")
ventral <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_ventral.csv")

## Joe lake from yoel's data doesn't have a GNSNM1 ID

lakes <- st_read("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/GIS data/Stuart 2017 and Kosciuch 2020 sample lakes/Stuart 2017 and Kosciuch 2020 sample lakes")
lakes_df <- lakes %>% st_drop_geometry()
kosciuch_lake_counts <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2020/Summer 2020 datasheets/Tristan count of fish by lake.csv")

sc_lateral <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_lateral_sc.csv")
sc_ventral <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_ventral_sc.csv")
kosciuch_gill_raker_n <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Gill Rakers/Gill Raker Counts + Measures & Gravid Scoring.xlsx", sheet = "gill.raker.counts") %>% 
  select(id, right.side.gillraker.number) %>% 
  mutate(right.side.gillraker.number = as.numeric(right.side.gillraker.number)) # NAs are ok
stuart_gill_raker_n <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Yoels 2017 project/OG files/01.univMorph.csv") %>% filter(habitat == "lake") %>% 
  select(fishID.univ, Right.side.gill.raker.number.insitu) %>% 
  rename(id = fishID.univ,
         right.side.gillraker.number = Right.side.gill.raker.number.insitu)
gill_raker_n <- bind_rows(kosciuch_gill_raker_n, stuart_gill_raker_n)

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

yoel_watershed <- unique(sc_lateral$watershed)[38:53] %>% str_replace("\\.", " ") %>% capwords() %>% paste("Lake") %>% cbind(unique(sc_lateral$watershed)[38:53]) %>% data.frame()
names(yoel_watershed) <- c("GNSNM1", "watershed")


#Joe watersheds doesnt make it and beaver is in there twice. Yoel's beaver lake is up north, not the VICT one. Joe lake doesnt have a GNSNM1 but it doesn't matter
yoel_lakes_df <- lakes_df %>% inner_join(yoel_watershed) %>% 
  filter(WTRBDKGRPC != "00247VICT") %>% # fix the duplicate beaver lake, only keep the one from yoel's samples
  bind_rows(lakes_df[lakes_df$WTRBDKGRPC == "00141HOLB",]) %>%  # add in joe lake
  mutate(tag_prefix = NA,
         watershed = case_when(WTRBDKGRPC == "00141HOLB" ~ "joe",
                               T ~ watershed))
# write.csv(yoel_lakes_df, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/yoel_lakes_df_joining_info.csv", row.names = F)

# warnings from read_excel are fine. Here I drop gill specimens that had an ID issue when the data was being collected. The collector used an excel formula that had an error at some unknown point so all formula IDs are dropped.
tk_parasite_gills_pre <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Parasite dissections/Parasite dissection data.xlsx", sheet = "Gills")
sort(table(tk_parasite_gills_pre$id), decreasing = T)[1:100]
tk_parasite_gills_pre <- tk_parasite_gills_pre[!tk_parasite_gills_pre$`OG ID was formula`,]
sort(table(tk_parasite_gills_pre$id), decreasing = T)[1:10] # duplicate IDs mostly fixed
tk_parasite_id_drop <- names(table(tk_parasite_gills_pre$id)[table(tk_parasite_gills_pre$id) > 1])
tk_parasite_gills <- tk_parasite_gills_pre %>% filter(!id %in% tk_parasite_id_drop) %>% 
  mutate(Unionidae.gill.clam.0.1 = case_when(Unoid > 0 ~ 1,
                               T ~ Unoid),
         Thersitina.gill.copepod.0.1 = case_when(Thersitina > 0 ~ 1,
                                                 T ~ Thersitina))

# warnings from read_excel are fine
tk_parasite_body_cavity <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Parasite dissections/Parasite dissection data.xlsx", sheet = "body.cavity") %>%  # fixed a few errors where glugea was recorded in the anasakis column, I only made changes when the notes said glugea but 0 glugea were recorded. When the note said too many glugea to count I put 10.
  rename(Schistocephalus.cestode.bodycavity.COUNT = schisto,
         Eustrongylides.nematode.bodycavity.COUNT = eustrongylides)
tk_raker_count <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/gill_n_size_correct.csv") %>% 
  mutate(SL.mm = round(SL.mm, 3))
yoel_parasite <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Yoels 2017 project/OG files/02.sexParas.csv") %>% 
  mutate(Unionidae.gill.clam.0.1 = case_when("Unionidae.gill.clam.0.1" == "Few" ~ 1,
                   "Unionidae.gill.clam.0.1" == "Many" ~ 1,
                   "Unionidae.gill.clam.Absent.0.1" == "Absent" ~ 0,
                   #TRUE ~ as.character(Thersitina.gill.copepod.Absent.Few.Many)
                   ),
         Thersitina.gill.copepod.0.1 = case_when("Thersitina.gill.copepod.Absent.0.1" == "Few" ~ 1,
                   "Thersitina.gill.copepod.Absent.0.1" == "Many" ~ 1,
                   "Thersitina.gill.copepod.Absent.0.1" == "Absent" ~ 0,
                   #TRUE ~ as.character(Thersitina.gill.copepod.Absent.Few.Many)
                   )) %>% 
  rename(id = fishID.univ)
merge_yoel_FWA <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/match_yoels_samples_to_FWA.csv")
merge_yoel_tristan <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/link_yoel_and_marius_traits.csv") %>% 
  mutate(Yoels.trait.sc = paste(Yoels.trait,".sc",sep = ""),
         Marius.trait.letter.sc = paste(Marius.trait.letter,".sc",sep = ""))

basin_area <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/basin_area.csv")
lake_and_shoreline_elevations <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/lake_and_shoreline_elevations.csv")
lake_maxtemp <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/lake_maxtemp.csv")
lake_mintemp <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/lake_mintemp.csv")
lake_order <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/lake_order.csv")
lake_precipitation_mean <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/lake_precipitation_mean.csv")
road_length <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/road_length.csv")
shoreline_development <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/shoreline_development.csv")
wtrshd_precipitation_vol <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/wtrshd_precipitation_vol.csv")
icp <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/ICP_wide_format.csv")
ysi_surface <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2021/Summer 2021 notebook data/All projects YSI measurements/analysis ready thesis ysi data.csv") %>% 
  filter(lake.surface.or.bottom == "surface") %>% 
  drop_na(Cond..uS.cm.)
ysi_bottom <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2021/Summer 2021 notebook data/All projects YSI measurements/analysis ready thesis ysi data.csv") %>% 
  filter(lake.surface.or.bottom == "bottom") %>% 
  drop_na(Cond..uS.cm.)

# ===  2. merge enviro predictors into one dataset ==================================================================

# just the GIS variables
gis_predictors <- inner_join(basin_area,lake_and_shoreline_elevations,by = "WTRBDPLD") %>% inner_join(lake_maxtemp,by = "WTRBDPLD") %>% 
  inner_join(lake_mintemp,by = "WTRBDPLD") %>% inner_join(lake_order,by = "WTRBDPLD") %>% inner_join(lake_precipitation_mean,by = "WTRBDPLD") %>%
  inner_join(road_length,by = "WTRBDPLD") %>% inner_join(shoreline_development,by = "WTRBDPLD") %>% inner_join(wtrshd_precipitation_vol,by = "WTRBDPLD") %>% 
  inner_join(lakes_df[,c("WTRBDPLD","AREA_SQM")],by = "WTRBDPLD")

# write.csv(gis_predictors,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/predictors.csv",row.names = FALSE)

ysi_bottom <- ysi_bottom[!ysi_bottom$WTRBDPLD %in% ysi_surface$WTRBDPLD,]
ysi_surface <- bind_rows(ysi_surface,ysi_bottom) # there are 3 lakes where the surface YSI measurement didn't log, so I use the bottom measurement

# checking match between ysi and icp WTRBDPLD
icp$Client.Sample.ID.[!icp$WTRBDPLD %in% ysi_surface$WTRBDPLD]
icp$WTRBDPLD[!icp$WTRBDPLD %in% ysi_surface$WTRBDPLD]
ysi_surface$sit[!ysi_surface$WTRBDPLD %in% icp$WTRBDPLD]
ysi_surface$WTRBDPLD[!ysi_surface$WTRBDPLD %in% icp$WTRBDPLD]


# ===  3. merge stuart 2017 morpho and tristan's collection morpho ==================================================================

# create yoel_data which has morpho and parasites

merge_yoel_tristan <- merge_yoel_tristan %>% filter(In.both.datasets == "y")
yoel_parasite <- yoel_parasite %>% select(fishID.univ, Thersitina.gill.copepod.0.1,
                                          Unionidae.gill.clam.0.1, 
                                          Schistocephalus.cestode.bodycavity.COUNT,
                                          Eustrongylides.nematode.bodycavity.COUNT) 
# yoel_data <- left_join(yoel_morpho, yoel_parasite,by = "fishID.univ") %>% 
#   mutate(id = fishID.univ,
#          dataset = "stuart.nee.2017")


# create tk_data which has morpho and parasites

watershed_WTRBDPLD <- yoel_lakes_df %>% 
  bind_rows(kosciuch_lake_joining_info) %>% 
  select(watershed, WTRBDPLD) 

merged_lengths_sc <- full_join(sc_lateral,sc_ventral, by = c("id","SL.mm","watershed")) %>% 
  mutate(id = str_remove(id,".JPG"),
         SL.mm = round(SL.mm, 3)) %>% 
  left_join(watershed_WTRBDPLD) %>% 
  left_join(gill_raker_n, by = "id")

# write.csv(merged_lengths_sc, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_size_corrected_trait_lengths.csv", row.names = F)
# write.csv(merged_lengths, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_trait_lengths.csv", row.names = F)

merged_lengths <- full_join(lateral,ventral, by = c("id","SL.mm","watershed")) %>% 
  mutate(id = str_remove(id,".JPG"),
         SL.mm = round(SL.mm, 3)) %>% 
  left_join(watershed_WTRBDPLD) %>% 
  left_join(gill_raker_n, by = "id")

mean_merged_lengths_sc <- merged_lengths_sc %>% 
  select(-c(id)) %>% 
  group_by(watershed) %>% 
  summarise_all(mean,na.rm = TRUE) %>% 
  left_join(watershed_WTRBDPLD)

mean_merged_lengths <- merged_lengths %>% 
  select(-c(id)) %>% 
  group_by(watershed) %>% 
  summarise_all(mean,na.rm = TRUE) %>% 
  left_join(watershed_WTRBDPLD)

n_merged_lengths_sc <- merged_lengths_sc %>% 
  select(-c(id)) %>% 
  group_by(watershed) %>% 
  summarise(n = n())

n_merged_lengths <- merged_lengths %>% 
  select(-c(id)) %>% 
  group_by(watershed) %>% 
  summarise(n = n())

mean_merged_lengths_sc <- left_join(mean_merged_lengths_sc, n_merged_lengths_sc)
mean_merged_lengths <- left_join(mean_merged_lengths, n_merged_lengths)

# write.csv(mean_merged_lengths_sc, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_mean_size_corrected_trait_lengths.csv", row.names = F)
# write.csv(mean_merged_lengths, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_mean_trait_lengths.csv", row.names = F)

# t <- left_join(mean_merged_lengths, kosciuch_lake_joining_info, by = "watershed")
# t <- left_join(t, yoel_lakes_df, by = "watershed")
# tk_data <- inner_join(tk_lengths, tk_parasite_body_cavity, by = "id") %>%
#   inner_join(tk_parasite_gills, by = "id") %>%
#   inner_join(tk_raker_count, by = c("id","watershed")) %>%
#   mutate(dataset = "kosciuch.masters")



# 
# # rename columns
# newcols <- colnames(tk_data)
# for(i in 1:length(newcols)){
#   if(newcols[i] %in% merge_yoel_tristan$Marius.trait.letter.sc)
#     newcols[i] <- with(merge_yoel_tristan,Yoels.trait.sc[Marius.trait.letter.sc==newcols[i]])
# } # next change colnames not in Marius.trait.letter.sc
# newcols[newcols =="schisto"] <- "Schistocephalus.cestode.bodycavity.COUNT"
# newcols[newcols =="eustrongylides"] <- "Eustrongylides.nematode.bodycavity.COUNT"
# newcols[newcols =="Unoid"] <- "Unionidae.gill.clam.0.1"
# newcols[newcols =="Thersitina"] <- "Thersitina.gill.copepod.0.1"
# colnames(tk_data) <- newcols

# create a dataset of tk_data with only the same columns as yoel_data
## I split it up because I have two different joining info datasets and there is are different beaver lakes in kosciuch and stuarts datasets
## I need to include parasites here. There is an issue with tk_parasite_gills. There are id duplicates. 
kosciuch_merged_sc <- merged_lengths_sc %>% inner_join(kosciuch_lake_joining_info, by = "watershed") %>% 
  select(!lake_name) %>% 
  left_join(tk_raker_count) %>% 
  left_join(tk_parasite_body_cavity) %>% 
  left_join(tk_parasite_gills, by = "id")
stuart_merged_sc <- merged_lengths_sc %>% inner_join(yoel_lakes_df, by = "watershed") %>% 
  left_join()
  select(colnames(kosciuch_merged_sc))

merged_sc <- bind_rows(kosciuch_merged_sc, stuart_merged_sc)
# %>% 
#   select(c(colnames(yoel_data)[colnames(yoel_data) %in% colnames(tk_data)]), "WTRBDPLD") %>% 
#   bind_rows(yoel_data) %>% 
#   select(-fishID.univ)

# write.csv(merged, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_kosciuch_stuart_morpho.csv", row.names = FALSE)


## lake means, add in parasite data first
lake_means_sc <- merged_sc %>% select(-c(id)) %>% 
  group_by(watershed, WTRBDPLD, tag_prefix, WTRBDKGRPC) %>% 
  summarise_all(mean,na.rm = TRUE)

# write.csv(lake_means, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_kosciuch_stuart_lake_mean_parasitenotsc.csv", row.names=FALSE) # I created merged_kosciuch_stuart_lake_mean_parasitenotsc.csv here when I didn't have the parasites size corrected. After talking with yoel I will leave the parasites not size corrected in my final analysis.

# ===  4. merge notebook and FWA lake data ==================================================================

ormond <- notebook %>% filter(X50k_atlas_waterbody_code == "00316SALM")
not_ormond <- notebook %>% filter(X50k_atlas_waterbody_code != "00316SALM")

## SLP Ormond has WTRBDPLD == 710017713, regular Ormond lake has WTRBDPLD == 710017394
ormond$WTRBDPLD[ormond$tag_prefix == "SLP"] <- 710017713
ormond$WTRBDPLD[ormond$tag_prefix == "OR"] <- 710017394

joined_ormond <- left_join(ormond,lakes,by = "WTRBDPLD")
joined_not_ormond <- left_join(not_ormond,lakes, by = c("X50k_atlas_waterbody_code"="WTRBDKGRPC"), keep = TRUE)

# this notebook_lakes shapefile can be used to link samples to predictors now
notebook_lakes <- bind_rows(joined_not_ormond, joined_ormond)

# drop the geometry column
notebook_lakes <- notebook_lakes %>% select(-geometry)

#write.csv(notebook_lakes,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/notebook_lakes.csv",row.names = FALSE)

## kosciuch notebook per specimen
notebook_WTRBDPLD <- notebook_lakes %>% select(lake_name, tag_prefix, X50k_atlas_waterbody_code, WTRBDPLD) %>% 
  unique() %>% 
  mutate(tag_prefix = toupper(tag_prefix))

# categories
substrate_categories <- unique(unlist(str_split(kosciuch_notebook_specimens$substrate, "")))
substrate_categories <- categories[!categories %in% c(NA, "v")] # not sure what V is (only 1 occurrence), drop NAs
substrate_categories <- paste0(substrate_categories, ".substrate")
bycatch_categories <- c("beetle", "caddisfly", "catfish", "crayfish", "leech", "newt", "odonata", 
                "salamander", "sculpin", "snail", "sunfish", "tadpole", "trout")

# 
notebook_enviro <- inner_join(kosciuch_notebook_specimens_WTRBDPLD, notebook_WTRBDPLD) %>% 
  group_by(WTRBDPLD) %>% 
  summarise(across(all_of(c(substrate_categories)), ~ max(.)), # maybe use max?
            across(all_of(bycatch_categories), ~ max(.)),
            number_of_stickleback = mean(number_of_stickleback, na.rm = T))

yoel_trap_data <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Yoels 2017 project/OG files/150911_EcolData_All.xlsx", sheet = "Trap") %>% 
  mutate(depth_meters_numeric = as.numeric(depth.conv.cm)/100)

 # stuart bycatch categories
table(unlist(str_split(yoel_trap_data$bycatch, "\\.")))
categories <- c("trout", "sculpin", "newt", "char", "barred", "crayfish", "tadpole", "") # I gotta combine some of these
salmonid_synonyms <- c("barred", "char", "salmon", "coho")
categories <- categories[!categories %in% c(NA, "v")] # not sure what V is (only 1 occurrence), drop NAs

one_hot_matrix <- matrix(0, nrow = length(kosciuch_notebook_specimens$bycatch), ncol = length(categories), dimnames = list(NULL, categories))

for (i in 1:length(categories)) {
  one_hot_matrix[, i] <- ifelse(grepl(categories[i], kosciuch_notebook_specimens$substrate), 1, 0)
}

one_hot_substrate <- as.data.frame(one_hot_matrix) %>% 
  rename_with(~ paste0(., ".substrate"), everything())

