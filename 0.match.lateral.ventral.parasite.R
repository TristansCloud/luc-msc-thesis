# === general ==================================================================

## 0.match.lateral.ventral.parasite.R
## See how many lateral, ventral, and parasite specimens have the same ID.
## Use this to find missing photos and incorrect data. The work done in 0.troubleshooting_photos.R
## is already done, this is the new script to find missing specimens.
#
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.0.2 (2020-06-22)"
#
#
# === table of contents ==================================================================

# 1. load libraries and read in data
#    - lateral lengths
#    - body & gill parasites
#    - ventral photos
#    - gill rakers
#    - armor overlap scoring
#
# 2. check for missing files
#    - In this pass through, I use the lateral images ID as a truth for which samples exist. The lateral images have been vetted.
#
#
#


# ===  1. load libraries and read in data ==================================================================

library(dplyr)
library(readxl)
library(stringr)

# true lake IDs
lake_id <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/notebook_lakes.csv") %>% 
  group_by(lake_name,tag_prefix,X50k_atlas_waterbody_code) %>% 
  summarise(n())

lateral_ids <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/size_corrected_lateral_lengths_mm.csv") %>% 
  select(id); lateral_ids <- str_remove(lateral_ids$id,".JPG")
lateral_ids <- c(lateral_ids, c("CC00144SALM45","CC00144SALM54","MD00254SALM06","XX00103CAMB11")) # these are lateral pictures I took after I found they were missing using this R script 

# warnings are ok
body_parasite <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Parasite dissections/Parasite dissection data.xlsx", sheet = "body.cavity") %>% 
  select(-c(dissector,Notes,date.dissected))
gill_parasite <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Parasite dissections/Parasite dissection data.xlsx", sheet = "Gills") %>% 
  select(-c(dissector, notes)) %>% 
  filter(`OG ID was formula` != TRUE)
gill_parasite_og_ids <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Parasite dissections/OLD parasite data don't use/DONT USE 7 07 2021 Parasite dissection data.xlsx", sheet = "Gills") %>% 
  select(-c(dissector, notes)) %>%
  filter(`OG ID was formula` != FALSE)

gill_parasite <- rbind(gill_parasite,gill_parasite_og_ids) # og_ids are the correct ID for rows that franklin used a formalu to set the ID. Over time extra rows were added which messed the formula ID rows up in Parasite dissection data.xlsx
# the og file reads from DONT USE 7 07 2021 Parasite dissection data.xlsx which has the correct formula IDs.

# parse through the photo files because no trustworthy file of IDs exist
ventral_ids <- vector()
path <- "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work"
for(i in 1:10){
  folder <- paste("ventral images set",i)
  ventral_ids <- c(ventral_ids,list.files(file.path(path,folder),pattern = ".JPG"))
}

# there is a folder called Missing ventral images and I want to see if these are being landmarked
ventral_potential_not_included <- list.files("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Missing ventral images",pattern = ".JPG")
# confirm all IDs are unique
length(ventral_potential_not_included) == length(unique(ventral_potential_not_included))

# gill rakers and plate overlap, which is on the same page
gill_rakers <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Gill Rakers/Gill Raker Counts + Measures & Gravid Scoring.xlsx", sheet = "gill.raker.counts") %>% 
  select(id,right.side.gillraker.number)

bp_overlap <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Gill Rakers/Gill Raker Counts + Measures & Gravid Scoring.xlsx", sheet = "bp.lp.overlap") %>% 
  select(id,ds1.bp.overlap.score,ds2.bp.overlap.score)

lp_overlap <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Gill Rakers/Gill Raker Counts + Measures & Gravid Scoring.xlsx", sheet = "lp.lp.overlap") %>% 
  select(id,overlap.scoring )

# any IDs in this list do not have to be checked because I am dropping them from 
# an analysis where fish need to be matched. There were too many ambiguous decisions about which sample goes with what tag for
# these fish. Its also possible these samples don't exist.
problem_id <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/lateral_specimens_that_shouldnt_exist.csv")

# ===  2. check for missing files ==================================================================
# the results I indicate in comments may change as I fix errors
# each section (parasite, gill, photo) goes through 4 main steps:
#    1. duplicated IDs
#    2. NA values
#    3. IDs not in the lateral photos (spelling error samples or missing lateral photos)
#    4. lateral photo IDs not in the dataset (missing samples)


##### Parasite data

# check body cavity (maybe gill raker too) dissections done by franklin on 2021-06-04: I checked some of 2021-06-04 and seems ok
# check body cavity FW00241SALM16

#    1. duplicated IDs
body_parasite_id_to_check <- sort(names(table(body_parasite$id)[table(body_parasite$id) > 1]))
gill_parasite_id_to_check <- sort(names(table(gill_parasite$id)[table(gill_parasite$id) > 1]))

#    2. NA values
body_parasite_id_to_check <- c(body_parasite_id_to_check, body_parasite[!complete.cases(body_parasite),]$id)
gill_parasite_id_to_check <- c(gill_parasite_id_to_check, gill_parasite[!complete.cases(gill_parasite),]$id)

#    3. IDs not in the lateral photos (spelling error samples or missing lateral photos)
potential_lateral_missing <- body_parasite$id[!body_parasite$id %in% lateral_ids]
potential_lateral_missing <- c(potential_lateral_missing, gill_parasite$id[!gill_parasite$id %in% lateral_ids])
gill_parasite_og_ids$id[!gill_parasite_og_ids$id %in% lateral_ids] # these will be contained in gill_parasite

View(lake_id) # use this to fix spelling errors

#    4. lateral photo IDs not in the dataset (missing samples)
body_parasite_id_to_check <- c(body_parasite_id_to_check,lateral_ids[!lateral_ids %in% body_parasite$id])
gill_parasite_id_to_check <- c(gill_parasite_id_to_check,lateral_ids[!lateral_ids %in% gill_parasite$id])

# print error samples to check
body_parasite_id_to_check <- sort(unique(body_parasite_id_to_check))
gill_parasite_id_to_check <- sort(unique(gill_parasite_id_to_check))

# 
#write.csv(cbind(body_parasite_id_to_check,gill_parasite_id_to_check),"C:/Users/trist/Documents/gill_body_parasite_fixing.csv")


##### Ventral photos

# check if that missing ventral folder is being landmarked
table(ventral_potential_not_included %in% ventral_ids)
# 40 FALSE it is not being landmarked!! I will have someone landmark this folder

all_ventral <- c(ventral_ids,ventral_potential_not_included); all_ventral <- str_remove(all_ventral,".JPG")

#    1. duplicated IDs
table(all_ventral)[table(all_ventral) > 1]

#    2. NA values # not applicable as I read in filenames that included .JPG

#    3. IDs not in the lateral photos (spelling error samples or missing lateral photos)
potential_lateral_missing <- c(potential_lateral_missing, all_ventral[!all_ventral %in% lateral_ids])

#    4. lateral photo IDs not in the dataset (spelling error samples or missing samples). These photos need to be taken
lateral_ids[!lateral_ids %in% all_ventral]


##### Lateral and basal plate overlap

#    1. duplicated IDs
lp_overlap_to_check <- table(lp_overlap$id)[table(lp_overlap$id) > 1]
# EC00126CAMB25 EC00126CAMB33 XX00400SALM01 XX00400SALM02 all have two occurrences
bp_overlap_to_check <- table(bp_overlap$id)[table(bp_overlap$id) > 1]
# BB00096SALM37 BK00304VICT08 CC00144SALM48 all have two occurrences

#    2. NA values
lp_overlap_to_check <- c(lp_overlap_to_check, lp_overlap[!complete.cases(lp_overlap),]$id)
bp_overlap_to_check <- c(bp_overlap_to_check, bp_overlap[!complete.cases(bp_overlap),]$id)

#    3. IDs not in the lateral photos (spelling error samples or missing lateral photos)
# there are no missing lateral images, but there are some spelling errors to fix
lp_overlap$id[!lp_overlap$id %in% lateral_ids]
bp_overlap$id[!bp_overlap$id %in% lateral_ids]


#    4. lateral photo IDs not in the dataset (missing samples)
lp_overlap_to_check <- c(lp_overlap_to_check,lateral_ids[!lateral_ids %in% lp_overlap$id])
bp_overlap_to_check <- c(bp_overlap_to_check,lateral_ids[!lateral_ids %in% bp_overlap$id])


sort(unique(lp_overlap_to_check))
sort(unique(bp_overlap_to_check))


##### gill raker counts
# Done with gill raker count corrections! Any IDs here can't be adressed unfortunately


# #    1. duplicated IDs
# gill_rakers_to_check <- names(table(gill_rakers$id)[table(gill_rakers$id) > 1])
# 
# #    2. NA values
# gill_rakers_to_check <- c(gill_rakers_to_check,gill_rakers[!complete.cases(gill_rakers),]$id)
# 
# #    3. IDs not in the lateral photos (spelling error samples or missing lateral photos)
# potential_lateral_missing <- c(potential_lateral_missing, gill_rakers$id[!gill_rakers$id %in% lateral_ids])
# 
# #    4. lateral photo IDs not in the dataset (spelling error samples or missing gill raker samples)
# gill_rakers_to_check <- c(gill_rakers_to_check,lateral_ids[!lateral_ids %in% gill_rakers$id])
# 
# sort(unique(gill_rakers_to_check))

#### potentially missing lateral photos
# After these may also be missing from other datasets
sort(unique(potential_lateral_missing))
# Fish from potential_lateral_missing I couldn't find when looking at our samples. These samples are likely lost or never existed:
# BB00096SALM01 -- damaged fish
# BB00096SALM05 -- damaged fish
# BB00096SALM06 -- damaged fish
# BB00096SALM50 -- damaged fish
# BK00304VICT43
# BW00292SALM25 -- damaged fish
# FR00078CAMB19 -- damaged fish
# LS00090SALM37
# MD00254SALM66
# MH00078SALM17 -- damaged fish
# MN00158SALM49 - Whatever has this probably should be MN00158SALM44
# XX00403SALM35 -- damaged fish







