# === general ==================================================================

## This script is done, new work finding missing specimens is in 0.match.lateral.ventral.parasite.R

## 0.troubleshooting_photos.r
## look for missing morphometric photos
#
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.0.2 (2020-06-22)"
#
#
# === table of contents ==================================================================

# 1. lateral images - already ran
#
# 2. both lateral and ventral
#
# 3. ventral images

# libraries
library(tidyverse)

# ===  1. lateral images ==================================================================

#### Check photo names for errors ####
lat_dir <- "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Lateral images"
lat_photos <- list.files(lat_dir, recursive = FALSE,pattern = ".JPG")

# remove ".JPG" from image names
fullname <- substr(lat_photos, 1, nchar(lat_photos)-4)
number <- substr(lat_photos, nchar(fullname)-1, nchar(fullname))
unique(number)
# fullname[number == "EL"]# removed this photo, it was called "XX00103CAMB UNCLEAR LABEL"
# fullname[number == "PG"]# renamed this photo from "BW00292SALM12JPG" to BW00292SALM12

# read in notebook data to compare site names
notebook <- read.csv("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Summer 2020/Summer 2020 datasheets/2020 tkosciuch fieldwork notebook data/2020 tkosciuch fieldwork notebook.csv")

true_lakes <- paste(notebook$tag_prefix, notebook$X50k_atlas_waterbody_code, sep = "")
photo_lakes <- substr(lat_photos, 1, nchar(lat_photos)-6)

# find out which names are errors
photo_error <- toupper(sort(unique(photo_lakes))) %in%
toupper(sort(unique(true_lakes))) 

error_photo_lakes <- toupper(sort(unique(photo_lakes)))[!photo_error]

# rename the photos
lat_photos[photo_lakes %in% error_photo_lakes] # "CC01044SALM33.JPG" "CG00144SALM14.JPG" "LS0090SALM23.JPG"  "XX0089SALM12.JPG" 
# file.rename(paste(lat_dir,"CC01044SALM33.JPG", sep = "/"),paste(lat_dir,"CC00144SALM33.JPG", sep = "/"))
# file.rename(paste(lat_dir,"CG00144SALM14.JPG", sep = "/"),paste(lat_dir,"CC00144SALM14.JPG", sep = "/"))
# file.rename(paste(lat_dir,"LS0090SALM23.JPG", sep = "/"),paste(lat_dir,"LS00090SALM23.JPG", sep = "/"))
# file.rename(paste(lat_dir,"XX0089SALM12.JPG", sep = "/"),paste(lat_dir,"XX00089SALM12.JPG", sep = "/"))

# for franklin to look for missing photos
lat_photos <- list.files(lat_dir, recursive = FALSE,pattern = ".JPG")
# write.csv(sort(lat_photos),"C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/lat_photo_names.csv")

# automating finding missing photos
fullname <- substr(lat_photos, 1, nchar(lat_photos)-4)
photo_lakes <- substr(lat_photos, 1, nchar(lat_photos)-6)
number <- as.numeric(substr(lat_photos, nchar(fullname)-1, nchar(fullname)))
fmphoto <- data.frame(photo_lakes,number) %>% 
  pivot_wider(names_from = number, values_from = number)

fmphoto2 <- fmphoto[ , order(as.numeric(names(fmphoto)))]
# write.csv(fmphoto2, "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/lat_photo_pivot_wide.csv")


# ===  2. both lateral and ventral  ==================================================================

## Check which ventral image names are not in the the lateral images, and which 
## lateral images are not in the ventral images
library(tidyverse)

lateral <- list.files("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Lateral images", recursive = FALSE,pattern = ".JPG")
ventral <- list.files("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Ventral images", recursive = FALSE,pattern = ".JPG")

missing_lat <- setdiff(ventral,lateral)
missing_vent <- setdiff(lateral,ventral)

print(missing_lat)
print(missing_vent)

# write.csv(missing_lat,"C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/more_lateral_missing_pics.csv")
# write.csv(missing_vent,"C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/more_ventral_missing_pics.csv")

## compare the ventral and lateral pivot csv (ventral csv created in section 3, run that first)
## to the notebook data. I created two columns in the notebook data manually based on the sample number notes
notebook_filtered <- read.csv("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Summer 2020/Summer 2020 datasheets/2020 tkosciuch fieldwork notebook data/2020 tkosciuch fieldwork notebook.csv") %>% 
  filter(!is.na(sample_low_number) & !is.na(sample_high_number) & !is.na(sample_numbers_also_included) & !is.na(sample_numbers_not_included) & !is.na(sample_notes)) # drop rows that did not have any samples in them
# creating a new data frame to mimic the pivot csvs
#up to 70 fish per lake
lake_names <- unique(toupper(paste(notebook_filtered$tag_prefix,notebook_filtered$X50k_atlas_waterbody_code, sep = "")))
all_fish <- cbind(data.frame(names = lake_names),matrix(nrow = length(lake_names), ncol = 70))
# iterate each row
for(i in 1:dim(notebook_filtered)[1]){
  temp <- notebook_filtered[i,]
  name <- toupper(paste(temp$tag_prefix, temp$X50k_atlas_waterbody_code, sep = ""))
  
  if(!is.na(temp$sample_high_number)){
    samples <- temp$sample_low_number:temp$sample_high_number
  } else {
    samples <- temp$sample_low_number
  }
  # add additional samples if needed, str_split works with numeric and integers as well as characters, no need to convert numeric to character for it to work
  if(temp$sample_numbers_also_included != ""){ # originally used !is.na() but blank cells are "", which doesn't register as NA
    add_samples <- str_split(temp$sample_numbers_also_included, pattern = "-")[[1]]
  if(length(add_samples == 1)){ # needed in case add samples is only one number
    add_samples <- as.numeric(add_samples[[1]][1])
  } else {
    add_samples <- as.numeric(add_samples[[1]][1]):as.numeric(add_samples[[1]][-1])
  }
  samples <- c(samples, add_samples)
  }
  # remove samples if needed
  if(!is.na(temp$sample_numbers_not_included)){
    drop_samples <- str_split(temp$sample_numbers_not_included, pattern = "-")[[1]]
    if(length(drop_samples == 1)){ # needed in case add samples is only one number
      drop_samples <- as.numeric(drop_samples[[1]][1])
    } else {
      drop_samples <- as.numeric(drop_samples[[1]][1]):as.numeric(drop_samples[[1]][-1])
    }
    samples <- samples[!samples %in% drop_samples]
  }
  
  all_fish[all_fish$names == name,samples+1] <- samples # +1 within the brackets to skip the names column
}

# write.csv(all_fish, "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/all_fish.csv")

# now compare the all_fish dataframe to the lat and ventral pivot_wider dataframes
# it's necessary to read in from CSV so the sample number column names have an "X"
# appended to make them readable as characters 
fmphoto2 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/lat_photo_pivot_wide.csv", row.names = 1)
vent_fmphoto2 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/vent_photo_pivot_wide.csv")
all_fish <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/all_fish.csv", row.names = 1)


# compare fmphoto2 to all_fish
check_missing <- function(x,y,id_col,sample_cols){ # y must be the all_fish dataframe. id_col is the integer that references the lake name column in the fmphoto2 dataframe. sample_cols is a vector giving the range of columns that contain the sample information, with the lowest number being sample 1.
  # this will hold the missing fish
  missing_fish <- c() 
  # make sure the first col of x is the ID col
  temp_data <- cbind(x[,id_col],x[,sample_cols])
  # find which lake names are not in the actual samples (spelling errors most likely)
  incorrect_lakes <- temp_data[,1][!temp_data[,1] %in% y$names]
  # find which lake names are missing in the temp_data
  missing_lakes <- y$names[!y$names %in% temp_data[,1]]
  # drop the missing lakes from all_fish to prevent errors with calling a missing column
  y <- y[!y$names %in% missing_lakes,]
  # iterate through each row in all_fish, which is each lake
  for(i in 1:dim(y)[1]){ 
    # drop any NA columns 
    temp_row <- y[i,colSums(is.na(y[i,]))<1] 
    
    ## record which samples are in temp_data but not all_fish to check for
    ## missing samples in all_fish, which is based on my 2020 notebook
    # keep only the same lake as the all_fish row
    temp_photo_row <- temp_data[temp_data[,1] %in% temp_row$names,] 
    # drop NA columns
    temp_photo_row <- temp_photo_row[,colSums(is.na(temp_photo_row[1,]))<1] 
    # store the photo sample numbers
    temp_photo_samples <- colnames(temp_photo_row)[2:length(temp_photo_row[1,])]
    # store the notebook sample numbers for the same lake
    nb_sample_num <- colnames(temp_row)[2:length(temp_row[1,])]
    # find samples in the photos but not in all_fish. The goal is to update my
    # notebook if some sample notes were missing but I imagine this will be rare.
    photo_extra_samples <- temp_photo_samples[!temp_photo_samples %in% nb_sample_num]
    # if i = 1, create a list to store sample numbers from the photos not present in my notebook (nb means notebook)
    if(i == 1){nb_missing_samples <- list(c(temp_row[1,1],photo_extra_samples))}
    else{nb_missing_samples[[i]] <- c(temp_row[1,1],photo_extra_samples)}
    # find samples in my notebook but not the photos
    nb_extra_samples <- nb_sample_num[!nb_sample_num %in% temp_photo_samples]
    # if i = 1, create a list to store sample numbers from the notebook which I do not have a photo for
    if(i == 1){photo_missing_samples <- list(c(temp_row[1,1],nb_extra_samples))}
    else{photo_missing_samples[[i]] <- c(temp_row[1,1],nb_extra_samples)}
    
  }

  return(list(incorrect_lakes, missing_lakes, nb_missing_samples, photo_missing_samples))
  }

check_missing(fmphoto2, all_fish, 65, c(1:64))
check_missing(vent_fmphoto2, all_fish, 1, c(2:65))

## After running this, I went and fixed the errors presented:
# 

# ===  3. ventral images ==================================================================

#### List photos in folder ####
vent_photo_dir <- "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Ventral images"

vent_photos <- list.files(vent_photo_dir, recursive = TRUE,pattern = ".JPG")

# remove ".JPG" from image names
vent_fullname <- substr(vent_photos, 1, nchar(vent_photos)-4)
vent_number <- substr(vent_photos, nchar(vent_fullname)-1, nchar(vent_fullname))
unique(vent_number)

# read in notebook data to compare site names
notebook <- read.csv("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Summer 2020/Summer 2020 datasheets/2020 tkosciuch fieldwork notebook data/2020 tkosciuch fieldwork notebook.csv")

true_lakes <- paste(notebook$tag_prefix, notebook$X50k_atlas_waterbody_code, sep = "")
vent_photo_lakes <- substr(vent_photos, 1, nchar(vent_photos)-6)

# find out which names are errors
vent_photo_error <- toupper(sort(unique(vent_photo_lakes))) %in%
  toupper(sort(unique(true_lakes))) 

vent_error_photo_lakes <- toupper(sort(unique(vent_photo_lakes)))[!vent_photo_error] # only two name errors and both are not spelling mistakes, the tags were unclear or they were potential duplicates and this info was included in the filename leading to a mismatch. The spelling errors were "POTENTIAL DUPLICATES/LM00232SALM UNKNO" and "XX00103CAMB(LABEL UNCLE" 

# automating finding missing photos
vent_fullname <- substr(vent_photos, 1, nchar(vent_photos)-4)
vent_photo_lakes <- substr(vent_photos, 1, nchar(vent_photos)-6)
vent_number <- as.numeric(substr(vent_photos, nchar(vent_fullname)-1, nchar(vent_fullname)))
vent_fmphoto <- data.frame(vent_photo_lakes,vent_number) %>% 
  pivot_wider(names_from = vent_number, values_from = vent_number)

vent_fmphoto2 <- vent_fmphoto[ , order(as.numeric(names(vent_fmphoto)))]
# I use this csv to compare with my notebook and see if any fish are missing
write.csv(vent_fmphoto2, "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/vent_photo_pivot_wide.csv")

