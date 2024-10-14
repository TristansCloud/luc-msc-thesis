w322222222222222222222222222222222222ccccccccccccccc                                                                #==== 0.rename.photos.R ======================

# Purpose: remove "-" and " " from photo names.
# Note: cannot have "/" as part of the filepath or filename

#=============================================

library(tidyverse)

#=============================================

# directory containing all the photos to rename. Photos can be in subdirectories
# within the main photo directory. Do not include a trailing "/".
# photo_dir <- "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Landscape_Project_Imaging/Ventral Images"
photo_dir <- "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Landscape_Project_Imaging/Lateral Images"

# list filenames and filepaths
files <- list.files(photo_dir, recursive = TRUE,pattern = ".JPG")
names <- str_split(files,"/")

# directory to record the name changes
name_changes <- data.frame()

for(i in 1:length(files)){
  #remove " " and "-" from filenames
  new_name <- gsub(" ","",tail(names[[i]],1))
  new_name <- gsub("_","",new_name)
  
  old_file <- paste(photo_dir,files[i],sep = "/")
  new_file <- paste(c(photo_dir, paste(names[[i]][1:length(names[[i]])-1]), new_name),collapse = "/")
  
  name_changes <- rbind(name_changes,cbind(old_file,new_file,new_name))
  
  # file.rename(old_file,new_file) # uncomment if you want to run this
}

write.csv(name_changes,paste(photo_dir,"filename_changes.csv",sep = "/"),row.names = FALSE)
