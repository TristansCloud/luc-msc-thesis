# === general ==================================================================

# 0.split.ventral.pics.R
# separate ventral pics for undergrads to landmark

# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.0.2 (2020-06-22)"
#
#

# ==============================================================================

folder_directory <- "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/ventral images set"
# directory containing all the photos to split up. 
photo_dir <- "C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Ventral images"

# list filenames and filepaths
files <- list.files(photo_dir,pattern = ".JPG")

split_vect <- rep(c(1:10),length.out = length(files))

for(i in 1:length(files)){
  file.copy(from = paste(photo_dir, files[i], sep = "/"),
            to = paste(folder_directory,split_vect[i]))
}
