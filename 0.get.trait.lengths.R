# === general ==================================================================

## 0.get.trait.lengths.R
## take landmark coordinates and output feature lengths. Write the feature
## lengths to a CSV file
#
#
# Note for author use: this script used to be called 0.format.ObjectJ.results.R
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
# 2. Check 15mm scale
#
# 3. Calculate feature lengths
#
# 4. Append notes on missing features and rigor mortis


# ===  1. load libraries and import data ==================================================================

library(tidyverse)

datapath <- "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/OJJ file backups"

## Lateral image data import. lateral_missing are samples that were found to be missing from the lateral dataset so need to be added in.
## there is more info in the README in the datapath.
lateral <- read.csv(file.path(datapath,"231029lateral_without_missing_fish.csv")) # , fileEncoding = 'UTF-8-BOM'
lateral_missing <- read.csv(file.path(datapath,"Dec10lateralmissingfish.csv")) # , fileEncoding = 'UTF-8-BOM'

# lateral landmarks used to test ML-morph's performance
mlmorph <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/ML-morph chapter/5foldCV/mlmorph.csv") %>% 
  mutate(id = stringi::stri_replace_all_regex(id, pattern = "./test/", replacement = "")) %>% 
  mutate(id = stringi::stri_replace_all_regex(id, pattern = "jpg", replacement = "JPG"))

# notes on rigor mortis and  missing features
notes_lateral <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Notes on lateral photos.csv")

# ventral landmarks
v1 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/Set 1 Results 1.19.22 - Sidney Ryan.csv")
v2 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/Set2Results.csv")
v3 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/Set 3 Results Sidney Ryan.csv")
v4 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/Set4Results.csv")
v5 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/Set5Results.csv")
v6 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/Set6Results.csv")
v7 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/Set7Results.csv")
v8 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/Set8Results.csv")
v9 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/Set9Results.csv")
v10 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/Tristan ventral set 10 Dec 17.csv")
v11 <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Landmarking my fish/images to landmark/Undergrads ventral work/Landmark results/missingimages_results.csv")
ventral <- bind_rows(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11) ; rm(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11)

# ===  2. Check 15mm scale ==================================================================

# as a correctly set scalebar is critical to correctly size correcting, I check
# the largest and smallest scale lengths (in pixels) to ensure they are properly
# placed.

lateral.scale.check <- lateral[,c("file","scale.15.mm")]
ventral.scale.check <- ventral[,c("file","scale.15.mm")]

# use the sort function in View and check the 5 largest and smallest scales
View(ventral.scale.check)
# smallest and largest scales by pixels
# LM00232SALM02.JPG 757.3520 # checked, ok
# MT00294VICT33.JPG 1417.401 # checked, ok

View(lateral.scale.check)

# I found some errors! Only in the largest scales, the smallest were correct. I checked each sample from largest 
# to smallest scale in pixels until I found no errors for at least 5 photos (I actually checked ~10 more pics)

# I fixed the scale errors I found, and lateral_except_missing_fish.csv and Dec11lateral_without_missing_fish.csv
# have the corrected scale bars. Mar1lateral_without_missing_fish also has corrected scalebars. I believe lateral_except_missing_fish.csv is the same as Dec11lateral_without_missing_fish.csv

# I checked all of the lateral_missing photos, no scale errors there

# ===  3. Calculate feature lengths ==================================================================

# Note for the author: adapted from 2.mlmorph.error.R section 6

## Lateral lengths

# Original order of lateral landmarks
#[1] "1.x"           "1.y"           "2.x"           "2.y"           "3.x"           "3.y"           "4.x"           "4.y"          
#[9] "5.x"           "5.y"           "6.x"           "6.y"           "7.x"           "7.y"           "8.x"           "8.y"          
#[17] "9.x"           "9.y"           "10.x"          "10.y"          "11.x"          "11.y"          "12.x"          "12.y"         
#[25] "13.x"          "13.y"          "14.x"          "14.y"          "15.x"          "15.y"          "16.x"          "16.y"         
#[33] "17.x"          "17.y"          "18.x"          "18.y"          "19.x"          "19.y"          "20.x"          "20.y"         
#[41] "21.x"          "21.y"          "22.x"          "22.y"          "23.x"          "23.y"          "24.x"          "24.y"         
#[49] "25.x"          "25.y"          "26.x"          "26.y"          "27.x"          "27.y"          "28.x"          "28.y"         
#[57] "29.x"          "29.y"          "30.x"          "30.y"          "31.x"          "31.y"          "d.posterior.x" "d.posterior.y"
#[65] "h.ventral.x"   "h.ventral.y"   "k.tip.x"       "k.tip.y"       "l.tip.x"       "l.tip.y"       "n.dorsal.x"    "n.dorsal.y"   
#[73] "n.ventral.x"   "n.ventral.y"   "r.ventral.x"   "r.ventral.y"   "s.anterior.x"  "s.anterior.y"  "s.posterior.x" "s.posterior.y"

# calculate the euclidean distance between two points
euclid_distance <- function(x1, y1, x2, y2){
  return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
}


# get.lengths() takes my manual or MLmorph landmarks and finds the trait lengths
# described in Marius' landmarking protocol.
# Note for functionality: lm.data must have a column called "id" containing the sample ID
# These are the subsections within get.lengths():
#    1. Rename OJJ results
#       - The colnames of the OJJ files are labelled differently than what MLMorph
#         outputs. This will match the MLmorph formatting. Note that MLMorph will
#         have 0 indexed landmarks but section 1 assumes 1 indexing, not 0.
#         The code under #ml morph's landmarks are 0 indexed, so subtract 1
#         accounts for this 0 index.
#    2. Iterate through each sample and trait length
#
#    3. Find the trait lengths  

## Note for future: I can remove n.lengths argument if I want to and calculate it in the function.
##    This function could be sped up considerably by following how I used euclid_distance() in 2.mlmorph.error.R
get.lengths <- function(n.lengths,lm.data,lm.view = "lateral",mlmorph = FALSE){
  
  colnames(lm.data)[1] <- "id"
  
  if(lm.view == "ventral"){
    
    # DF to store trait lengths, add two columns for photo id and scalebar
    trait_lengths <- data.frame(matrix(ncol = n.lengths +2, nrow = dim(lm.data)[1])) 
    
    trait_lengths[,1:2] <- c(lm.data$id,lm.data$scale.15.mm)
    
    colnames(trait_lengths[,1:2]) <- c("id","scale.15.mm")
    
    for(i.letter in letters[1:n.lengths]){
      print(i.letter)
      
      if(i.letter == "a"){
        a <- with(lm.data, euclid_distance(a.anterior.x, a.anterior.y, a.posterior.x, a.posterior.y))
      }
      if(i.letter == "b"){
        b <- with(lm.data, euclid_distance(b.anterior.x, b.anterior.y, b.posterior.x, b.posterior.y))
      }
      if(i.letter == "c"){
        c <- with(lm.data, euclid_distance(c.fishs.left.x,c.fishs.left.y,c.fishs.right.x,c.fishs.right.y))
      }
      if(i.letter == "d"){
        d <- with(lm.data, euclid_distance(d.fishs.left.x,d.fishs.left.y,d.fishs.right.x,d.fishs.right.y))
      }
      if(i.letter == "e"){
        e <- with(lm.data, euclid_distance(e.fishs.left.x,e.fishs.left.y,e.fishs.right.x,e.fishs.right.y))
      }
      if(i.letter == "f"){
        f <- with(lm.data, euclid_distance(f.fishs.left.x,f.fishs.left.y,f.fishs.right.x,f.fishs.right.y))
      }
      if(i.letter == "g"){
        g <- with(lm.data, euclid_distance(g.fishs.left.x,g.fishs.left.y,g.fishs.right.x,g.fishs.right.y))
      }
      if(i.letter == "h"){
        h <- with(lm.data, euclid_distance(h.anterior.x, h.anterior.y, h.posterior.x, h.posterior.y))
      }
      if(i.letter == "i"){
        i <- with(lm.data, euclid_distance(i.anterior.x, i.anterior.y, h.posterior.x, h.posterior.y)) # uses the posterior end of h
      }
      if(i.letter == "j"){
        j <- with(lm.data, euclid_distance(j.fishs.left.x,j.fishs.left.y,j.fishs.right.x,j.fishs.right.y))
      }
      if(i.letter == "k"){
        k <- with(lm.data, euclid_distance(k.fishs.left.x,k.fishs.left.y,k.fishs.right.x,k.fishs.right.y))
      }
      if(i.letter == "l"){
        l <- with(lm.data, euclid_distance(l.anterior.x, l.anterior.y, l.posterior.x, l.posterior.y))
      }
      if(i.letter == "m"){
        m <- with(lm.data, euclid_distance(m.anterior.x, m.anterior.y, m.posterior.x, m.posterior.y))
      }
      
    }
    
    trait_lengths[,-c(1:2)] <- cbind(a,b,c,d,e,f,g,h,i,j,k,l,m)
    colnames(trait_lengths) <- c("id","scale.15.mm",letters[1:n.lengths])
  }
  
  
  
  # I wrote the lateral trait length calculations to work for both ML-morph
  # and ObjectJ and the script is a bit convoluted as a result. The results are
  # correct, I've checked! All of the `lm.view == lateral` if statements
  # are unnecessary, earlier I intended the script below to work for both lateral
  # and ventral, but realized I could write it more efficiently (see ventral) and
  # didn't want to remove the extra if statements.
  if(lm.view == "lateral"){
    
    if(mlmorph == FALSE){
      
      if(lm.view == "lateral"){
        # vector of armor plate count and scale bar length in pixels, will be appended to the results
        armor.plate.count <- lm.data$t.armor.plate.count
        scale.15.mm <- lm.data$scale.15.mm
        # drop the armor and scale column
        lm.data <- lm.data[,!colnames(lm.data) %in% c("t.armor.plate.count","scale.15.mm")]
        n_landmarks <- dim(lm.data)[2] - 1 # -1 because the id column. This assumes every other column is a landmark!
      }
    }
    
    if(mlmorph == TRUE){
      # vectors to be appended back later
      box.width <- lm.data$box_width
      box.height <- lm.data$box_height
      # drop the non landmark or id columns
      lm.data <- lm.data[,!colnames(lm.data) %in% c("box_id","box_top","box_left","box_width","box_height")]
      
      n_landmarks <- dim(lm.data)[2] - 1 # -1 because the id column. This assumes every other column is a landmark!
    }
    
    # divide by two because each landmark has an X and Y
    n_landmarks <- n_landmarks/2
    
    XYs <- rep(c("X","Y"),n_landmarks)
    
    # iterate each landmark
    for(i in 1:(n_landmarks)){
      
      iY <- i * 2
      iX <- iY -1
      
      XYs[iY] <- paste(XYs[iY],i,sep="")
      XYs[iX] <- paste(XYs[iX],i,sep="")
      
    }
    # rename every column except the first. The first column is the id
    colnames(lm.data)[-1] <- XYs
    
    
    ### 2. Iterate through each sample and trait length
    
    # DF to store trait lengths, add one column for photo id
    trait_lengths <- data.frame(matrix(ncol = n.lengths +1, nrow = 0)) 
    colnames(trait_lengths) <- c("id",LETTERS[1:n.lengths])
    
    # replace the last column name with SL if the photos are lateral
    if(lm.view == "lateral"){
      colnames(trait_lengths)[n.lengths+1] <- "SL" # n.lengths + 1 because the id column is still in trait_lengths
    }
    # output to judge runtime
    print(c("samples to process:",length(lm.data$id))) ; n.sample <- 0
    # for each sample. # could probably rewrite this to vectorize and run faster like was done for ventral lengths
    for(j in lm.data$id){
      n.sample <- n.sample + 1
      cat(n.sample, " ")
      # single sample
      temp_lm.data <- lm.data[lm.data$id==j,]
      # empty vector to store trait lengths
      temp_lengths <- rep(0,n.lengths)
      
      z <- 0 
      # Use LETTERS as the trait lengths are identified by capital letters
      for(i in LETTERS[1:n.lengths]){ # 
        z <- z+1
        # these if statements set the landmarks used to get trait lengths
        if(lm.view == "lateral"){
          if(i == "A"){
            lm.1 <- 1
            lm.2 <- 22
          } else if(i == "B"){
            lm.1 <- 3
            lm.2 <- 23
          } else if(i == "C"){
            lm.1 <- 1
            lm.2 <- 6
          } else if(i == "D"){
            lm.1 <- 1
            lm.2 <- 32
          } else if(i == "E"){
            lm.1 <- 1
            lm.2 <- 8
          } else if(i == "F"){
            lm.1 <- 6
            lm.2 <- 7
          } else if(i == "G"){
            lm.1 <- 22
            lm.2 <- 21
          } else if(i == "H"){
            lm.1 <- 5
            lm.2 <- 33
          } else if(i == "I"){
            lm.1 <- 8
            lm.2 <- 21
          } else if(i == "J"){
            lm.1 <- 9
            lm.2 <- 20
          } else if(i == "K"){
            lm.1 <- 9
            lm.2 <- 34
          } else if(i == "L"){
            lm.1 <- 10
            lm.2 <- 35
          } else if(i == "M"){
            lm.1 <- 11
            lm.2 <- 18
          } else if(i == "N"){
            lm.1 <- 36
            lm.2 <- 37
          } else if(i == "O"){
            lm.1 <- 11
            lm.2 <- 12
          } else if(i == "P"){
            lm.1 <- 17
            lm.2 <- 16
          } else if(i == "Q"){
            lm.1 <- 20
            lm.2 <- 18
          } else if(i == "R"){
            lm.1 <- 28
            lm.2 <- 38
          } else if(i == "S"){
            lm.1 <- 39
            lm.2 <- 40
          } else if(i == "T"){ # this is the standard length
            lm.1 <- 1
            lm.2 <- 14
          }
        }
        
        
        ### 3. Find the trait lengths 
        
        # create landmark character to reference the lm.data's colnames
        lm.1.x <- paste("X", lm.1, sep = "")
        lm.1.y <- paste("Y", lm.1, sep = "")
        lm.2.x <- paste("X", lm.2, sep = "")
        lm.2.y <- paste("Y", lm.2, sep = "")
        
        
        # pull out x and y coords of each landmark needed
        lm.1.x <- unname(unlist(temp_lm.data %>% select(all_of(lm.1.x))))
        lm.1.y <- unname(unlist(temp_lm.data %>% select(all_of(lm.1.y))))
        lm.2.x <- unname(unlist(temp_lm.data %>% select(all_of(lm.2.x))))
        lm.2.y <- unname(unlist(temp_lm.data %>% select(all_of(lm.2.y))))
        
        # correct for the adjustments to the Y coordinate made by Simple-ML-morph
        # which sets the bounding box two pixels in from the edge of the photo
        # This might not be needed now that I removed the image height (in pixels) - the original Y coordinate
        # transformation in line 104 of utils.py in ML-morph. For more info refer to
        # `Log of Tristan's work on ML-morph.Rmd.` As both coordinates are equally transformed
        # I'm not sure it makes a difference if I do ml_bb_height - lm when calculating 
        # ML-morph trait lengths.
        # if(mlmorph == TRUE){
        #   ml_bb_height <- temp_lm.data$box_height 
        #   lm.1.y <- ml_bb_height - lm.1.y + 2
        #   lm.2.y <- ml_bb_height - lm.2.y + 2
        # }
        
        # calculate the trait length
        temp_lengths[z] <- euclid_distance(lm.1.x,lm.1.y,lm.2.x,lm.2.y)
        
      }
      
      # bind the fishes trait length to the trait length dataframe. J is an iterator of the sample's id
      trait_lengths[nrow(trait_lengths)+1,] <-  c(j, temp_lengths)
      
    }
    
    # ensure all trait length columns are numeric
    trait_lengths[,-1] <- trait_lengths[,-1] %>% mutate_if(is.character,as.numeric)
    
    # add the scale and armor columns that were dropped earlier
    if(mlmorph == FALSE){
      if(lm.view == "lateral"){
        trait_lengths <- cbind(trait_lengths,armor.plate.count,scale.15.mm)
      } else {
        trait_lengths <- cbind(trait_lengths,scale.15.mm)
      }}
    if(mlmorph == TRUE){
      trait_lengths <- cbind(trait_lengths,box.width,box.height)
    }
    
  }
  return(trait_lengths)
}

# ===  3. Calculate feature lengths ==================================================================

outputpath <- "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data"


#### Ventral

ventral_lengths <- get.lengths(n.lengths = 13, lm.data = ventral, lm.view = "ventral",mlmorph = FALSE)

# write.csv(ventral_lengths,file.path(outputpath,"ventral_lengths.csv"), row.names = FALSE)

# write landmark files
# write.csv(ventral, file.path(outputpath, "ventral_landmarks.csv"), row.names = FALSE)

#### Lateral

# write lateral landmarks.
# lateral <- rename(lateral, id = "file")
# lateral_landmarks <- rbind(lateral, lateral_missing)
# write.csv(lateral_landmarks, file.path(outputpath, "lateral_landmarks.csv"), row.names = FALSE)

## For manually marked images, ouput from ObjectJ

# the number of trait lengths to calculate
n.lateral.lengths <- 20 # 19 lateral lengths to calculate from marius, plus the standard length of the fish

lateral_missing_lengths <- get.lengths(n.lateral.lengths,lateral_missing)
lateral_lengths <- get.lengths(n.lateral.lengths,lateral)

all_lateral_lengths <- rbind(lateral_lengths,lateral_missing_lengths)

# write.csv(all_lateral_lengths, file.path(outputpath,"lateral_lengths.csv"), row.names = FALSE)

all_lateral_lengths <- read.csv(file.path(outputpath,"lateral_lengths.csv"))

## For ML-morph's testing photos
## ERROR: Error in data.frame(..., check.names = FALSE) :
## arguments imply differing number of rows: 1287, 0
# This error is probably bc the mlmorph file doesn't have the
# box_width, etc data. Not sure why that isn't there but it's not important
mlmorph_lengths <- get.lengths(n.lateral.lengths,mlmorph,mlmorph = TRUE)

# add scalebar measurements
scalebar <- all_lateral_lengths %>% select(id,scale.15.mm)

mlmorph_lengths <- left_join(mlmorph_lengths,scalebar,by = c("id"))

#write.csv(mlmorph_lengths, file.path(outputpath,"mlmorph_lateral_lengths.csv"), row.names = FALSE)
