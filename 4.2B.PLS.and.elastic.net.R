# === general ==================================================================

# Pretty sure this is old, the saved RDS files in the _figures R file 
# of the same name has whats in my thesis.


## 4.2B.PLS.R 
## perform a two block partial least squares analysis between morpho and enviro data
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.1.2 (2021-11-01)"
#

# === table of contents ==================================================================

# 1. load libraries and import data
#
# 2. create morpho and enviro blocks
#
# 3. perform 2B-PLS and save results and plots
#
# 4. Elastic net regression
#
#    - All terms togther in the same model
#
#    - Separate analysis for trophic, morpho, armor, and parasite

## There is a fair amount of commented out code here, this is from earlier analysis and can be ignored. I kept it here in case its needed again.

# ===  1. load libraries and import data ==================================================================
library(dplyr)
library(geomorph)
library(glmnet)


# tk only lateral length data, none of stuarts, so its commented out
# lateral_lengths <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/size_corrected_lateral_lengths_mm.csv")
# lateral_lengths$lake_id <- substring(lateral_lengths$id,1,nchar(lateral_lengths$id)-6)
# lateral_lengths$tag_prefix <- gsub('.{9}$', '', lateral_lengths$lake_id)
# lateral_lengths$WTRBDKGRPC <- substring(lateral_lengths$lake_id,nchar(lateral_lengths$lake_id)-8,nchar(lateral_lengths$lake_id))

# traits <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_kosciuch_stuart_lake_mean_parasitenotsc.csv")
# Until I add parasite data, use this dataset: These traits are size corrected.
traits <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_mean_size_corrected_trait_lengths.csv")

# mean per lake, not size corrected, suffix is _m
traits_m <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_mean_trait_lengths.csv")

#  individual traits sized corrected, suffix is _i.sc
traits_i.sc <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_size_corrected_trait_lengths.csv")

#  individual traits not sized corrected, suffix is _i
traits_i <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_trait_lengths.csv")


enviro <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/predictors.csv")

merge_yoel_tristan <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/link_yoel_and_marius_traits.csv") %>% 
  mutate(Yoels.trait.sc = paste(Yoels.trait,".sc",sep = ""),
         Marius.trait.letter.sc = paste(Marius.trait.letter,".sc",sep = ""))

# WTRBDPLD added to trait data so lake_joining_info no longer needed
# lake_joining_info <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/notebook_lakes.csv") %>% 
#   select(tag_prefix,WTRBDKGRPC,WTRBDPLD) %>% 
#   mutate(tag_prefix = toupper(tag_prefix)) %>% 
#   distinct()

# ===  2. create morpho and enviro blocks ==================================================================

# block <- left_join(lateral_lengths,lake_joining_info)
# block <- left_join(block,enviro)

prep_block <- function(enviro, traits, individual = F){

block <- left_join(traits, enviro, by = "WTRBDPLD")
lake_id <- traits$WTRBDPLD
# block_m <- left_join(enviro, traits_m, by = "WTRBDPLD")
# block_i.sc <- left_join(enviro, traits_i.sc, by = "WTRBDPLD")
# block_i <- left_join(enviro, traits_i, by = "WTRBDPLD")


# drop lakes with only two samples
# block <- block %>% filter(!lake_id %in% c("MH00385SALM","PY00098SALM")) # lake_id no longer is called watershed in the merged trait data
# block <- block %>% filter(!watershed %in% c("MH00385SALM","PY00098SALM"))

# extract environmental and morphomotric data as separate blocks
# morpho_cols <- c("dorsal.spine.1.len.mm.lateral.sc", "dorsal.spine.2.len.mm.lateral.sc",
#                  "dorsal.fin.len.mm.lateral.sc",             "caudal.depth.mm.lateral.sc",              
#                  "anal.fin.len.mm.lateral.sc",               "body.depth.mm.lateral.sc",                
#                  "mouth.len.mm.lateral.sc",                  "snout.len.mm.lateral.sc",                 
#                  "eye.diam.mm.lateral.sc",                   "head.len.mm.lateral.sc",                  
#                  "buccal.cavity.length.mm.sc",               "gape.width.mm.sc",                        
#                  "body.width.eye.mm.sc",                     "body.width.midbody.mm.sc",                
#                  "pelvic.girdle.width.mm.sc",                "pelvic.girdle.diamond.width.mm.sc",       
#                  "pelvic.girdle.length.mm.sc",               "pelvic.girdle.diamond.length.mm.sc",      
#                  "body.width.anal.1.mm.sc",                  "body.width.anal.2.mm.sc",                 
#                  "Left.Side.Pelvic.Spine.Length.mm.sc",      "Left.Side.Plate.Count.sc",                
#                  "Right.Side.Pelvic.Spine.Length.mm.sc",     "Right.side.gill.raker.number.insitu.sc",  
#                  "Thersitina.gill.copepod.0.1",              "Unionidae.gill.clam.0.1",                 
#                  "Schistocephalus.cestode.bodycavity.COUNT", "Eustrongylides.nematode.bodycavity.COUNT")

if(individual){
  morpho_cols <- traits %>% select(-c(id, WTRBDPLD, watershed)) %>% names()
  n_samples <- rep(1, dim(traits)[1])
} else {
  morpho_cols <- traits %>% select(-c(WTRBDPLD, watershed, n)) %>% names()
  n_samples <- block[,"n"]
}

enviro_cols <- enviro %>% select(-WTRBDPLD) %>% names()

# I join then subset to ensure the rows are in the same order
block_enviro <- data.frame(block[,enviro_cols])
block_morpho <- data.frame(block[,morpho_cols])

# PCA of monthly lake precipitation, keep the first PC as it explains 0.979 of the variance
precip_cols <- c("January.lake.precipitation..mm.","February.lake.precipitation..mm.",
                 "March.lake.precipitation..mm.","April.lake.precipitation..mm.",
                 "May.lake.precipitation..mm.","June.lake.precipitation..mm.",
                 "July.lake.precipitation..mm.","August.lake.precipitation..mm.",
                 "September.lake.precipitation..mm.","October.lake.precipitation..mm.",
                 "November.lake.precipitation..mm.","December.lake.precipitation..mm.")
precip <- block_enviro[,colnames(block_enviro) %in% precip_cols] %>% distinct()
precip_pca <- prcomp(precip)
summary(precip_pca)

precip <- cbind(precip,precip_pca$x[,"PC1"])
precip <- precip %>% mutate(precip_PC1 = `precip_pca$x[, "PC1"]`) %>% select(-c(`precip_pca$x[, "PC1"]`))

# PCA of monthly watershed basin precipitation, keep the first PC1 as it explains >0.999 of the variance
basin_precip_cols <- c("January.basin.precipitation..m.3.","February.basin.precipitation..m.3.","March.basin.precipitation..m.3.",
                       "April.basin.precipitation..m.3.","May.basin.precipitation..m.3.","June.basin.precipitation..m.3.",
                       "July.basin.precipitation..m.3.","August.basin.precipitation..m.3.","September.basin.precipitation..m.3.",
                       "October.basin.precipitation..m.3.","November.basin.precipitation..m.3.","December.basin.precipitation..m.3.")

basin_precip <- block_enviro[,colnames(block_enviro) %in% basin_precip_cols] %>% distinct()
basin_precip_pca <- prcomp(basin_precip)
summary(basin_precip_pca)

basin_precip <- cbind(basin_precip,basin_precip_pca$x[,"PC1"])
basin_precip <- basin_precip %>% mutate(basin_precip_PC1 = `basin_precip_pca$x[, "PC1"]`) %>% select(-c(`basin_precip_pca$x[, "PC1"]`))

# PCA of monthly temperature, keep first two PCs to get to >.95 variance explained
temperature_cols <- c("January.max.temp","February.max.temp","March.max.temp","April.max.temp",                     
  "May.max.temp","June.max.temp","July.max.temp","August.max.temp",
  "September.max.temp","October.max.temp","November.max.temp","December.max.temp",                  
  "January.min.temp","February.min.temp","March.min.temp","April.min.temp",                     
  "May.min.temp","June.min.temp","July.min.temp","August.min.temp",                    
  "September.min.temp","October.min.temp","November.min.temp","December.min.temp")
temperature <- block_enviro[,colnames(block_enviro) %in% temperature_cols] %>% distinct()
temperature_pca <- prcomp(temperature)
summary(temperature_pca)

temperature <- cbind(temperature,temperature_pca$x[,c("PC1","PC2")])
temperature <- temperature %>% mutate(temperature_PC1 = PC1,temperature_PC2 = PC2) %>% select(-c("PC1","PC2"))

# add the PCs to the enviro_block then drop the original precipitation and temperature columns
block_enviro <- left_join(block_enviro,temperature) %>% left_join(precip) %>% left_join(basin_precip) %>% select(-c(precip_cols,temperature_cols,basin_precip_cols))

# inspect data for how many orders of magnitude it spans
apply(block_morpho, 2, range) # Left.Side.Plate.Count.sc spans 1 order of magnitude
apply(block_enviro, 2, range) # basin area and paved and unpaved road length all span >2 orders of magnitude from the smallest to the largest values

# inspect the distribution of basin area, paved and unpaved road length, and left side plate count
hist(block_enviro$basin.area..m.2.)
hist(block_enviro$paved_road_length.m.)
hist(block_enviro$unpaved_road_length.m.)
# hist(block_morpho$armor.plate.count.sc)

# these variables all have a few observations with very high values, so log transform the variables to bring them closer together
# The minimum value of paved_road_length is 0, so add 1 to avoid errors
block_enviro$basin.area..m.2. <- log(block_enviro$basin.area..m.2.)
block_enviro$paved_road_length.m. <- log(block_enviro$paved_road_length.m. + 1)
block_enviro$unpaved_road_length.m. <- log(block_enviro$unpaved_road_length.m.)
# block_morpho$armor.plate.count.sc <- log(block_morpho$armor.plate.count.sc) # I don't think counts should be log transformed?

# check the distributions after log transformation, its much better
hist(block_enviro$basin.area..m.2.)
hist(block_enviro$paved_road_length.m.)
hist(block_enviro$unpaved_road_length.m.)
# hist(block_morpho$armor.plate.count.sc)



# scale and center blocks
block_enviro_scaled <- scale(block_enviro)
block_morpho_scaled <- scale(block_morpho)

return(list(block_enviro, block_morpho, block_enviro_scaled, block_morpho_scaled, n_samples, lake_id))
}

block_i.sc <- prep_block(enviro, traits_i.sc, individual = T)
block_i <- prep_block(enviro, traits_i, individual = T)
block_m <- prep_block(enviro, traits_m)
block <- prep_block(enviro, traits)

# ===  3. perform 2B-PLS and save results and plots ==================================================================

result <- two.b.pls(block_morpho_scaled,block_enviro_scaled)
result <- two.b.pls(block[[4]],block[[3]])
result_i.sc <- two.b.pls(block_i.sc[[4]],block_i.sc[[3]])
result_i <- two.b.pls(block_i[[4]],block_i[[3]])
result_m <- two.b.pls(block_m[[4]],block_m[[3]])


plot(result)
summary(result)
plot(result_i)
summary(result_i)
plot(result_i.sc)
summary(result_i.sc)
plot(result_m)
summary(result_m)
result$right.pls.vectors
result$left.pls.vectors

# ===  4. Elastic net regression ==================================================================

#### setup variables ####

armor_traits <- merge_yoel_tristan %>% filter(Yoels.trait.sc %in% c("dorsal.spine.1.len.mm.lateral.sc","dorsal.spine.2.len.mm.lateral.sc",
                                                                    "Left.Side.Pelvic.Spine.Length.mm.sc","Right.Side.Pelvic.Spine.Length.mm.sc",
                                                                    "Left.Side.Plate.Count.sc","pelvic.girdle.width.mm.sc",
                                                                    "pelvic.girdle.diamond.width.mm.sc","pelvic.girdle.length.mm.sc",
                                                                    "pelvic.girdle.diamond.length.mm.sc")) %>% 
  select(Marius.trait.letter, Marius.trait.letter.sc)

# eye diameter, snout length, head length, mouth length are not defined as trophic traits in Stuart 2017 NEE but maybe are trophic?
trophic_traits <- merge_yoel_tristan %>% filter(Yoels.trait.sc %in% c("buccal.cavity.length.mm.sc","gape.width.mm.sc","Right.side.gill.raker.number.insitu.sc")) %>% 
  select(Marius.trait.letter, Marius.trait.letter.sc) %>% 
  mutate(Marius.trait.letter.sc = case_when(Marius.trait.letter.sc == "right.side.gillraker.number.sc" ~ "right.side.gillraker.number",
                                            T ~ Marius.trait.letter.sc))

# missing standard length
swimming_traits <- merge_yoel_tristan %>% filter(Yoels.trait.sc %in% c("dorsal.fin.len.mm.lateral.sc","caudal.depth.mm.lateral.sc",
                                                                       "anal.fin.len.mm.lateral.sc","body.depth.mm.lateral.sc",
                                                                       "body.width.eye.mm.sc","body.width.midbody.mm.sc",
                                                                       "body.width.anal.1.mm.sc","body.width.anal.2.mm.sc")) %>% 
  select(Marius.trait.letter, Marius.trait.letter.sc)




#### allometric scaled blocks ####

# model suffix is _ub (unscaled blocks)
model_ub <- cv.glmnet(as.matrix(block[[1]]), as.matrix(block[[2]]), family = "mgaussian", nfolds = 5, relax = TRUE, gamma = c(0.5), weights = block[[5]])
coef(model_ub, s = "lambda.min")
coef(model_ub, s = "lambda.1se")
model_ub
plot(model_ub)

# do many values of lambda near the optimal value from cv.glmnet to see if accuracy can be improved
model_ub_fine <- cv.glmnet(as.matrix(block[[1]]), as.matrix(block[[2]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5),
                                   lambda = seq(0.1,model_ub$relaxed$lambda.1se+0.1,0.01), weights = block[[5]])
model_ub_fine
coef(model_ub_fine, s = "lambda.min")
coef(model_ub_fine, s = "lambda.1se")
plot(model_ub_fine)

trophic_ub <- cv.glmnet(as.matrix(block[[1]]),
                         as.matrix(block[[2]][,trophic_traits$Marius.trait.letter.sc]), # yoel trait names: 
                         nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])
trophic_ub
coef(trophic_ub, s = "lambda.min")
coef(trophic_ub, s = "lambda.1se")
plot(trophic_ub)


swimming_ub <- cv.glmnet(as.matrix(block[[1]]),
                          as.matrix(block[[2]][,swimming_traits$Marius.trait.letter.sc]),
                          nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])

swimming_ub
coef(swimming_ub, s = "lambda.min")
coef(swimming_ub, s = "lambda.1se")
plot(swimming_ub)


armor_ub <- cv.glmnet(as.matrix(block[[1]]),
                       as.matrix(block[[2]][,armor_traits$Marius.trait.letter.sc]),
                       nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights =  block[[5]])
armor_ub
coef(armor_ub, s = "lambda.min")
coef(armor_ub, s = "lambda.1se")
plot(armor_ub)


#### mean 0, sd 1 allometric scaled blocks ####

# model suffix is _sd1

model_sd1 <- cv.glmnet(as.matrix(block[[3]]), as.matrix(block[[4]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])
coef(model_sd1, s = model_sd1$relaxed$lambda.min)
model_sd1
plot(model_sd1)

# do many values of lambda near the optimal value from cv.glmnet to see if accuracy can be improved
model_sd1_scaled_fine <- cv.glmnet(as.matrix(block[[3]]), as.matrix(block[[4]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5),
                               lambda = seq(model_sd1$relaxed$lambda.min-.5,model_sd1$relaxed$lambda.1se,0.01), weights = block[[5]])
model_sd1_scaled_fine
coef(model_sd1_scaled_fine, s = "lambda.min")
coef(model_sd1_scaled_fine, s = "lambda.1se")
plot(model_sd1_scaled_fine)

trophic_sd1 <- cv.glmnet(as.matrix(block[[3]]),
                           as.matrix(block[[4]][,trophic_traits$Marius.trait.letter.sc]), # yoel trait names: 
                           nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])
trophic_sd1
coef(trophic_sd1, s = "lambda.min")
coef(trophic_sd1, s = "lambda.1se")
plot(trophic_sd1)


swimming_sd1 <- cv.glmnet(as.matrix(block[[3]]),
                            as.matrix(block[[4]][,swimming_traits$Marius.trait.letter.sc]),
                            nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])

swimming_sd1
coef(swimming_sd1, s = "lambda.min")
coef(swimming_sd1, s = "lambda.1se")
plot(swimming_sd1)


armor_sd1 <- cv.glmnet(as.matrix(block[[3]]),
                         as.matrix(block[[4]][,armor_traits$Marius.trait.letter.sc]),
                         nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])
armor_sd1
coef(armor_sd1, s = "lambda.min")
coef(armor_sd1, s = "lambda.1se")
plot(armor_sd1)


#### no allometric scaling blocks ####

# model suffix is _m 
model_m <- cv.glmnet(as.matrix(block_m[[1]]), as.matrix(block_m[[2]]), family = "mgaussian", nfolds = 5, relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
coef(model_m, s = "lambda.min")
coef(model_m, s = "lambda.1se")
model_m
plot(model_m)

# do many values of lambda near the optimal value from cv.glmnet to see if accuracy can be improved
model_m_fine <- cv.glmnet(as.matrix(block_m[[1]]), as.matrix(block_m[[2]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5),
                           lambda = seq(0.1,model_m$relaxed$lambda.1se+0.1,0.01), weights = block_m[[5]])
model_m_fine
coef(model_m_fine, s = "lambda.min")
coef(model_m_fine, s = "lambda.1se")
plot(model_m_fine)

trophic_m <- cv.glmnet(as.matrix(block_m[[1]]),
                        as.matrix(block_m[[2]][,trophic_traits$Marius.trait.letter]), # yoel trait names: 
                        nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
trophic_m
coef(trophic_m, s = "lambda.min")
coef(trophic_m, s = "lambda.1se")
plot(trophic_m)


swimming_m <- cv.glmnet(as.matrix(block_m[[1]]),
                         as.matrix(block_m[[2]][,swimming_traits$Marius.trait.letter]),
                         nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])

swimming_m
coef(swimming_m, s = "lambda.min")
coef(swimming_m, s = "lambda.1se")
plot(swimming_m)


armor_m <- cv.glmnet(as.matrix(block_m[[1]]),
                      as.matrix(block_m[[2]][,armor_traits$Marius.trait.letter]),
                      nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
armor_m
coef(armor_m, s = "lambda.min")
coef(armor_m, s = "lambda.1se")
plot(armor_m)


#### no allometric scaling mean 0, sd 1 blocks ####

# model suffix is _m.sd 
model_m.sd <- cv.glmnet(as.matrix(block_m[[3]]), as.matrix(block_m[[4]]), family = "mgaussian", nfolds = 5, relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
coef(model_m.sd, s = "lambda.min")
coef(model_m.sd, s = "lambda.1se")
model_m.sd
plot(model_m.sd)

# do many values of lambda near the optimal value from cv.glmnet to see if accuracy can be improved
model_m.sd_fine <- cv.glmnet(as.matrix(block_m[[3]]), as.matrix(block_m[[4]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5),
                          lambda = seq(0.1,model_m.sd$relaxed$lambda.1se+0.1,0.01), weights = block_m[[5]])
model_m.sd_fine
coef(model_m.sd_fine, s = "lambda.min")
coef(model_m.sd_fine, s = "lambda.1se")
plot(model_m.sd_fine)

trophic_m.sd <- cv.glmnet(as.matrix(block_m[[3]]),
                       as.matrix(block_m[[4]][,trophic_traits$Marius.trait.letter]), # yoel trait names: 
                       nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
trophic_m.sd
coef(trophic_m.sd, s = "lambda.min")
coef(trophic_m.sd, s = "lambda.1se")
plot(trophic_m.sd)


swimming_m.sd <- cv.glmnet(as.matrix(block_m[[3]]),
                        as.matrix(block_m[[4]][,swimming_traits$Marius.trait.letter]),
                        nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])

swimming_m.sd
coef(swimming_m.sd, s = "lambda.min")
coef(swimming_m.sd, s = "lambda.1se")
plot(swimming_m.sd)


armor_m.sd <- cv.glmnet(as.matrix(block_m[[3]]),
                     as.matrix(block_m[[4]][,armor_traits$Marius.trait.letter]),
                     nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
armor_m.sd
coef(armor_m.sd, s = "lambda.min")
coef(armor_m.sd, s = "lambda.1se")
plot(armor_m.sd)

#### individual values, SL as a predictor, unscaled ####

# for individual values, use foldid in cv.glmnet, make folds based on lake of origin
# Use match to get the integer values
integer_vector <- match(block_i[[6]], unique(block_i[[6]]))

# model suffix is _i. # Having issues with weights being used for some reason.
# make predictor and response by moving SL.mm. Ithink NAs are the issue

i_sl <- data.frame(block_i[[4]]) %>% 
  select(SL.mm)
i_response <- data.frame(block_i[[4]]) %>% 
  select(-SL.mm) %>% 
  as.matrix()
i_predictor <- as.matrix(data.frame(block_i[[3]], i_sl))
i_predictor <- as.matrix(data.frame(block_i[[3]]))


# drop NAs
i_predictor <- i_predictor[rowSums(is.na(i_response)) <= 0,]
i_response <- i_response[rowSums(is.na(i_response)) <= 0,]

model_i <- cv.glmnet(i_predictor, i_response, family = "mgaussian", relax = TRUE, gamma = c(0.5),
                     nfolds = max(nfoldid = ceiling(integer_vector/11)) , nfoldid = ceiling(integer_vector/11)) # I can use the ceiling() to change how many lakes per fold. Results dont really change much though
#, family = "mgaussian", nfolds = 5, relax = TRUE, gamma = c(0.5))
coef(model_i, s = "lambda.min")
coef(model_i, s = "lambda.1se")
model_i
plot(model_i)


model_i <- cv.glmnet(as.matrix(block_i[[1]]), as.matrix(block_i[[2]]), family = "mgaussian", relax = TRUE, gamma = c(0.5),
                     nfolds = max(nfoldid = ceiling(integer_vector/11)) , nfoldid = ceiling(integer_vector/11)) 

# # do many values of lambda near the optimal value from cv.glmnet to see if accuracy can be improved
# model_i_fine <- cv.glmnet(as.matrix(block_i[[3]]), as.matrix(block_i[[4]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5),
#                              lambda = seq(0.1,model_i$relaxed$lambda.1se+0.1,0.01), weights = block_i[[5]])
# model_i_fine
# coef(model_i_fine, s = "lambda.min")
# coef(model_i_fine, s = "lambda.1se")
# plot(model_i_fine)

trophic_i <- cv.glmnet(i_predictor, i_response[,trophic_traits$Marius.trait.letter], # yoel trait names: 
                          nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5))
trophic_i
coef(trophic_i, s = "lambda.min")
coef(trophic_i, s = "lambda.1se")
plot(trophic_i)


swimming_i <- cv.glmnet(i_predictor, i_response[,swimming_traits$Marius.trait.letter],
                           nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5))

swimming_i
coef(swimming_i, s = "lambda.min")
coef(swimming_i, s = "lambda.1se")
plot(swimming_i)


armor_i <- cv.glmnet(i_predictor, i_response[,armor_traits$Marius.trait.letter],
                        nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5))
armor_i
coef(armor_i, s = "lambda.min")
coef(armor_i, s = "lambda.1se")
plot(armor_i)

#### individual values, SL as a predictor, mean 0, sd 1 ####

#### interactions, best mean value method ####
# might be interesting to try glmnet with interaction terms: library(glinternet), tutorial: https://strakaps.github.io/post/glinternet/#fnref1

library(glinternet)


t <- glinternet.cv(i_predictor, i_response[,1], numLevels = rep(1, ncol(i_predictor)))
t2 <- glinternet.cv(as.matrix(block[[1]]), as.vector(block[[2]][,1]), numLevels = rep(1, ncol(block[[1]])))
plot(t)
plot(t2)

#### interactions, best individual values method ####






# # scaling just morpho trait is still not great, not converging on good predictions
# model_scaled_morpho <- cv.glmnet(as.matrix(block_enviro), as.matrix(block_morpho_scaled), family = "mgaussian", nfolds = dim(block)[1], relax = TRUE, gamma = c(0.5))
# coef(model_scaled_morpho, s = model_scaled_morpho$lambda.min)

# both scaled



# not ready to do parasite yet
# parasite_sd1 <- cv.glmnet(as.matrix(block_enviro_scaled), as.matrix(block_morpho_scaled[,c("Thersitina.gill.copepod.0.1","Unionidae.gill.clam.0.1",
#                                                                                              "Schistocephalus.cestode.bodycavity.COUNT","Eustrongylides.nematode.bodycavity.COUNT")]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5))
# parasite_sd1
# coef(parasite_sd1, s = parasite_sd1$relaxed$lambda.min)
# coef(parasite_sd1, s = parasite_sd1$relaxed$lambda.1se)
# plot(parasite_sd1)


#### Debugging cv.glmnet difference between lambda.min and what is shown when the model is called directly
# pretty sure its because multiple values are calculated in each fold and when called the single model trained on all the data is used
a <- runif(100)
b <- runif(100)
c <- runif(100)
d <- b + runif(100)/10
e <- a + runif(100)/10

test <- cv.glmnet(cbind(a,b,c), cbind(d,e), family = "mgaussian", relax = TRUE, gamma = 0.5)
test
test$gamma
