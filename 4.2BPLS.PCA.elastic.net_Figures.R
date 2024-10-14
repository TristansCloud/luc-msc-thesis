# === general ==================================================================

# This is the analysis used in my thesis

## 4.2BPLS.PCA.elastic.net.R 
## perform a two block partial least squares analysis between morpho and enviro data
## PCA and pearson correlation coefficient with first axis of PCA and 2BPLS trait ordination
## Then elastic net regression enviro ~ trait
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.3.1 (2021-11-01)"
#

# === table of contents ==================================================================

# 1. load libraries and import data
#
# 2. create morpho and enviro blocks
#
# 3. perform 2B-PLS and save results and plots
#
# 4. PCA and pearson correlation with 2BPLS
#
# 5. Elastic net regression
#
#    - All terms togther in the same model
#
#    - Separate analysis for trophic, morpho, armor, and parasite

## There is a fair amount of commented out code here, this is from earlier analysis and can be ignored. I kept it here in case its needed again.

# ===  1. load libraries and import data ==================================================================
library(dplyr)
library(geomorph)
library(glmnet)
library(ggfortify)
library(geomorph)
library(abind)
library(flextable)
library(foreach)
library(doParallel)
library(tidyr)

# tk only lateral length data, none of stuarts, so its commented out
# lateral_lengths <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/size_corrected_lateral_lengths_mm.csv")
# lateral_lengths$lake_id <- substring(lateral_lengths$id,1,nchar(lateral_lengths$id)-6)
# lateral_lengths$tag_prefix <- gsub('.{9}$', '', lateral_lengths$lake_id)
# lateral_lengths$WTRBDKGRPC <- substring(lateral_lengths$lake_id,nchar(lateral_lengths$lake_id)-8,nchar(lateral_lengths$lake_id))

# traits <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_kosciuch_stuart_lake_mean_parasitenotsc.csv")
# Until I add parasite data, use this dataset: These traits are size corrected.
traits <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_mean_size_corrected_trait_lengths.csv")



## this is code for other "versions" of data I could test. I did not do much investigating how 
## the results differ from the size corrected traits I used because of time constraints.

# # mean per lake, not size corrected, suffix is _m
# traits_m <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_mean_trait_lengths.csv")
# 
# #  individual traits sized corrected, suffix is _i.sc
# traits_i.sc <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_size_corrected_trait_lengths.csv")
# 
# #  individual traits not sized corrected, suffix is _i
# traits_i <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/kosciuch_stuart_trait_lengths.csv")


enviro <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/predictors.csv")

enviro_boots <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/on_the_ground_env.csv")

kosciuch_substrate <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Summer 2020/Summer 2020 datasheets/2020 tkosciuch fieldwork notebook data/2020 tkosciuch habitat categories.csv") %>% 
  mutate(substrate = paste(tolower(symbol), "substrate", sep = ".")) %>% 
  filter(substrate %in% names(enviro_boots))


rename_vec <- setNames(kosciuch_substrate$substrate, kosciuch_substrate$category)

# Rename columns using rename() function
enviro_boots <- enviro_boots %>% rename(!!!rename_vec)

enviro_both <- inner_join(enviro, enviro_boots)

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

prep_block <- function(enviro, traits, individual = F, boots_only = F){

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

if(!boots_only){
  
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
  
  # PCA of monthly temperature, keep first two PCs to get to >.95 variance explained.
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
  
}

# scale and center blocks
block_enviro_scaled <- scale(block_enviro)
block_morpho_scaled <- scale(block_morpho)

return(list(block_enviro, block_morpho, block_enviro_scaled, block_morpho_scaled, n_samples, lake_id))
}

# block_i.sc <- prep_block(enviro, traits_i.sc, individual = T)
# block_i <- prep_block(enviro, traits_i, individual = T)
# block_m <- prep_block(enviro, traits_m)
block <- prep_block(enviro, traits)
block_both <- prep_block(enviro_both, traits)
block_boots <- prep_block(enviro_boots, traits, boots_only = T)

# ===  3. perform 2B-PLS and save results and plots ==================================================================

# result <- two.b.pls(block_morpho_scaled,block_enviro_scaled)
result <- two.b.pls(block[[4]],block[[3]])
result_both <- two.b.pls(block_both[[4]],block_both[[3]])
result_boots <- two.b.pls(block_boots[[4]],block_boots[[3]])
# result_not_0 <- two.b.pls(block[[1]],block[[2]])
# result_i.sc <- two.b.pls(block_i.sc[[4]],block_i.sc[[3]])
# result_i <- two.b.pls(block_i[[4]],block_i[[3]])
# result_m <- two.b.pls(block_m[[4]],block_m[[3]])


plot(result)
summary(result)
plot(result_both)
summary(result_both)
plot(result_boots)
summary(result_boots)
# plot(result_not_0) # ew
# summary(result_not_0) # lower R^2, crazy thats got a .3 pearson correlation
# plot(result_i)
# summary(result_i)
# plot(result_i.sc)
# summary(result_i.sc)
# plot(result_m)
# summary(result_m)
result$right.pls.vectors
result$left.pls.vectors
result_both$right.pls.vectors
result_both$left.pls.vectors
result_boots$right.pls.vectors
result_boots$left.pls.vectors

# ===  4. PCA and correlation with 2BPLS ==================================================================

PCA <- prcomp(block[[2]], scale. = TRUE, center = TRUE)
summary(PCA)
# plot first two principal components and outline by lake, saved as a png
# png(file="C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/PCA.png",
#     width=600, height=350)
autoplot(PCA, data = traits)
# dev.off()

PCA$rotation

# correlation, pretty high correlation but 2BPLS is still trying to maximize the proportion of variance in the first axis
cor(result$left.pls.vectors[,1], PCA$rotation[,1], method = "pearson")
# try with the env coord
env_PCA <- prcomp(block[[1]], scale. = TRUE, center = TRUE)
summary(env_PCA)
cor(result$right.pls.vectors[,1], env_PCA$rotation[,1], method = "pearson")

# ===  5. Elastic net regression ==================================================================

#### setup variables ####

armor_traits <- merge_yoel_tristan %>% filter(Yoels.trait.sc %in% c("dorsal.spine.1.len.mm.lateral.sc","dorsal.spine.2.len.mm.lateral.sc",
                                                                    "Left.Side.Pelvic.Spine.Length.mm.sc","Right.Side.Pelvic.Spine.Length.mm.sc",
                                                                    "Left.Side.Plate.Count.sc","pelvic.girdle.width.mm.sc",
                                                                    "pelvic.girdle.diamond.width.mm.sc","pelvic.girdle.length.mm.sc",
                                                                    "pelvic.girdle.diamond.length.mm.sc")) %>% 
  select(Marius.trait.letter, Marius.trait.letter.sc)

# eye diameter, snout length, head length, mouth length are not defined as trophic traits in Stuart 2017 NEE but maybe are trophic?
# some traits from yoel NEE that are trophic and I might have: Opercular 4-bar diagonal length, Neurocranium length (not sure if yoel measured this with the mouth open). Not sure where this data is from Yoel.
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

# benthic/limnetic traits (caudal peduncle and body depth). 
bio_relevance_traits <-  merge_yoel_tristan %>% filter(Yoels.trait.sc %in% c("caudal.depth.mm.lateral.sc","body.depth.mm.lateral.sc")) %>% 
  select(Marius.trait.letter, Marius.trait.letter.sc)

bootstrap_glmnet <- function(predictors, responses, weights, traits = NULL, bootstraps = 50, core = 1) {
  output <- list()
  
  if(!is.null(traits)){
    responses <- responses %>% select(traits)
  }
  
  if(as.integer(core) == 1){
    
    for(i in 1:bootstraps){
      output[[i]] <- cv.glmnet(x = predictors, y = responses, family = "mgaussian", gamma = c(0.5), weights = weights, intercept = T, standardize.response = T, nfolds = 10) #, keep = TRUE
    }
    
  }
  
  if(as.integer(core) > 1){
    
    cl <- makeCluster(core, type = "PSOCK")
    registerDoParallel(cl)

    for(i in 1:bootstraps){
      output[[i]] <- cv.glmnet(x = predictors, y = responses, family = "mgaussian", gamma = c(0.5), weights = weights, intercept = T, standardize.response = T, nfolds = 10, parallel = T) #, keep = TRUE
    }

    # Stop parallel backend
    stopCluster(cl)
    
  }

  
  return(output)
}

summarize_bootstraps <- function(glmnet_list, lambda){
  fit <- coef(glmnet_list[[1]], s = "lambda.1se")
  
  # setup outputs
  parameters <- names(fit)
  output <- vector(mode = "list", length = length(fit))
  names(output) <- parameters

  for(i in 1:length(glmnet_list)){
    
    if(lambda == "lambda.1se") {
      fit <- coef(glmnet_list[[i]], s = "lambda.1se")
    }
    if(lambda == "lambda.min"){
      fit <- coef(glmnet_list[[i]], s = "lambda.min")
    }
    
    for(j in 1:length(fit)){
      output[[j]] <- bind_rows(output[[j]],fit[[j]][,])
    }
  }
  
  return(output)
}

summarize_range_bootstep <- function(glmnet_list){
  fit <- coef(glmnet_list[[1]], s = "lambda.1se")
  # 
  parameters <- names(fit)
  # output <- vector(mode = "list", length = length(fit))
  # names(output) <- parameters
  
  
  for(i in 1:length(glmnet_list)){
    
    lambda <- glmnet_list[[i]]$lambda[glmnet_list[[i]]$lambda.min <= glmnet_list[[i]]$lambda] #  & glmnet_list[[i]]$lambda <= glmnet_list[[i]]$lambda.1se # if I want to restrict it between min and 1se
    
    fit <- coef(glmnet_list[[i]], s = lambda)
    
    if(i == 1){
      fit.coefficients <- fit
    } else {
      
      for(j in parameters){
        
        fit.coefficients[[j]] <- cbind(fit.coefficients[[j]], fit[[j]])
        
      }
      
    }
    
  }
  
  # vote count predictors
  vote_count_proportion <- rowSums(abs(fit.coefficients[[j]]) > 0) / ncol(fit.coefficients[[j]]) # any value of j doesn't change the vote count so I left it as j
  
  output <- vector(mode = "list", length = length(parameters) + 1)
  names(output) <- c("vote count",parameters)
  
  output[["vote count"]] <- vote_count_proportion
  
  for(j in parameters){
    
    sum_stats <- data.frame(vector(mode = "character", length = dim(fit.coefficients[[1]])[1]),
                            vector(mode = "numeric", length = dim(fit.coefficients[[1]])[1]),
                            vector(mode = "numeric", length = dim(fit.coefficients[[1]])[1]),
                            vector(mode = "numeric", length = dim(fit.coefficients[[1]])[1]),
                            vector(mode = "numeric", length = dim(fit.coefficients[[1]])[1]),
                            vector(mode = "numeric", length = dim(fit.coefficients[[1]])[1]))
    
    sum_stats[,1] <- names(fit.coefficients[[1]][,1])
    colnames(sum_stats) <- c("variable", "mean", "median", "sd", "min", "max")
    
    sum_stats[,"mean"] <- rowMeans(fit.coefficients[[j]])
    sum_stats[,"median"] <- apply(fit.coefficients[[j]], 1, median)
    sum_stats[,"sd"] <- apply(fit.coefficients[[j]], 1, sd)
    sum_stats[,"min"] <- apply(fit.coefficients[[j]], 1, min)
    sum_stats[,"max"] <- apply(fit.coefficients[[j]], 1, max)
    
    output[[j]] <- sum_stats
    
  }
  
  return(output)
  
}

summary_tables <- function(bootstrap_summary) {
  # use flextable() to turn dataframes into visual tables
  
  # vote count predictors
  vote_count_proportion <- colSums(abs(bootstrap_summary[[1]]) > 0) / nrow(bootstrap_summary[[1]])
  
  # sumstats per trait, I need a dataframe per trait. Each dataframe will have the predictor variables as rows and the sumstats as columns
  # sum_stats <- vector(mode = "list", length = length(bootstrap_summary))
  sum_stats <- data.frame(vector(mode = "character", length = dim(bootstrap_summary[[1]])[2]),
                          vector(mode = "numeric", length = dim(bootstrap_summary[[1]])[2]),
                          vector(mode = "numeric", length = dim(bootstrap_summary[[1]])[2]),
                          vector(mode = "numeric", length = dim(bootstrap_summary[[1]])[2]),
                          vector(mode = "numeric", length = dim(bootstrap_summary[[1]])[2]),
                          vector(mode = "numeric", length = dim(bootstrap_summary[[1]])[2]))
  
  sum_stats[,1] <- names(bootstrap_summary[[1]])
  colnames(sum_stats) <- c("variable", "mean", "median", "sd", "min", "max")
  
  output <- list(vote_count_proportion)
  
  for(i in 1:length(bootstrap_summary)){
    sum_stats_i <- sum_stats
    sum_stats_i[,"mean"] <- colMeans(bootstrap_summary[[i]])
    sum_stats_i[,"median"] <- apply(bootstrap_summary[[i]], 2, median)
    sum_stats_i[,"sd"] <- apply(bootstrap_summary[[i]], 2, sd)
    sum_stats_i[,"min"] <- apply(bootstrap_summary[[i]], 2, min)
    sum_stats_i[,"max"] <- apply(bootstrap_summary[[i]], 2, max)
    # sum_stats_i <- round(sum_stats_i, 3)
    # sum_stats_i[,"mean"] <- round(colMeans(bootstrap_summary[[i]]), 3)
    # sum_stats_i[,"median"] <- round(apply(bootstrap_summary[[i]], 2, median), 3)
    # sum_stats_i[,"sd"] <- round(apply(bootstrap_summary[[i]], 2, sd), 3)
    # sum_stats_i[,"min"] <- round(apply(bootstrap_summary[[i]], 2, min), 3)
    # sum_stats_i[,"max"] <- round(apply(bootstrap_summary[[i]], 2, max), 3)
    output <- append(output,list(sum_stats_i))
  }
  names(output) <- c("vote count proportion", names(bootstrap_summary))
  return(output)
}

cv_error <- function(glmnet_list, lambda){ # return the cv MSE at lambda
  output <- c()
  
  for(i in 1:length(glmnet_list)){
    if(lambda == "lambda.1se") {
      index <- glmnet_list[[i]]$index["1se",]
      error <- glmnet_list[[i]]$cvm[index]
    }
    if(lambda == "lambda.min"){
      index <- glmnet_list[[i]]$index["min",]
      error <- glmnet_list[[i]]$cvm[index]
    }
    output <- append(output,error)
  }
  return(output)
}

cv_means <- function(traits, bootstraps = 50, num_folds = 10, method = "mean"){ # Calculate cv MSE across folds using trait means, bootstrap this
  
  # # mean squared Frobenius norm of the error:
  # # 1. Get elementwise difference between prediction and actual then square it(y_pred - y_actual)^2
  # # Assuming each row is an individual, take the row sum
  # # take the mean of all the row sums to get your answer
  # 
  # # proof that this function does calculate the mean squared Frobenius norm of the error
  # y <- data.frame(runif(100),runif(100),runif(100))
  # z <- colMeans(y +rnorm(100))
  # z_mat <- matrix(rep(z,100), ncol = 3, byrow = T)
  # 
  # #Trevor Hastie's way
  # hastie.errormin = (y - z_mat)^2
  # mean(rowSums(hastie.errormin))
  # 
  # # The way used in this function
  # func.errormin <- sweep(y, 2, z, FUN = "-")^2
  # mean(rowSums(func.errormin)) # It matches!
  
  # for each bootstrap, assign foldID (1 to 10) to each row.
  # calculate the mean then test on the fold, repeat 10 times, keep the MSE
  
  output <- c() # MSEs, should be 50 long
  
  for(i in 1:bootstraps){
    
    # create folds
    num_rows <- nrow(traits)
    balanced_counts <- rep(num_rows %/% num_folds, num_folds)
    remainder <- num_rows %% num_folds
    if(remainder > 0){
      balanced_counts[1:remainder] <- balanced_counts[1:remainder] + 1
    }
    fold <- rep(1:num_folds, balanced_counts)
    fold <- sample(fold)
    
    mse <- c()
    
    for(j in 1:num_folds){
      train <- traits[fold!=j,]
      test <- traits[fold==j,]
      errormin <- sweep(test, 2, colMeans(train), FUN = "-")^2
      mse <- append(mse,mean(rowSums(errormin)))
    }
    output <- append(output, mean(mse))
  }
  
  return(output)
}

# nice this works
round_with_zeros <- function(x, digits, as_character = T) {
  formatted_numbers <- sprintf(paste0("%.", digits, "g"), x)
  # formatted_numbers <- ifelse(abs(x) > 0.01, round(x, digits = digits), sprintf("%0.2g", x))
  if(as_character){
    return(as.character(formatted_numbers))
  } else {
    return(as.numeric(formatted_numbers))
  }
}

# # Test the function
# x <- c(1.002342, 0.4234523, 0.000023423)
# rounded_numbers <- round_with_zeros(x, 2)
# print(rounded_numbers)

pretty_tables <- function(table_list, directory = "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/tables", 
                          concatenate = F, save_file = T, digits = 2, keep_cols = NULL){
  
  # this gives the list object name
  # print(deparse(substitute(table_list)))
  
  list_name <- deparse(substitute(table_list))
  
  numeric_cols <- names(table_list[[2]])[names(table_list[[2]]) != "variable"]
  
  for (table_name in names(table_list)) {
    table <- table_list[[table_name]]
    for (i in seq_along(table)) {
      if (is.numeric(table[[i]]) | is.double(table[[i]])) {
        table[[i]] <- round_with_zeros(table[[i]], 2)  # Adjust the number of significant digits as needed
      }
    }
    table_list[[table_name]] <- table
    
    # print(table_name)
    
    # 
    # print(output_path)
    # print("")
    # 
  }
  
  if(save_file) {
    
    if(concatenate){
      
      output_path <- file.path(directory, paste(list_name, ".docx", sep = "_"))
      
      output <- data.frame(table_list)
      
      print(names(output))
      
      if(!isFALSE(keep_cols)) {
        
        output <- output[,grepl(keep_cols,names(output))]
        
      }
      
      output <- cbind(rownames(output), data.frame(output, row.names=NULL))
      
      names(output)[names(output) == 'rownames(output)'] <- 'Predictor'
      
      flextable(output) %>% 
        # colformat_double(digits = 2, scientific = TRUE) %>%
        # set_formatter(digits = 2, scientific = T)
        # set_formatter(numbers = function(x) {
        #   formatC(x, format = "e", digits = 3)
        # }) %>%
        save_as_docx(path = output_path, align = "left")
      # save_as_image(output_path)
      
      return(output)
      
    } else {
      
      for (table_name in names(table_list)[!names(table_list) %in% "vote count proportion"]) {
        output_path <- file.path(directory, paste(list_name, table_name, ".docx", sep = "_"))
        
        output <- table_list[[table_name]]
        
        # output <- cbind(rownames(output), data.frame(output, row.names=NULL))
        # 
        # names(output)[names(output) == 'rownames(output)'] <- 'Predictor'
        
        flextable(output) %>% 
          # colformat_double(digits = 2, scientific = TRUE) %>%
          # set_formatter(digits = 2, scientific = T)
          # set_formatter(numbers = function(x) {
          #   formatC(x, format = "e", digits = 3)
          # }) %>%
          save_as_docx(path = output_path, align = "left")
        # save_as_image(output_path)
      }
      return(table_list)
    }
      
    }
}


vote_count_table <- function(table_list, name_order, varname_from = NULL, varname_to = NULL){
  
  if(typeof(varname_from) != typeof(varname_to)){
    stop("Missing a varname to or from!")
  }
  
  
  for (i in seq_along(table_list)) {
    
    if(i == 1){
      
      output <- data.frame(table_list[[i]]$"vote count proportion")
      
    } else {
      
      output <- data.frame(output, table_list[[i]]$"vote count proportion")
      
    }
    
  }
  
  names(output) <- name_order
  
  if(typeof(varname_from) == "character") {
    output_names_i <- match(varname_from, rownames(output))
    
    rownames(output)[output_names_i] <- varname_to
  }
  
  return(output)
  
}

# t <- vote_count_table(list(table_min, table_trophic_min, table_swimming_min, table_armor_min),
#                       name_order = c("All traits", "Trophic", "Swimming", "Armor"))

# show parallel speedup, 6 cores gives ~50% speedup
# t <- Sys.time()
# bootstrap_glmnet(as.matrix(block_both[[1]]),as.matrix(block_both[[2]]), weights = block_both[[5]], core = 6)
# u <- Sys.time()
# bootstrap_glmnet(as.matrix(block_both[[1]]),as.matrix(block_both[[2]]), weights = block_both[[5]])
# v <- Sys.time()

#### rdata files of 500 bootstraps: ####

# saveRDS(bootstrap_both, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_both.rds")
# saveRDS(bootstrap_trophic_both, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_trophic_both.rds")
# saveRDS(bootstrap_swimming_both, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_swimming_both.rds")
# saveRDS(bootstrap_armor_both, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_armor_both.rds")
# 
# saveRDS(bootstrap_boots, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_boots.rds")
# saveRDS(bootstrap_trophic_boots, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_trophic_boots.rds")
# saveRDS(bootstrap_swimming_boots, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_swimming_boots.rds")
# saveRDS(bootstrap_armor_boots, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_armor_boots.rds")
# 
# saveRDS(bootstrap, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap.rds")
# saveRDS(bootstrap_trophic, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_trophic.rds")
# saveRDS(bootstrap_swimming, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_swimming.rds")
# saveRDS(bootstrap_armor, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_armor.rds")
# 
# saveRDS(all_cv_means, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/all_cv_means.rds")
# saveRDS(trophic_cv_means, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/trophic_cv_means.rds")
# saveRDS(swimming_cv_means, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/swimming_cv_means.rds")
# saveRDS(armor_cv_means, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/armor_cv_means.rds")

# read in saved data, the original functions to generate this data is commented out below

bootstrap_both <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_both.rds")
bootstrap_trophic_both <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_trophic_both.rds")
bootstrap_swimming_both <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_swimming_both.rds")
bootstrap_armor_both <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_armor_both.rds")

bootstrap_boots <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_boots.rds")
bootstrap_trophic_boots <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_trophic_boots.rds")
bootstrap_swimming_boots <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_swimming_boots.rds")
bootstrap_armor_boots <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_armor_boots.rds")

bootstrap <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap.rds")
bootstrap_trophic <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_trophic.rds")
bootstrap_swimming <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_swimming.rds")
bootstrap_armor <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/bootstrap_armor.rds")

all_cv_means <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/all_cv_means.rds")
trophic_cv_means <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/trophic_cv_means.rds")
swimming_cv_means <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/swimming_cv_means.rds")
armor_cv_means <- readRDS("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Thesis chapters/Exploratory chapter/rds files/armor_cv_means.rds")

bootstrap_bio_relevance_both <- readRDS()
bootstrap_bio_relevance_boots <- readRDS()
bootstrap_bio_relevance <- readRDS()


#### boots and GIS ####
# bootstrap_both <- bootstrap_glmnet(as.matrix(block_both[[1]]),as.matrix(block_both[[2]]), weights = block_both[[5]], core = 6, bootstraps = 500)
sum_min_both <- summarize_bootstraps(bootstrap_both, lambda = "lambda.min")
sum_1se_both <- summarize_bootstraps(bootstrap_both, lambda = "lambda.1se")
sum_range_both <- summarize_range_bootstep(bootstrap_both)
table_min_both <- summary_tables(sum_min_both)
table_1se_both <- summary_tables(sum_1se_both)
error_min_both <- cv_error(bootstrap_both, lambda = "lambda.min")
error_1se_both <- cv_error(bootstrap_both, lambda = "lambda.1se")
# all_cv_means_both <- cv_means(as.matrix(block_both[[2]]), bootstraps = 500) 
t.test(error_min_both, all_cv_means)
t.test(error_1se_both, all_cv_means)

# trophic
# bootstrap_trophic_both <- bootstrap_glmnet(as.matrix(block_both[[1]]),as.matrix(block_both[[2]][,trophic_traits$Marius.trait.letter.sc]), weights = block_both[[5]], core = 6, bootstraps = 500)
sum_trophic_min_both <- summarize_bootstraps(bootstrap_trophic_both, lambda = "lambda.min")
sum_trophic_1se_both <- summarize_bootstraps(bootstrap_trophic_both, lambda = "lambda.1se")
sum_range_trophic_both <- summarize_range_bootstep(bootstrap_trophic_both)
table_trophic_min_both <- summary_tables(sum_trophic_min_both)
table_trophic_1se_both <- summary_tables(sum_trophic_1se_both)
error_trophic_min_both <- cv_error(bootstrap_trophic_both, lambda = "lambda.min")
error_trophic_1se_both <- cv_error(bootstrap_trophic_both, lambda = "lambda.1se")
# trophic_cv_means_both <- cv_means(as.matrix(block_both[[2]][,trophic_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_trophic_min_both, trophic_cv_means)
t.test(error_trophic_1se_both, trophic_cv_means)

# swimming
# bootstrap_swimming_both <- bootstrap_glmnet(as.matrix(block_both[[1]]),as.matrix(block_both[[2]][,swimming_traits$Marius.trait.letter.sc]), weights = block_both[[5]], core = 6, bootstraps = 500)
sum_swimming_min_both <- summarize_bootstraps(bootstrap_swimming_both, lambda = "lambda.min")
sum_swimming_1se_both <- summarize_bootstraps(bootstrap_swimming_both, lambda = "lambda.1se")
sum_range_swimming_both <- summarize_range_bootstep(bootstrap_swimming_both)
table_swimming_min_both <- summary_tables(sum_swimming_min_both)
table_swimming_1se_both <- summary_tables(sum_swimming_1se_both)
error_swimming_min_both <- cv_error(bootstrap_swimming_both, lambda = "lambda.min")
error_swimming_1se_both <- cv_error(bootstrap_swimming_both, lambda = "lambda.1se")
# swimming_cv_means_both <- cv_means(as.matrix(block_both[[2]][,swimming_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_swimming_min_both, swimming_cv_means)
t.test(error_swimming_1se_both, swimming_cv_means)

# armor
# bootstrap_armor_both <- bootstrap_glmnet(as.matrix(block_both[[1]]),as.matrix(block_both[[2]][,armor_traits$Marius.trait.letter.sc]), weights = block_both[[5]], core = 6, bootstraps = 500)
sum_armor_min_both <- summarize_bootstraps(bootstrap_armor_both, lambda = "lambda.min")
sum_armor_1se_both <- summarize_bootstraps(bootstrap_armor_both, lambda = "lambda.1se")
sum_range_armor_both <- summarize_range_bootstep(bootstrap_armor_both)
table_armor_min_both <- summary_tables(sum_armor_min_both)
table_armor_1se_both <- summary_tables(sum_armor_1se_both)
error_armor_min_both <- cv_error(bootstrap_armor_both, lambda = "lambda.min")
error_armor_1se_both <- cv_error(bootstrap_armor_both, lambda = "lambda.1se")
# armor_cv_means_both <- cv_means(as.matrix(block_both[[2]][,armor_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_armor_min_both, armor_cv_means)
t.test(error_armor_1se_both, armor_cv_means)
# flextable(table_armor_min[[2]]) %>% save_as_docx( path = "C:/Users/trist/OneDrive/Desktop/test_flextable.docx")

# bio relevance
bootstrap_bio_relevance_both <- bootstrap_glmnet(as.matrix(block_both[[1]]),as.matrix(block_both[[2]][,bio_relevance_traits$Marius.trait.letter.sc]), weights = block_both[[5]], core = 6, bootstraps = 500)
sum_bio_relevance_min_both <- summarize_bootstraps(bootstrap_bio_relevance_both, lambda = "lambda.min")
sum_bio_relevance_1se_both <- summarize_bootstraps(bootstrap_bio_relevance_both, lambda = "lambda.1se")
sum_bio_relevance_range_both <- summarize_range_bootstep(bootstrap_bio_relevance_both)
table_bio_relevance_min_both <- summary_tables(sum_bio_relevance_min_both)
table_bio_relevance_1se_both <- summary_tables(sum_bio_relevance_1se_both)
error_bio_relevance_min_both <- cv_error(bootstrap_bio_relevance_both, lambda = "lambda.min")
error_bio_relevance_1se_both <- cv_error(bootstrap_bio_relevance_both, lambda = "lambda.1se")
bio_relevance_cv_means <- cv_means(as.matrix(block_boots[[2]][,bio_relevance_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_bio_relevance_min_both, bio_relevance_cv_means)
t.test(error_bio_relevance_1se_both, bio_relevance_cv_means)

vote_count_min_both <- vote_count_table(list(table_min_both, table_trophic_min_both, table_swimming_min_both, table_armor_min_both, table_bio_relevance_min_both),
                      name_order = c("All traits", "Trophic", "Swimming", "Armor", "Benthic - Limnetic"),
                      )

#### just boots on the ground enviro ####
# bootstrap_boots <- bootstrap_glmnet(as.matrix(block_boots[[1]]),as.matrix(block_boots[[2]]), weights = block_boots[[5]], core = 6, bootstraps = 500)
sum_min_boots <- summarize_bootstraps(bootstrap_boots, lambda = "lambda.min")
sum_1se_boots <- summarize_bootstraps(bootstrap_boots, lambda = "lambda.1se")
sum_range_boots <- summarize_range_bootstep(bootstrap_boots)
table_min_boots <- summary_tables(sum_min_boots)
table_1se_boots <- summary_tables(sum_1se_boots)
error_min_boots <- cv_error(bootstrap_boots, lambda = "lambda.min")
error_1se_boots <- cv_error(bootstrap_boots, lambda = "lambda.1se")
# all_cv_means_boots <- cv_means(as.matrix(block_boots[[2]]), bootstraps = 500) 
t.test(error_min_boots, all_cv_means)
t.test(error_1se_boots, all_cv_means)

# trophic
# bootstrap_trophic_boots <- bootstrap_glmnet(as.matrix(block_boots[[1]]),as.matrix(block_boots[[2]][,trophic_traits$Marius.trait.letter.sc]), weights = block_boots[[5]], core = 6, bootstraps = 500)
sum_trophic_min_boots <- summarize_bootstraps(bootstrap_trophic_boots, lambda = "lambda.min")
sum_trophic_1se_boots <- summarize_bootstraps(bootstrap_trophic_boots, lambda = "lambda.1se")
sum_range_trophic_boots <- summarize_range_bootstep(bootstrap_trophic_boots)
table_trophic_min_boots <- summary_tables(sum_trophic_min_boots)
table_trophic_1se_boots <- summary_tables(sum_trophic_1se_boots)
error_trophic_min_boots <- cv_error(bootstrap_trophic_boots, lambda = "lambda.min")
error_trophic_1se_boots <- cv_error(bootstrap_trophic_boots, lambda = "lambda.1se")
# trophic_cv_means_boots <- cv_means(as.matrix(block_boots[[2]][,trophic_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_trophic_min_boots, trophic_cv_means)
t.test(error_trophic_1se_boots, trophic_cv_means)

# swimming
# bootstrap_swimming_boots <- bootstrap_glmnet(as.matrix(block_boots[[1]]),as.matrix(block_boots[[2]][,swimming_traits$Marius.trait.letter.sc]), weights = block_boots[[5]], core = 6, bootstraps = 500)
sum_swimming_min_boots <- summarize_bootstraps(bootstrap_swimming_boots, lambda = "lambda.min")
sum_swimming_1se_boots <- summarize_bootstraps(bootstrap_swimming_boots, lambda = "lambda.1se")
sum_range_swimming_boots <- summarize_range_bootstep(bootstrap_swimming_boots)
table_swimming_min_boots <- summary_tables(sum_swimming_min_boots)
table_swimming_1se_boots <- summary_tables(sum_swimming_1se_boots)
error_swimming_min_boots <- cv_error(bootstrap_swimming_boots, lambda = "lambda.min")
error_swimming_1se_boots <- cv_error(bootstrap_swimming_boots, lambda = "lambda.1se")
# swimming_cv_means_boots <- cv_means(as.matrix(block_boots[[2]][,swimming_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_swimming_min_boots, swimming_cv_means)
t.test(error_swimming_1se_boots, swimming_cv_means)

# armor
# bootstrap_armor_boots <- bootstrap_glmnet(as.matrix(block_boots[[1]]),as.matrix(block_boots[[2]][,armor_traits$Marius.trait.letter.sc]), weights = block_boots[[5]], core = 6, bootstraps = 500)
sum_armor_min_boots <- summarize_bootstraps(bootstrap_armor_boots, lambda = "lambda.min")
sum_armor_1se_boots <- summarize_bootstraps(bootstrap_armor_boots, lambda = "lambda.1se")
sum_range_armor_boots <- summarize_range_bootstep(bootstrap_armor_boots)
table_armor_min_boots <- summary_tables(sum_armor_min_boots)
table_armor_1se_boots <- summary_tables(sum_armor_1se_boots)
error_armor_min_boots <- cv_error(bootstrap_armor_boots, lambda = "lambda.min")
error_armor_1se_boots <- cv_error(bootstrap_armor_boots, lambda = "lambda.1se")
# armor_cv_means_boots <- cv_means(as.matrix(block_boots[[2]][,armor_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_armor_min_boots, armor_cv_means)
t.test(error_armor_1se_boots, armor_cv_means)
# flextable(table_armor_min[[2]]) %>% save_as_docx( path = "C:/Users/trist/OneDrive/Desktop/test_flextable.docx")

# bio relevance
bootstrap_bio_relevance_boots <- bootstrap_glmnet(as.matrix(block_boots[[1]]),as.matrix(block_boots[[2]][,bio_relevance_traits$Marius.trait.letter.sc]), weights = block_boots[[5]], core = 6, bootstraps = 500)
sum_bio_relevance_min_boots <- summarize_bootstraps(bootstrap_bio_relevance_boots, lambda = "lambda.min")
sum_bio_relevance_1se_boots <- summarize_bootstraps(bootstrap_bio_relevance_boots, lambda = "lambda.1se")
sum_bio_relevance_range_boots <- summarize_range_bootstep(bootstrap_bio_relevance_boots)
table_bio_relevance_min_boots <- summary_tables(sum_bio_relevance_min_boots)
table_bio_relevance_1se_boots <- summary_tables(sum_bio_relevance_1se_boots)
error_bio_relevance_min_boots <- cv_error(bootstrap_bio_relevance_boots, lambda = "lambda.min")
error_bio_relevance_1se_boots <- cv_error(bootstrap_bio_relevance_boots, lambda = "lambda.1se")
# bio_relevance_cv_means <- cv_means(as.matrix(block_boots[[2]][,bio_relevance_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_bio_relevance_min_boots, bio_relevance_cv_means)
t.test(error_bio_relevance_1se_boots, bio_relevance_cv_means)

vote_count_min_boots <- vote_count_table(list(table_min_boots, table_trophic_min_boots, table_swimming_min_boots, table_armor_min_boots, table_bio_relevance_min_boots),
                                   name_order = c("All traits", "Trophic", "Swimming", "Armor", "Benthic - Limnetic"))


#### Just GIS (allometric scaled blocks, the same as both and boots) ####

# bootstrap <- bootstrap_glmnet(as.matrix(block[[1]]),as.matrix(block[[2]]), weights = block[[5]], core = 6, bootstraps = 500)
sum_min <- summarize_bootstraps(bootstrap, lambda = "lambda.min")
sum_1se <- summarize_bootstraps(bootstrap, lambda = "lambda.1se")
sum_range <- summarize_range_bootstep(bootstrap)
table_min <- summary_tables(sum_min)
table_1se <- summary_tables(sum_1se)
error_min <- cv_error(bootstrap, lambda = "lambda.min")
error_1se <- cv_error(bootstrap, lambda = "lambda.1se")
# all_cv_means <- cv_means(as.matrix(block[[2]]), bootstraps = 500)
t.test(error_min, all_cv_means)
t.test(error_1se, all_cv_means)
t.test(error_min, error_min_both)
t.test(error_min, error_min_boots)
t.test(error_min_both, error_min_boots)
# GIS has the lowest error, both is middle, and boots has the highest. All are lower than using means

# trophic
# bootstrap_trophic <- bootstrap_glmnet(as.matrix(block[[1]]),as.matrix(block[[2]][,trophic_traits$Marius.trait.letter.sc]), weights = block[[5]], core = 6, bootstraps = 500)
sum_trophic_min <- summarize_bootstraps(bootstrap_trophic, lambda = "lambda.min")
sum_trophic_1se <- summarize_bootstraps(bootstrap_trophic, lambda = "lambda.1se")
sum_range_trophic <- summarize_range_bootstep(bootstrap_trophic)
table_trophic_min <- summary_tables(sum_trophic_min)
table_trophic_1se <- summary_tables(sum_trophic_1se)
error_trophic_min <- cv_error(bootstrap_trophic, lambda = "lambda.min")
error_trophic_1se <- cv_error(bootstrap_trophic, lambda = "lambda.1se")
# trophic_cv_means <- cv_means(as.matrix(block[[2]][,trophic_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_trophic_min, trophic_cv_means)
t.test(error_trophic_1se, trophic_cv_means)
t.test(error_trophic_min, error_trophic_min_both)
t.test(error_trophic_min, error_trophic_min_boots)
t.test(error_trophic_min_both, error_trophic_min_boots)
# No difference between boots and both, GIS has the worst fit

# swimming
# bootstrap_swimming <- bootstrap_glmnet(as.matrix(block[[1]]),as.matrix(block[[2]][,swimming_traits$Marius.trait.letter.sc]), weights = block[[5]], core = 6, bootstraps = 500)
sum_swimming_min <- summarize_bootstraps(bootstrap_swimming, lambda = "lambda.min")
sum_swimming_1se <- summarize_bootstraps(bootstrap_swimming, lambda = "lambda.1se")
sum_range_swimming <- summarize_range_bootstep(bootstrap_swimming)
table_swimming_min <- summary_tables(sum_swimming_min)
table_swimming_1se <- summary_tables(sum_swimming_1se)
error_swimming_min <- cv_error(bootstrap_swimming, lambda = "lambda.min")
error_swimming_1se <- cv_error(bootstrap_swimming, lambda = "lambda.1se")
# swimming_cv_means <- cv_means(as.matrix(block[[2]][,swimming_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_swimming_min, swimming_cv_means)
t.test(error_swimming_1se, swimming_cv_means)
t.test(error_swimming_min, error_swimming_min_both)
t.test(error_swimming_min, error_swimming_min_boots)
t.test(error_swimming_min_both, error_swimming_min_boots)
# Both has the least error, GIS and boots are equal

# armor
# bootstrap_armor <- bootstrap_glmnet(as.matrix(block[[1]]),as.matrix(block[[2]][,armor_traits$Marius.trait.letter.sc]), weights = block[[5]], core = 6, bootstraps = 500)
sum_armor_min <- summarize_bootstraps(bootstrap_armor, lambda = "lambda.min")
sum_armor_1se <- summarize_bootstraps(bootstrap_armor, lambda = "lambda.1se")
sum_range_armor <- summarize_range_bootstep(bootstrap_armor)
table_armor_min <- summary_tables(sum_armor_min)
table_armor_1se <- summary_tables(sum_armor_1se)
error_armor_min <- cv_error(bootstrap_armor, lambda = "lambda.min")
error_armor_1se <- cv_error(bootstrap_armor, lambda = "lambda.1se")
# armor_cv_means <- cv_means(as.matrix(block[[2]][,armor_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_armor_min, armor_cv_means)
t.test(error_armor_1se, armor_cv_means)
t.test(error_armor_min, error_armor_min_both)
t.test(error_armor_min, error_armor_min_boots)
t.test(error_armor_min_both, error_armor_min_boots)
# GIS is the lowest, boots and both are the same and worse than means

# bio relevance
# bootstrap_bio_relevance <- bootstrap_glmnet(as.matrix(block[[1]]),as.matrix(block[[2]][,bio_relevance_traits$Marius.trait.letter.sc]), weights = block[[5]], core = 6, bootstraps = 500)
sum_bio_relevance_min <- summarize_bootstraps(bootstrap_bio_relevance, lambda = "lambda.min")
sum_bio_relevance_1se <- summarize_bootstraps(bootstrap_bio_relevance, lambda = "lambda.1se")
sum_bio_relevance_range <- summarize_range_bootstep(bootstrap_bio_relevance)
table_bio_relevance_min <- summary_tables(sum_bio_relevance_min)
table_bio_relevance_1se <- summary_tables(sum_bio_relevance_1se)
error_bio_relevance_min <- cv_error(bootstrap_bio_relevance, lambda = "lambda.min")
error_bio_relevance_1se <- cv_error(bootstrap_bio_relevance, lambda = "lambda.1se")
# bio_relevance_cv_means <- cv_means(as.matrix(block[[2]][,bio_relevance_traits$Marius.trait.letter.sc]), bootstraps = 500) 
t.test(error_bio_relevance_min, bio_relevance_cv_means)
t.test(error_bio_relevance_1se, bio_relevance_cv_means)
t.test(error_bio_relevance_min, error_bio_relevance_min_both)
t.test(error_bio_relevance_min, error_bio_relevance_min_boots)
t.test(error_bio_relevance_min_both, error_bio_relevance_min_boots)

# flextable(table_armor_min[[2]]) %>% save_as_docx( path = "C:/Users/trist/OneDrive/Desktop/test_flextable.docx")
# make pretty 
# formatC(table_armor_min[[2]]$max, format = "e", digits = 2)
# formatC(sum_range_armor, format = "e", digits = 2)

vote_count_min <- vote_count_table(list(table_min, table_trophic_min, table_swimming_min, table_armor_min, table_bio_relevance_min),
                                   name_order = c("All traits", "Trophic", "Swimming", "Armor", "Benthic - Limnetic"))


#### make tables ####
# pretty_tables(vote_count_min, concatenate = T)
# pretty_tables(table_swimming_min)
# pretty_tables(table_armor_min)
# pretty_tables(table_trophic_min)
# pretty_tables(table_min)
# pretty_tables(table_min, concatenate = T, keep_cols = c("mean"))
# pretty_tables(table_swimming_min, concatenate = T, keep_cols = c("mean"))
# pretty_tables(table_trophic_min, concatenate = T, keep_cols = c("mean"))
# pretty_tables(table_armor_min, concatenate = T, keep_cols = c("mean"))

# GGplots of CV error


# n <- 200  # number of data points
# 
# data <- data.frame(
#   Trait = rep(c("All Traits", "Swimming", "Armor", "Trophic"), each = n/4),
#   Method = sample(rep(c("Manual", "Automated"), each = n/2), n),
#   Metric = c(rnorm(n/2, mean = 50, sd = 10), rnorm(n/2, mean = 60, sd = 15))
# )

# data$x <- 1  # Adding a constant x value for all data points

# Updated plot using geom_dotplot
# p <- ggplot(data, aes(x = x, y = Metric)) +
#   geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 1.5, dotsize = 1.0, stackratio = 1, method = "dotdensity") +
#   stat_summary(fun = mean, geom = "hline", aes(yintercept = ..y..), 
#                color = "red", linetype = "dashed", size = 1.5) +
#   facet_grid(Trait ~ Method) +
#   theme_minimal() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   labs(title = "Metric Distribution by Trait and Method",
#        y = "Metric Value")

# p <- ggplot(data, aes(x = Method, y = Metric, fill = Method)) +
#   geom_boxplot() +
#   facet_grid(Trait ~ .) +  # Facets by Trait, with Method as x-axis within each facet
#   theme_minimal() +
#   labs(title = "Metric Distribution by Trait and Method",
#        x = "Method",
#        y = "Metric Value") +
#   scale_fill_brewer(palette = "Pastel1")  # Optional: nicer colors
# 
# # Display the plot
# print(p)

# add third factor
# 
# data$Source <- rep(c("Trait Mean", "GIS", "Both", "Local"), length.out = nrow(data))
# data <- data[order(data$Trait, data$Source), ] 
# 
# 
# p <- ggplot(data, aes(x = "", y = Metric, fill = Source)) +
#   geom_boxplot() +
#   facet_grid(Source ~ Trait) +  # Facets by Source and Trait
#   theme_minimal() +
#   labs(title = "Metric Distribution by Trait and Source",
#        x = "",  # No x-axis title needed
#        y = "Metric Value") +
#   scale_fill_brewer(palette = "Pastel1")  # Pleasant color scheme
# 
# # Display the plot
# print(p)

n <- 500
data1 <- data.frame(
  Trait = c(rep("All Traits",n*4), rep("Trophic",n*4), rep("Swimming",n*4), rep("Armor",n*4)),
  Metric = c(all_cv_means,error_min,error_min_both,error_min_boots,
             trophic_cv_means,error_trophic_min,error_trophic_min_both,error_trophic_min_boots,
             swimming_cv_means,error_swimming_min,error_swimming_min_both,error_swimming_min_boots,
             armor_cv_means,error_armor_min,error_armor_min_both,error_armor_min_boots), 
  Source = rep(c(rep("Trait Mean",n), rep("GIS",n), rep("Both",n), rep("Boots",n)), 4)
)

p <- ggplot(data1, aes(x = Source, y = Metric, fill = Source)) +
  geom_boxplot() +
  facet_wrap(~ Trait, scales = "free_y") +  # Nested facets
  theme_minimal() +
  labs(title = "Error Distribution by Trait and Data Source",
       x = "Data Source",
       y = "Error") +
  scale_fill_brewer(palette = "Pastel1")  # Use a pleasant color scheme

# Display the plot
print(p)


# Test standardized coefficients

X_sd <- apply(as.matrix(block[[1]]), 2, sd)

diag(X_sd) %*% bootstrap_armor[[2]]$glmnet.fit$beta[[1]][,10]



#### Figures ####

t <- cbind(block_both[[1]], block_both[[2]], block_both[[6]])
pivot_longer(, cols = )

ggplot(data = block[[2]]) +
  geom_point()
  geom_smooth(aes(), method = "lm") +
  facet_grid(rows = )

# 
# 
# #### old glmnets, not bootstrapped ####
# 
# # x <- data.frame(runif(100),runif(100),runif(100),runif(100),runif(100))
# # y <- data.frame(runif(100),runif(100),runif(100))
# # 
# # cvfit <- cv.glmnet(x = as.matrix(x), y = as.matrix(y), keep = T, family = "mgaussian")
# # mean(unlist(abs(cvfit$fit.preval[,,cvfit$index["min",]] - y)^2))
# # mean(unlist(abs(cvfit$fit.preval[,,41] - y)^2))
# # ## [1] 0.08803571
# # cvfit
# # plot(cvfit)
# #  # keep prevalidated array
# #    cvf1 <- cv.glmnet(x = as.matrix(mtcars[, c("disp", "hp", "mpg")]), 
# #                       y = mtcars$am, family = "binomial", keep = TRUE)
# #  dim(mtcars)
# # # [1] 32 11
# #  length(cvf1$lambda)
# # # [1] 84
# #  # leave-n out fitted predictions
# #    # 84 columns, 2 columns padded with NAs
# #    dim(cvf1$fit.preval)
# # # [1] 32 86
# #  # performance of cross-validated model predictions
# #    round(mtcars$am - cvf1$fit.preval[, cvf1$lambda == cvf1$lambda.min])
# # #  [1]  1  1  0  0  0  0  0  0 -1  0  0  0  0  0  0
# # # [16]  0  0  0  0  0 -1  0  0  0  0  0  0  0  1  0
# # # [31]  0  0
# #  cvf1$foldid
# 
# # I should replace the following code with the bootstrap and summarizing functions 
# 
# # model suffix is _ub (unscaled blocks)
# model_ub <- cv.glmnet(as.matrix(block[[1]]), as.matrix(block[[2]]), family = "mgaussian", nfolds = 5, gamma = c(0.5), weights = block[[5]], keep = TRUE)
# coef(model_ub, s = "lambda.min")
# coef(model_ub, s = "lambda.1se")
# model_ub
# plot(model_ub)
# # I think these are the fitted values
# model_ub$fit.preval
# # this is the index of the best lambda
# model_ub$index["min",]
# # this is the best fitted values
# model_ub$fit.preval[,,model_ub$index["min",]]
# # if relax = true you need the [[1]], i.e. model_ub$fit.preval[[1]][,,model_ub$index["min",]]
# # compare them
# model_ub$fit.preval[10:23,10:23,model_ub$index["min",]]
# block[[2]][10:23,10:23]
# model_ub$fit.preval[1:3,1:3,73]
# mean(unlist((model_ub$fit.preval[,,model_ub$index["min",]] - block[[2]])^2)) # this doesn't match the measure
# # 
# # # Multinomial
# # n = 500
# # p = 30
# # nzc = trunc(p/10)
# # x = matrix(rnorm(n * p), n, p)
# # beta3 = matrix(rnorm(30), 10, 3)
# # beta3 = rbind(beta3, matrix(0, p - 10, 3))
# # f3 = x %*% beta3
# # p3 = exp(f3)
# # p3 = p3/apply(p3, 1, sum)
# # g3 = glmnet:::rmult(p3)
# # set.seed(10101)
# # cvfit = cv.glmnet(x, g3, family = "multinomial", keep = TRUE)
# # plot(cvfit)
# # title("Multinomial Family", line = 2.5)
# # cvfit$fit.preval
# # # my stackoverflow example
# # x <- data.frame(runif(100),runif(100),runif(100),runif(100),runif(100))
# # y <- data.frame(runif(100),runif(100),runif(100))
# # 
# # cvfit <- cv.glmnet(x = as.matrix(x), y = as.matrix(y), keep = T, family = "mgaussian")
# # mean(unlist((cvfit$fit.preval[,,cvfit$index["min",]] - y)^2))
# 
# 
# # do many values of lambda near the optimal value from cv.glmnet to see if accuracy can be improved
# model_ub_fine <- cv.glmnet(as.matrix(block[[1]]), as.matrix(block[[2]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5),
#                                    lambda = seq(0.1,model_ub$relaxed$lambda.1se+0.1,0.01), weights = block[[5]])
# model_ub_fine
# coef(model_ub_fine, s = "lambda.min")
# coef(model_ub_fine, s = "lambda.1se")
# plot(model_ub_fine)
# 
# trophic_ub <- cv.glmnet(as.matrix(block[[1]]),
#                          as.matrix(block[[2]][,trophic_traits$Marius.trait.letter.sc]), # yoel trait names: 
#                          nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])
# trophic_ub
# coef(trophic_ub, s = "lambda.min")
# coef(trophic_ub, s = "lambda.1se")
# plot(trophic_ub)
# 
# 
# swimming_ub <- cv.glmnet(as.matrix(block[[1]]),
#                           as.matrix(block[[2]][,swimming_traits$Marius.trait.letter.sc]),
#                           nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])
# 
# swimming_ub
# coef(swimming_ub, s = "lambda.min")
# coef(swimming_ub, s = "lambda.1se")
# plot(swimming_ub)
# 
# 
# armor_ub <- cv.glmnet(as.matrix(block[[1]]),
#                        as.matrix(block[[2]][,armor_traits$Marius.trait.letter.sc]),
#                        nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights =  block[[5]])
# armor_ub
# coef(armor_ub, s = "lambda.min")
# coef(armor_ub, s = "lambda.1se")
# plot(armor_ub)
# 
# 
# ### mean 0, sd 1 allometric scaled blocks ###
# 
# # model suffix is _sd1
# 
# model_sd1 <- cv.glmnet(as.matrix(block[[3]]), as.matrix(block[[4]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])
# coef(model_sd1, s = model_sd1$relaxed$lambda.min)
# model_sd1
# plot(model_sd1)
# 
# # do many values of lambda near the optimal value from cv.glmnet to see if accuracy can be improved
# model_sd1_scaled_fine <- cv.glmnet(as.matrix(block[[3]]), as.matrix(block[[4]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5),
#                                lambda = seq(model_sd1$relaxed$lambda.min-.5,model_sd1$relaxed$lambda.1se,0.01), weights = block[[5]])
# model_sd1_scaled_fine
# coef(model_sd1_scaled_fine, s = "lambda.min")
# coef(model_sd1_scaled_fine, s = "lambda.1se")
# plot(model_sd1_scaled_fine)
# 
# trophic_sd1 <- cv.glmnet(as.matrix(block[[3]]),
#                            as.matrix(block[[4]][,trophic_traits$Marius.trait.letter.sc]), # yoel trait names: 
#                            nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])
# trophic_sd1
# coef(trophic_sd1, s = "lambda.min")
# coef(trophic_sd1, s = "lambda.1se")
# plot(trophic_sd1)
# 
# 
# swimming_sd1 <- cv.glmnet(as.matrix(block[[3]]),
#                             as.matrix(block[[4]][,swimming_traits$Marius.trait.letter.sc]),
#                             nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])
# 
# swimming_sd1
# coef(swimming_sd1, s = "lambda.min")
# coef(swimming_sd1, s = "lambda.1se")
# plot(swimming_sd1)
# 
# 
# armor_sd1 <- cv.glmnet(as.matrix(block[[3]]),
#                          as.matrix(block[[4]][,armor_traits$Marius.trait.letter.sc]),
#                          nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block[[5]])
# armor_sd1
# coef(armor_sd1, s = "lambda.min")
# coef(armor_sd1, s = "lambda.1se")
# plot(armor_sd1)
# 
# 
# ### no allometric scaling blocks ###
# 
# # model suffix is _m 
# model_m <- cv.glmnet(as.matrix(block_m[[1]]), as.matrix(block_m[[2]]), family = "mgaussian", nfolds = 5, relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
# coef(model_m, s = "lambda.min")
# coef(model_m, s = "lambda.1se")
# model_m
# plot(model_m)
# 
# # do many values of lambda near the optimal value from cv.glmnet to see if accuracy can be improved
# model_m_fine <- cv.glmnet(as.matrix(block_m[[1]]), as.matrix(block_m[[2]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5),
#                            lambda = seq(0.1,model_m$relaxed$lambda.1se+0.1,0.01), weights = block_m[[5]])
# model_m_fine
# coef(model_m_fine, s = "lambda.min")
# coef(model_m_fine, s = "lambda.1se")
# plot(model_m_fine)
# 
# trophic_m <- cv.glmnet(as.matrix(block_m[[1]]),
#                         as.matrix(block_m[[2]][,trophic_traits$Marius.trait.letter]), # yoel trait names: 
#                         nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
# trophic_m
# coef(trophic_m, s = "lambda.min")
# coef(trophic_m, s = "lambda.1se")
# plot(trophic_m)
# 
# 
# swimming_m <- cv.glmnet(as.matrix(block_m[[1]]),
#                          as.matrix(block_m[[2]][,swimming_traits$Marius.trait.letter]),
#                          nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
# 
# swimming_m
# coef(swimming_m, s = "lambda.min")
# coef(swimming_m, s = "lambda.1se")
# plot(swimming_m)
# 
# 
# armor_m <- cv.glmnet(as.matrix(block_m[[1]]),
#                       as.matrix(block_m[[2]][,armor_traits$Marius.trait.letter]),
#                       nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
# armor_m
# coef(armor_m, s = "lambda.min")
# coef(armor_m, s = "lambda.1se")
# plot(armor_m)
# 
# 
# #### no allometric scaling mean 0, sd 1 blocks ###
# 
# # model suffix is _m.sd 
# model_m.sd <- cv.glmnet(as.matrix(block_m[[3]]), as.matrix(block_m[[4]]), family = "mgaussian", nfolds = 5, relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
# coef(model_m.sd, s = "lambda.min")
# coef(model_m.sd, s = "lambda.1se")
# model_m.sd
# plot(model_m.sd)
# 
# # do many values of lambda near the optimal value from cv.glmnet to see if accuracy can be improved
# model_m.sd_fine <- cv.glmnet(as.matrix(block_m[[3]]), as.matrix(block_m[[4]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5),
#                           lambda = seq(0.1,model_m.sd$relaxed$lambda.1se+0.1,0.01), weights = block_m[[5]])
# model_m.sd_fine
# coef(model_m.sd_fine, s = "lambda.min")
# coef(model_m.sd_fine, s = "lambda.1se")
# plot(model_m.sd_fine)
# 
# trophic_m.sd <- cv.glmnet(as.matrix(block_m[[3]]),
#                        as.matrix(block_m[[4]][,trophic_traits$Marius.trait.letter]), # yoel trait names: 
#                        nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
# trophic_m.sd
# coef(trophic_m.sd, s = "lambda.min")
# coef(trophic_m.sd, s = "lambda.1se")
# plot(trophic_m.sd)
# 
# 
# swimming_m.sd <- cv.glmnet(as.matrix(block_m[[3]]),
#                         as.matrix(block_m[[4]][,swimming_traits$Marius.trait.letter]),
#                         nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
# 
# swimming_m.sd
# coef(swimming_m.sd, s = "lambda.min")
# coef(swimming_m.sd, s = "lambda.1se")
# plot(swimming_m.sd)
# 
# 
# armor_m.sd <- cv.glmnet(as.matrix(block_m[[3]]),
#                      as.matrix(block_m[[4]][,armor_traits$Marius.trait.letter]),
#                      nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5), weights = block_m[[5]])
# armor_m.sd
# coef(armor_m.sd, s = "lambda.min")
# coef(armor_m.sd, s = "lambda.1se")
# plot(armor_m.sd)
# 
# #### individual values, SL as a predictor, unscaled ###
# 
# # for individual values, use foldid in cv.glmnet, make folds based on lake of origin
# # Use match to get the integer values
# integer_vector <- match(block_i[[6]], unique(block_i[[6]]))
# 
# # model suffix is _i. # Having issues with weights being used for some reason.
# # make predictor and response by moving SL.mm. Ithink NAs are the issue
# 
# i_sl <- data.frame(block_i[[4]]) %>% 
#   select(SL.mm)
# i_response <- data.frame(block_i[[4]]) %>% 
#   select(-SL.mm) %>% 
#   as.matrix()
# i_predictor <- as.matrix(data.frame(block_i[[3]], i_sl))
# i_predictor <- as.matrix(data.frame(block_i[[3]]))
# 
# 
# # drop NAs
# i_predictor <- i_predictor[rowSums(is.na(i_response)) <= 0,]
# i_response <- i_response[rowSums(is.na(i_response)) <= 0,]
# 
# model_i <- cv.glmnet(i_predictor, i_response, family = "mgaussian", relax = TRUE, gamma = c(0.5),
#                      nfolds = max(nfoldid = ceiling(integer_vector/11)) , nfoldid = ceiling(integer_vector/11)) # I can use the ceiling() to change how many lakes per fold. Results dont really change much though
# #, family = "mgaussian", nfolds = 5, relax = TRUE, gamma = c(0.5))
# coef(model_i, s = "lambda.min")
# coef(model_i, s = "lambda.1se")
# model_i
# plot(model_i)
# 
# 
# model_i <- cv.glmnet(as.matrix(block_i[[1]]), as.matrix(block_i[[2]]), family = "mgaussian", relax = TRUE, gamma = c(0.5),
#                      nfolds = max(nfoldid = ceiling(integer_vector/11)) , nfoldid = ceiling(integer_vector/11)) 
# 
# # # do many values of lambda near the optimal value from cv.glmnet to see if accuracy can be improved
# # model_i_fine <- cv.glmnet(as.matrix(block_i[[3]]), as.matrix(block_i[[4]]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5),
# #                              lambda = seq(0.1,model_i$relaxed$lambda.1se+0.1,0.01), weights = block_i[[5]])
# # model_i_fine
# # coef(model_i_fine, s = "lambda.min")
# # coef(model_i_fine, s = "lambda.1se")
# # plot(model_i_fine)
# 
# trophic_i <- cv.glmnet(i_predictor, i_response[,trophic_traits$Marius.trait.letter], # yoel trait names: 
#                           nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5))
# trophic_i
# coef(trophic_i, s = "lambda.min")
# coef(trophic_i, s = "lambda.1se")
# plot(trophic_i)
# 
# 
# swimming_i <- cv.glmnet(i_predictor, i_response[,swimming_traits$Marius.trait.letter],
#                            nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5))
# 
# swimming_i
# coef(swimming_i, s = "lambda.min")
# coef(swimming_i, s = "lambda.1se")
# plot(swimming_i)
# 
# 
# armor_i <- cv.glmnet(i_predictor, i_response[,armor_traits$Marius.trait.letter],
#                         nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5))
# armor_i
# coef(armor_i, s = "lambda.min")
# coef(armor_i, s = "lambda.1se")
# plot(armor_i)
# 
# #### individual values, SL as a predictor, mean 0, sd 1 ###
# 
# #### interactions, best mean value method ###
# # might be interesting to try glmnet with interaction terms: library(glinternet), tutorial: https://strakaps.github.io/post/glinternet/#fnref1
# 
# library(glinternet)
# 
# 
# t <- glinternet.cv(i_predictor, i_response[,1], numLevels = rep(1, ncol(i_predictor)))
# t2 <- glinternet.cv(as.matrix(block[[1]]), as.vector(block[[2]][,1]), numLevels = rep(1, ncol(block[[1]])))
# plot(t)
# plot(t2)
# t2
# #### interactions, best individual values method ###
# 
# 
# 
# 
# 
# 
# # # scaling just morpho trait is still not great, not converging on good predictions
# # model_scaled_morpho <- cv.glmnet(as.matrix(block_enviro), as.matrix(block_morpho_scaled), family = "mgaussian", nfolds = dim(block)[1], relax = TRUE, gamma = c(0.5))
# # coef(model_scaled_morpho, s = model_scaled_morpho$lambda.min)
# 
# # both scaled
# 
# 
# 
# # not ready to do parasite yet
# # parasite_sd1 <- cv.glmnet(as.matrix(block_enviro_scaled), as.matrix(block_morpho_scaled[,c("Thersitina.gill.copepod.0.1","Unionidae.gill.clam.0.1",
# #                                                                                              "Schistocephalus.cestode.bodycavity.COUNT","Eustrongylides.nematode.bodycavity.COUNT")]), nfolds = 5, family = "mgaussian", relax = TRUE, gamma = c(0.5))
# # parasite_sd1
# # coef(parasite_sd1, s = parasite_sd1$relaxed$lambda.min)
# # coef(parasite_sd1, s = parasite_sd1$relaxed$lambda.1se)
# # plot(parasite_sd1)
# 
# 
# #### Debugging cv.glmnet difference between lambda.min and what is shown when the model is called directly
# # pretty sure its because multiple values are calculated in each fold and when called the single model trained on all the data is used
# a <- runif(100)
# b <- runif(100)
# c <- runif(100)
# d <- b + runif(100)/10
# e <- a + runif(100)/10
# 
# test <- cv.glmnet(cbind(a,b,c), cbind(d,e), family = "mgaussian", relax = TRUE, gamma = 0.5)
# test
# test$gamma
