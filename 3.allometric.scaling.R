# === general ==================================================================

## 3.allometric.scaling.R
## perform allometric scaling of traits
#
# authors: T. Kosciuch adapting code provided by Krista Oke
# 
# R version
R.version.string
# "R version 4.1.2 (2021-11-01)"
#
#

# === table of contents ==================================================================

# 1. load libraries and import data
#    - Use the scale measure to convert pixels into cm
#
# 2. Create NAs for tristan missing lateral and dorsal spines, select yoel traits to scale
#
# 3. Size correction

# ===  1. load libraries and import data ==================================================================
library(stringr)
library(lme4)
library(dplyr)
library(readxl)

lateral <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/lateral_lengths.csv") %>% rename(SL.mm = SL)
ventral <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/ventral_lengths.csv")
yoel_lengths <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Yoels 2017 project/OG files/01.univMorph.csv") %>% filter(habitat == "lake")
lake_stream_yoel_lengths <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Yoels 2017 project/OG files/01.univMorph.csv")
merge_yoel_tristan <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/link_yoel_and_marius_traits.csv") %>% filter(In.both.datasets == "y")
merge_yoel_FWA <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Enviro data/match_yoels_samples_to_FWA.csv")
fish_to_remove <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/lateral_specimens_that_shouldnt_exist.csv")
kosciuch_gill_raker_n <- read_excel("C:/Users/trist/OneDrive - Loyola University Chicago/Tristan_LandscapeSampling/Gill Rakers/Gill Raker Counts + Measures & Gravid Scoring.xlsx", sheet = "gill.raker.counts") %>% 
  select(id, right.side.gillraker.number) %>% 
  mutate(right.side.gillraker.number = as.numeric(right.side.gillraker.number)) # NAs are ok
stuart_gill_raker_n <- read.csv("C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/Yoels 2017 project/OG files/01.univMorph.csv") %>% filter(habitat == "lake") %>% 
  select(fishID.univ, Right.side.gill.raker.number.insitu) %>% 
  rename(id = fishID.univ,
         right.side.gillraker.number = Right.side.gill.raker.number.insitu)
gill_raker_n <- bind_rows(kosciuch_gill_raker_n, stuart_gill_raker_n)

# save armor.plate.count as a vector
# armor.plate.count <- lateral$armor.plate.count
# convert lengths to mm from pixels
lateral[c(LETTERS[1:19],"SL.mm")] <- (lateral[c(LETTERS[1:19],"SL.mm")] / lateral$scale.15.mm) * 15
ventral[c(letters[1:13])] <- (ventral[letters[1:13]] / ventral$scale.15.mm) * 15
# add armor plate counts back in
# lateral$armor.plate.count <- armor.plate.count
# add a watershed column. Watershed here is actually the lake of origin, not the true watershed
lateral$watershed <- substr(lateral$id,1,nchar(lateral$id)-6)
ventral$watershed <- substr(ventral$id,1,nchar(ventral$id)-6)
gill_raker_n$watershed <- substr(gill_raker_n$id,1,nchar(gill_raker_n$id)-2)

# select traits from yoel dataset
traits_to_get <- merge_yoel_tristan %>% filter(In.both.datasets == "y")
yoel_lengths <- yoel_lengths %>% select(c(fishID.univ,watershed,standard.len.mm.lateral,Left.Side.Plate.Count,traits_to_get$Yoels.trait))
lake_stream_yoel_lengths <- lake_stream_yoel_lengths %>% select(c(fishID.univ,watershed,habitat,standard.len.mm.lateral,Left.Side.Plate.Count,traits_to_get$Yoels.trait)) %>% 
  mutate(ls.watershed = paste0(habitat, watershed))


# ===  2. Create NAs for tristan missing spines, select yoel traits to scale ==================================================================
# NAs are needed for tristans dataset because when a spine was missing, I still placed a landmark
# because ObjectJ needed one to move to the next landmark. I placed the landmark right next to the 
# base of the spine so the calculated spine length will be very small.

#### lateral

# find missing first dorsal spines.View these dataframes and sort by trait length
dsp1 <- lateral[lateral$K < 2,] %>% select(id, K)
dsp1 <- dsp1[order(dsp1$K),]
# all individual up to and including LS00090SALM31 have a missing first dorsal spine
# also these have broken first dorsal spines: SC00792SALM25
# I checked between LS00090SALM31 and ST00097SALM17 and only SC00792SALM25 had a very obviously broken spine
dsp1_NAs <- dsp1$id[1:match("LS00090SALM31.JPG",dsp1$id)]
dsp1_NAs <- c(dsp1_NAs,"SC00792SALM25.JPG")

# second dorsal spine
dsp2 <- lateral[lateral$L < 2.5,] %>% select(id, L)
dsp2 <- dsp2[order(dsp2$L),]
# up to GR00039CAMB03 have broken second dorsal spines except for MT00294VICT18 and BB00096SALM24 which do not have clearly broken spines.
# CB00058CAMB11 also has a broken second dorsal spine
# up to ST00097SALM11 was checked
dsp2_NAs <- dsp2$id[1:match("GR00039CAMB03.JPG",dsp2$id)]
dsp2_NAs <- dsp2_NAs[!dsp2_NAs %in% c("MT00294VICT18.JPG","BB00096SALM24.JPG")]
dsp2_NAs <- c(dsp2_NAs,"CB00058CAMB11.JPG")

# change values to NA
lateral$K[lateral$id %in% dsp1_NAs] <- NA
lateral$L[lateral$id %in% dsp2_NAs] <- NA

#### ventral

# ventral

# find missing first dorsal spines.View these dataframes and sort by trait length
vsp1 <- ventral[ventral$l < 2,] %>% select(id, l) # one with NA. this NA is there because the spine is missing
vsp2 <- ventral[ventral$m < 2,] %>% select(id, m) # no broken or missing spines


# ===  3. size correction ==================================================================

# Size correction code provided by Krista Oke

# OKE ET AL COMMON GARDEN MANUSCRIPT METHOD
#"All relative warps and univariate shape traits were allometrically standardized to a common body size (Reist 1985; Lleonart et al. 2000) based on Ms = M0 (Ls / L0)^b, where Ms is the standardized trait value (mm), M0 is the non-standardized trait value (mm), Ls is the overall mean centroid size (for RWs) or standard length (mm, for univariate shape traits), and L0 is the centroid size or standard length (mm) of the individual. The common within-group slope, b, was calculated from a linear mixed model of log10(M0) regressed on log10(L0), with group included as a random factor (Reist 1985; Lleonart et al. 2000)." - Oke et al.    The size correction formula described above: Ms = M0(Ls/L0)^b

#variables:
#standard.length = data.1.ps.3$stl #the trait to be used for size correction
#watershed = substr(data.1.ps.3$univ.specimenID, 1, 4) #the random factor in the linear mixed model. This is the section
#dataframe.to.use = data.1.ps.3 #the data
#vector.of.columns = 2 #the columns in dataframe.to.use that are to be size corrected
#to.length = provide a body length to scale traits to 

#The output will have all the original columns from dataframe.to.use with the size-corrected columns added on at the end.
f.size.correct <- function(standard.length, vector.of.columns, dataframe.to.use, watershed, to.length = NA) {
  if(is.na(to.length)){
    #Calculate overall mean standard length (Ls)
    Ls <- mean(standard.length, na.rm = TRUE) 
  } else {
    Ls <- to.length
  }
  #Call individual standard length
  L0 <- standard.length
  #Calculate common, within-group slope, b, for each trait. 
  #Here, I am treating each section as a group for random factor
  b.vector.lmm <- vector()
  for (i in vector.of.columns) {
    abcd <- (dataframe.to.use[i])
    #b.model <- lmer((abcd[,])~(standard.length) + (1|watershed))
    #On 16 November, realized that this is wrong. Calculation of Beta needs log(M0) ~ Log(L0). This wasn't logged.
    #Re-running now with log10 values.
    b.model <- lmer(log10(abcd[,])~log10(standard.length) + (1|watershed))
    b <- coef(summary(b.model))[2,1]
    b.vector.lmm <- c(b.vector.lmm, b)
  }
  # size correct
  xx <- dataframe.to.use  
  columnnames <- colnames(xx)
  j=1
  for (i in vector.of.columns) {
    M0 <- xx[,i] #grab the appropriate column of data
    Ms = M0 * ((Ls/L0)^b.vector.lmm[j]) #size correction formula
    j=j+1
    columnnames <- c(columnnames, paste(colnames(xx[i]), "sc", sep = "."))
    xx <- cbind(xx, Ms)
  }
  colnames(xx) <- columnnames # Rename the columns in the temporary dataframe xx
  return(xx) #Output a new dataframe with the name provided in "outputfilename"
}
# end 4 f.size.correct

to.length <- mean(c(lateral$SL.mm,yoel_lengths$standard.len.mm.lateral))

## Merge yoel and tk datasets
# create yoel_data which has morpho and parasites

# use merge_yoel_tristan and filter for in both datasets and lateral
merge_yoel_tristan_lateral <- merge_yoel_tristan %>% filter(lateral.or.ventral == "lateral") %>% 
  filter(Yoels.trait %in% colnames(yoel_lengths))

yoel_lateral <- yoel_lengths %>% select("fishID.univ","watershed",merge_yoel_tristan_lateral$Yoels.trait)
matching_indices_lateral <- match(names(yoel_lateral), merge_yoel_tristan_lateral$Yoels.trait)
names(yoel_lateral) <- c("id", "watershed", merge_yoel_tristan_lateral$Marius.trait.letter[matching_indices_lateral[-c(1:2)]])

merged_lateral <- lateral %>% select("id", "watershed", merge_yoel_tristan_lateral$Marius.trait.letter) %>% 
  bind_rows(yoel_lateral)

lateral_traits_to_sc <- merge_yoel_tristan_lateral$Marius.trait.letter[!merge_yoel_tristan_lateral$Marius.trait.letter %in% c("id", "watershed", "SL.mm")]
merged_lateral_sc <- f.size.correct(merged_lateral$SL.mm, vector.of.columns = lateral_traits_to_sc, dataframe.to.use = merged_lateral, watershed = merged_lateral$watershed) %>% 
  select("id", "watershed", "SL.mm", paste0(lateral_traits_to_sc, ".sc"))

# ventral
merge_yoel_tristan_ventral <- merge_yoel_tristan %>% filter(lateral.or.ventral == "ventral") %>% 
  filter(Yoels.trait %in% colnames(yoel_lengths))

yoel_ventral <- yoel_lengths %>% select("fishID.univ","watershed",merge_yoel_tristan_ventral$Yoels.trait)
matching_indices_ventral <- match(names(yoel_ventral), merge_yoel_tristan_ventral$Yoels.trait)
names(yoel_ventral) <- c("id", "watershed", merge_yoel_tristan_ventral$Marius.trait.letter[matching_indices_ventral[-c(1:2)]])
SL.mm <- merged_lateral %>% select(id, SL.mm)
merged_ventral <- ventral %>% select("id", "watershed", merge_yoel_tristan_ventral$Marius.trait.letter) %>% 
  bind_rows(yoel_ventral) %>% 
  inner_join(SL.mm)
ventral_traits_to_sc <- merge_yoel_tristan_ventral$Marius.trait.letter[!merge_yoel_tristan_ventral$Marius.trait.letter %in% c("id", "watershed", "SL.mm")]
merged_ventral_sc <- f.size.correct(merged_ventral$SL.mm, vector.of.columns = ventral_traits_to_sc, dataframe.to.use = merged_ventral, watershed = merged_ventral$watershed) %>% 
  select("id", "watershed", "SL.mm", paste0(ventral_traits_to_sc, ".sc"))

# # before writing filter out the fish from "lateral_specimens_that_shouldnt_exist.csv" from "C:\Users\trist\OneDrive - Loyola University Chicago\Thesis\R\Morphological data"
merged_lateral_sc <- merged_lateral_sc %>% filter(!id %in% fish_to_remove$id)
merged_ventral_sc <- merged_ventral_sc %>% filter(!id %in% fish_to_remove$id)
merged_lateral <- merged_lateral %>% filter(!id %in% fish_to_remove$id)
merged_ventral <- merged_ventral %>% filter(!id %in% fish_to_remove$id)
# these specimens were already removed apparently

apply(X = is.na(merged_lateral_sc), MARGIN = 2, FUN = sum)
apply(X = is.na(merged_ventral_sc), MARGIN = 2, FUN = sum) # 95 l.sc missing spines but thats what the data is, maybe some were spine clips. They come from yoels data
apply(X = is.na(merged_lateral), MARGIN = 2, FUN = sum)
apply(X = is.na(merged_ventral), MARGIN = 2, FUN = sum) # 95 l.sc missing spines but thats what the data is, maybe some were spine clips. They come from yoels data

# 
# write.csv(merged_lateral_sc,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_lateral_sc.csv", row.names = F)
# write.csv(merged_ventral_sc,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_ventral_sc.csv", row.names = F)
# write.csv(merged_lateral,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_lateral.csv", row.names = F)
# write.csv(merged_ventral,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_ventral.csv", row.names = F)









# old stuff 

# yoel_lateral <- yoel_lengths %>% select(c("fishID.univ","watershed",merge_yoel_tristan_lateral$Yoels.trait)) %>%
#   rename_with(~(c("id","watershed",merge_yoel_tristan_lateral$Marius.trait.letter)), everything())
# yoel_parasite <- yoel_parasite %>% select(fishID.univ, Thersitina.gill.copepod.0.1,
#                                           Unionidae.gill.clam.0.1, 
#                                           Schistocephalus.cestode.bodycavity.COUNT,
#                                           Eustrongylides.nematode.bodycavity.COUNT) 
# yoel_data <- left_join(yoel_morpho, yoel_parasite,by = "fishID.univ") %>% 
#   mutate(id = fishID.univ,
#          dataset = "stuart.nee.2017")

# create tk_data which has morpho and parasites
# 
# tk_lengths <- inner_join(sc_lateral,sc_ventral, by = c("id","SL.mm","watershed")) %>% 
#   mutate(id = str_remove(id,".JPG"),
#          SL.mm = round(SL.mm, 3))
# 
# tk_data <- inner_join(tk_lengths, tk_parasite_body_cavity, by = "id") %>%
#   inner_join(tk_parasite_gills, by = "id") %>% 
#   inner_join(tk_raker_count, by = c("id","watershed")) %>% 
#   mutate(dataset = "kosciuch.masters")

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

# # create a dataset of tk_data with only the same columns as yoel_data
# merged <- tk_data %>% inner_join(lake_joining_info, by = "watershed") %>% 
#   select(c(colnames(yoel_data)[colnames(yoel_data) %in% colnames(tk_data)]), "WTRBDPLD") %>% 
#   bind_rows(yoel_data) %>% 
#   select(-fishID.univ)
# 
# # write.csv(merged, "C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/merged_kosciuch_stuart_morpho.csv", row.names = FALSE)
# 
# tk_ys_lateral <- lateral
# tk_ys_ventral <- ventral
# 
# lateral_size_correct <- f.size.correct(lateral$SL.mm, c(LETTERS[1:19],"armor.plate.count"), lateral, lateral$watershed, to.length = to.length)
# lateral_size_correct <- lateral_size_correct[,!(names(lateral_size_correct) %in% c("scale.15.mm",LETTERS[1:19],"armor.plate.count"))]
# colnames(lateral_size_correct)[2] <- "SL.mm"

# SL <- lateral_size_correct %>% select(id, SL.mm)
# ventral <- inner_join(ventral,SL)
# ventral_size_correct <- f.size.correct(ventral$SL.mm,letters[1:13],ventral,ventral$watershed, to.length = to.length)
# ventral_size_correct <- ventral_size_correct[,!(names(ventral_size_correct) %in% c("scale.15.mm",letters[1:13]))]
# 
# # size_correct_lake_stream_yoel_lengths <- f.size.correct(lake_stream_yoel_lengths$standard.len.mm.lateral, traits_to_get$Yoels.trait[!traits_to_get$Yoels.trait %in% c("standard.len.mm.lateral", "ls.watershed")], lake_stream_yoel_lengths, lake_stream_yoel_lengths$ls.watershed, to.length = to.length)
# 
# yoel_size_correct <- f.size.correct(yoel_lengths$standard.len.mm.lateral, traits_to_get$Yoels.trait[!traits_to_get$Yoels.trait %in% "standard.len.mm.lateral"], yoel_lengths, yoel_lengths$watershed, to.length = to.length) #
# ysc_to_drop <- traits_to_get$Yoels.trait[traits_to_get$Yoels.trait != "standard.len.mm.lateral"]
# yoel_size_correct <- yoel_size_correct[,!names(yoel_size_correct) %in% ysc_to_drop]
# 

# # Not sure why but cannot size correct gill raker number. There is hardly a correlation with SL though
# SL$id_trim <- substr(SL$id,1,nchar(SL$id)-4)
# gill_raker_n_sl <- inner_join(gill_raker_n,SL[,c("id_trim", "SL.mm")], by=c("id" = "id_trim"))
# gill_raker_n_sl$right.side.gillraker.number <- as.numeric(gill_raker_n_sl$right.side.gillraker.number)
# gill_raker_n_sl <- gill_raker_n_sl %>% filter(!is.na(right.side.gillraker.number))
# # gill_n_size_correct <- f.size.correct(gill_raker_n_sl$SL.mm, "right.side.gillraker.number", gill_raker_n_sl, gill_raker_n_sl$watershed, to.length = to.length)
# # show lack of correlation with SL, across the range of most fish there is a 0.4 difference in count according to slope of LM
# with(gill_raker_n_sl, plot(right.side.gillraker.number ~ SL.mm))
# with(gill_raker_n_sl, summary(lm(right.side.gillraker.number ~ SL.mm + watershed)))

# # ===  4. set NAs to the within group (lake) population mean ==================================================================
# 
# # columns that have NAs
# colnames(yoel_size_correct)[apply(yoel_size_correct, 2, anyNA)]
# colnames(lateral_size_correct)[apply(lateral_size_correct, 2, anyNA)]
# 
# 
# yoel_sc_no_NA <- yoel_size_correct %>%
#   group_by(watershed) %>%
#   mutate(dorsal.spine.1.len.mm.lateral.sc = case_when(is.na(dorsal.spine.1.len.mm.lateral.sc) ~ mean(dorsal.spine.1.len.mm.lateral.sc, na.rm=TRUE),
#                                                       TRUE ~ as.numeric(dorsal.spine.1.len.mm.lateral.sc)),
#          dorsal.spine.2.len.mm.lateral.sc = case_when(is.na(dorsal.spine.2.len.mm.lateral.sc) ~ mean(dorsal.spine.2.len.mm.lateral.sc, na.rm=TRUE),
#                                                       TRUE ~ as.numeric(dorsal.spine.2.len.mm.lateral.sc)),
#          pelvic.girdle.width.mm.sc = case_when(is.na(pelvic.girdle.width.mm.sc) ~ mean(pelvic.girdle.width.mm.sc, na.rm=TRUE),
#                                                TRUE ~ as.numeric(pelvic.girdle.width.mm.sc)),
#          body.width.anal.2.mm.sc = case_when(is.na(body.width.anal.2.mm.sc) ~ mean(body.width.anal.2.mm.sc, na.rm=TRUE),
#                                              TRUE ~ as.numeric(body.width.anal.2.mm.sc)),
#          Left.Side.Pelvic.Spine.Length.mm.sc = case_when(is.na(Left.Side.Pelvic.Spine.Length.mm.sc) ~ mean(Left.Side.Pelvic.Spine.Length.mm.sc, na.rm=TRUE),
#                                                          TRUE ~ as.numeric(Left.Side.Pelvic.Spine.Length.mm.sc)),
#          Right.Side.Pelvic.Spine.Length.mm.sc = case_when(is.na(Right.Side.Pelvic.Spine.Length.mm.sc) ~ mean(Right.Side.Pelvic.Spine.Length.mm.sc, na.rm=TRUE),
#                                                           TRUE ~ as.numeric(Right.Side.Pelvic.Spine.Length.mm.sc))
#          )
# 
# lateral_sc_no_NA <- lateral_size_correct %>%
#   group_by(watershed) %>%
#   mutate(K.sc = case_when(is.na(K.sc) ~ mean(K.sc, na.rm=TRUE),
#                           TRUE ~ as.numeric(K.sc)),
#          L.sc = case_when(is.na(L.sc) ~ mean(L.sc, na.rm=TRUE),
#                           TRUE ~ as.numeric(L.sc)))
# 
# # No more NAs
# colnames(yoel_sc_no_NA)[apply(yoel_sc_no_NA, 2, anyNA)]
# colnames(lateral_sc_no_NA)[apply(lateral_sc_no_NA, 2, anyNA)]
# 

# # ===  5. write data to csv ==================================================================
# 
# write.csv(lateral_sc_no_NA,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/size_corrected_lateral_lengths_mm.csv", row.names = FALSE)
# write.csv(ventral_size_correct,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/size_corrected_ventral_lengths_mm.csv", row.names = FALSE)
# write.csv(yoel_sc_no_NA,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/yoels_size_corrected_trait_lengths_mm.csv", row.names = FALSE)
# # write.csv(gill_n_size_correct,"C:/Users/trist/OneDrive - Loyola University Chicago/Thesis/R/Morphological data/gill_n_size_correct.csv", row.names = FALSE)

