# === general ==================================================================

## 1.join.FWA_lakes.to.site_data.R 
## Join FWA lake polygon information to the list of lakes I actually sampled
# authors: T. Kosciuch
# 
# R version
R.version.string
# "R version 4.1.2 (2021-11-01)"
#

# =====================================================================


lakes <- st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/2020 Sampling/2020 sample lakes/2020 sample lakes.shp")


already_sampled <- read.csv("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Summer 2020/Summer 2020 datasheets/Tristan already sampled lakes.csv") %>% 
  mutate(lake = toupper(ï..lake)) %>%
  group_by(lake,id) %>%
  summarise(n_tsb = sum(X.tsb.sampled))

# rename 'id' 'into WATERBODY_KEY_GROUP_CODE_50K' to allow joining with lakes
names(already_sampled)[names(already_sampled) == 'id'] <- 'WTRBDKGRPC'


for(i in 1:length(already_sampled$WTRBDKGRPC)){
  
  if(i ==1){joined <- data.frame()}
  
  # "00316SALM" has an SLP lake with the same id for two different polygons. I'm using the WATERBODY_POLY_ID's described in the Feb 4 log.
  if(already_sampled$WTRBDKGRPC[i] == "00316SALM"){
    # print(grouped$WTRBDKGRPC[i])
    
    if(already_sampled$lake[i] == "SLP"){
      
      temp_lakes <- lakes %>% filter(WTRBDPLD == 710017713)
      temp_grouped <- already_sampled[i,] %>% filter(lake == "SLP")
      
      joined <- rbind(joined,left_join(temp_grouped,temp_lakes,by ="WTRBDKGRPC"))
      
    } else {
      
      temp_lakes <- lakes %>% filter(WTRBDPLD == 710017394)
      temp_grouped <- already_sampled[i,] %>% filter(lake == "ORMOND")
      
      joined <- rbind(joined,left_join(temp_grouped,temp_lakes,by ="WTRBDKGRPC"))
    }
  } else {
    joined <- rbind(joined,left_join(already_sampled[i,],lakes,by ="WTRBDKGRPC"))
  }
  
} ; rm(temp_grouped,temp_lakes,i)

st_write(joined,"C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/2020 Sampling/2020 sample lakes/2020 sample lakes FWA join.shp",
         driver = "ESRI Shapefile",append = FALSE)

st_read("C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Thesis/GIS data/2020 Sampling/2020 sample lakes/2020 sample lakes FWA join.shp")