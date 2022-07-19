#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#        Colinearity of predictors          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load environmental space
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km_normalized.grd"))

# remove SoilTemp and Aridity to remove vif of MAT
Env_norm <- raster::subset(Env_norm, c(
                           #"Aridity", 
                           "MAP", "MAP_Seas", "MAT", 
                           "MAT_Seas", #"Snow", 
                            "Agriculture", "Dist_Urban", 
                           "Forest_Coni", "Forest_Deci", "NDVI", "Pastures", 
                            "Pop_Dens", "Shrubland", "Aspect", "Dist_Coast", 
                            "Dist_River", "Elev", "Slope", "CEC", "Clay.Silt", 
                            "Cu", "Hg", "Moisture", "N", "P", "pH", "SOC" #,"SoilT"
))

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate variable inflation factor (VIF) ####
# VIF is the extent of correlation between one predictor and all others.
# The lower VIF, the better we can tell what predictor contributed (most) to the model

## VIF basen on raw data (explanatory raster stack)
env_vif <- usdm::vif(Env_norm)

# which predictors should be excluded?
vif_cor <- usdm::vifcor(Env_norm, th=0.8)  #th = threshold vif for exclusion
# how: first find a pair of variables which has the maximum linear correlation 
# (greater than th), and exclude one of them which has greater VIF. The 
# procedure is repeated untill no variable with a high corrrelation coefficient 
# (grater than threshold) with other variables remains.

# merge both data.frames
env_vif <- env_vif %>% rename("VIF_raw" = VIF) %>% full_join(vif_cor@results) %>%
  full_join(as.data.frame(vif_cor@corMatrix) %>% mutate("Variables"=rownames(vif_cor@corMatrix)))

env_vif

write.csv(env_vif, file=paste0(here::here(), "/results/VIF_predictors.csv"), row.names = F)

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate correlations between predictors ####
# https://github.com/joaofgoncalves/GoncalvesAna_et_al_2021/tree/master/RCODE/PostModelAnalyses
# load predictors as dataframe
load(paste0(here::here(),"/results/EnvPredictor_2km_df_normalized.RData")) #Env_norm_df

# remove excluded variables
Env_norm_df <- Env_norm_df[,c("x", "y", names(Env_norm))] 

corMatSpearman <- cor(Env_norm_df, use="complete.obs", method="spearman") %>% round(2)
corMatPearson <- cor(Env_norm_df, use="complete.obs", method="pearson") %>% round(2)

write.csv(corMatSpearman, paste0("./results/corMatSpearman_predictors.csv"), row.names = F)
write.csv(corMatPearson, paste0("./results/corMatPearson_predictors.csv"), row.names = F)

rm(env_vif, vif_cor, corMatPearson, corMatSpearman, Env_norm_df)


#- - - - - - - - - - - - - - - - - - - - - -
# 5km #
#- - - - - - - - - - - - - - - - - - - - - -
# load environmental space
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_5km_normalized.grd"))

# remove SoilTemp and Aridity to remove vif of MAT
Env_norm <- raster::subset(Env_norm, c(
  #"Aridity", 
  "MAP", "MAP_Seas", "MAT", 
  "MAT_Seas", #"Snow", 
  "Agriculture", "Dist_Urban", 
  "Forest_Coni", "Forest_Deci", "NDVI", "Pastures", 
  "Pop_Dens", "Shrubland", "Aspect", "Dist_Coast", 
  "Dist_River", "Elev", "Slope", "CEC", "Clay.Silt", 
  "Cu", "Hg", "Moisture", "N", "P", "pH", "SOC" #,"SoilT"
))

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate variable inflation factor (VIF) ####
# VIF is the extent of correlation between one predictor and all others.
# The lower VIF, the better we can tell what predictor contributed (most) to the model

## VIF basen on raw data (explanatory raster stack)
env_vif <- usdm::vif(Env_norm)

# which predictors should be excluded?
vif_cor <- usdm::vifcor(Env_norm, th=0.8)  #th = threshold vif for exclusion
# how: first find a pair of variables which has the maximum linear correlation 
# (greater than th), and exclude one of them which has greater VIF. The 
# procedure is repeated untill no variable with a high corrrelation coefficient 
# (grater than threshold) with other variables remains.

# merge both data.frames
env_vif <- env_vif %>% rename("VIF_raw" = VIF) %>% full_join(vif_cor@results) %>%
  full_join(as.data.frame(vif_cor@corMatrix) %>% mutate("Variables"=rownames(vif_cor@corMatrix)))

env_vif

write.csv(env_vif, file=paste0(here::here(), "/results/VIF_predictors_5km.csv"), row.names = F)

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate correlations between predictors ####
# https://github.com/joaofgoncalves/GoncalvesAna_et_al_2021/tree/master/RCODE/PostModelAnalyses
# load predictors as dataframe
load(paste0(here::here(),"/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df

# remove excluded variables
Env_norm_df <- Env_norm_df[,c("x", "y", names(Env_norm))] 

corMatSpearman <- cor(Env_norm_df, use="complete.obs", method="spearman") %>% round(2)
corMatPearson <- cor(Env_norm_df, use="complete.obs", method="pearson") %>% round(2)

write.csv(corMatSpearman, paste0("./results/corMatSpearman_predictors_5km.csv"), row.names = F)
write.csv(corMatPearson, paste0("./results/corMatPearson_predictors_5km.csv"), row.names = F)

rm(env_vif, vif_cor, corMatPearson, corMatSpearman, Env_norm_df)

