#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#        Define extent of predictors        #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(raster)

#- - - - - - - - - - - - - - - - - - - - - 
## Select uncorrelated variables ####
#- - - - - - - - - - - - - - - - - - - - -

# covariates
corMatPearson <- as.matrix(read.csv(file=paste0(here::here(), "/results/corMatPearson_predictors.csv")))
dimnames(corMatPearson)[[1]] <- dimnames(corMatPearson)[[2]]
# based on Valavi et al. 2021: Pearson 0.8
env_exclude <- caret::findCorrelation(corMatPearson, cutoff = 0.8, names=TRUE)
covarsNames <- dimnames(corMatPearson)[[1]][!(dimnames(corMatPearson)[[1]] %in% env_exclude)]
covarsNames <- covarsNames[covarsNames != "x" & covarsNames != "y"]
# exclude based on VIF
env_vif <- read.csv(file=paste0(here::here(), "/results/VIF_predictors.csv"))
env_exclude <- env_vif %>% filter(is.na(VIF)) %>% dplyr::select(Variables) %>% as.character()
covarsNames <- covarsNames[!(covarsNames %in% env_exclude)]
# excluded:
print("=== We excluded the following variables based on VIF and Pearson correlation: ===")
setdiff(env_vif$Variables, covarsNames)

# final predictor variables
print("=== And we kept the following, final predictor variables: ===")
covarsNames

#- - - - - - - - - - - - - - - - - - - - - 
## Clip 2km predictors ####
#- - - - - - - - - - - - - - - - - - - - -

# load environmental variables
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km_normalized.grd"))

Env_norm <- raster::subset(Env_norm, covarsNames)

# clip to exclude all cells with at least one layer being NA
Env_clip <- raster::mask(Env_norm, raster::calc(Env_norm, fun = sum))

# save Env_norm
raster::writeRaster(Env_clip, file=paste0(here::here(), "/results/EnvPredictor_2km_clipped.grd"), overwrite=T)

# same for dataframe
Env_clip <- terra::rast(paste0(here::here(), "/results/EnvPredictor_2km_clipped.grd"))
Env_clip_df <- as.data.frame(Env_clip, xy=TRUE, row.names=FALSE)
save(Env_clip_df, file=paste0(here::here(), "/results/EnvPredictor_2km_df_clipped.RData"))
rm(Env_clip, Env_norm, Env_clip_df)


#- - - - - - - - - - - - - - - - - - - - - 
## Clip 5km predictors ####
#- - - - - - - - - - - - - - - - - - - - -

# load environmental variables
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_5km_normalized.grd"))

Env_norm <- raster::subset(Env_norm, covarsNames)

# clip to exclude all cells with at least one layer being NA
Env_clip <- raster::mask(Env_norm, raster::calc(Env_norm, fun = sum))

# save Env_norm
raster::writeRaster(Env_clip, file=paste0(here::here(), "/results/EnvPredictor_5km_clipped.grd"), overwrite=T)

# same for dataframe
Env_clip <- terra::rast(paste0(here::here(), "/results/EnvPredictor_5km_clipped.grd"))
Env_clip_df <- as.data.frame(Env_clip, xy=TRUE, row.names=FALSE)
save(Env_clip_df, file=paste0(here::here(), "/results/EnvPredictor_5km_df_clipped.RData"))
rm(Env_clip, Env_norm, Env_clip_df)

