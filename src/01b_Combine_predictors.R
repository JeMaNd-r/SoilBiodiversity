#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Combine predictor variables          #
#        into one raster stack              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")
gc()

library(terra)
library(raster)
library(tidyverse)
library(here)

# change temporary directory
terraOptions(tempdir =  paste0(here::here(), "/Trash"))

#- - - - - - - - - - - - - - - - - - - - - 
## Merge raster files, 2km ####
#- - - - - - - - - - - - - - - - - - - - -

# load 2km grid
grid2k <- raster::raster(paste0(here::here(), "/data_environment/grid_2k_0p016.tif"))

# get names of the 2km files
files <- list.files(paste0(here::here(), "/data_environment"), 
                    include.dirs = F, recursive=F, full.names = T)
files <- files[str_detect(files, "2km")]
files

# create empty stack
Env <- raster::stack() 
#unfortunately, "terra" uses too much memory while creating the stack

# load and merge them
for(i in 1:length(files)){
  temp_raster <- raster::raster(files[i])
  
  names(temp_raster) <- gsub("_2.*", "", basename(files[i]))
  
  temp_raster <- raster::mask(temp_raster, grid2k)
  
  Env <- raster::stack(Env, temp_raster)
  
  print(paste0("Stacked file ", names(temp_raster)))
}

Env

# save raster
raster::writeRaster(Env, file=paste0(here::here(), "/results/EnvPredictor_2km.grd"), overwrite=T)
# You may need to create the "results" folder first.

# as dataframe
Env <- terra::rast(paste0(here::here(), "/results/EnvPredictor_2km.grd"))
Env_df <- as.data.frame(Env, xy=TRUE, row.names=FALSE)
save(Env_df, file=paste0(here::here(), "/results/EnvPredictor_2km_df.RData"))
rm(Env_df)

## cut to grid (should not be necessary anymore...)
#Env <- raster::mask(Env, grid1k)
# 
# # trim unnecessary margins
# Env_trimmed <- trim(Env, values = NA)
# Env_trimmed[is.na(Env_trimmed)] <- 0
# 
# Env <- raster::stack(Env)

#raster::plot(Env, maxnl=35)
#dev.off()

#- - - - - - - - - - - - - - - - - - - - -
## Scale predictors ####
Env <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km.grd"))

# scale rasterStack
Env_norm <- raster::scale(Env)
Env_norm <- raster::stack(Env_norm)

# save Env_norm
raster::writeRaster(Env_norm, file=paste0(here::here(), "/results/EnvPredictor_2km_normalized.grd"), overwrite=T)

#- - - - - - - - - - - - - - - - - - - - - 
## Merge raster files, 5km ####
#- - - - - - - - - - - - - - - - - - - - -

# load 5km grid
grid5k <- raster::raster(paste0(here::here(), "/data_environment/grid_5k_0p041.tif"))

# get names of the 5km files
files <- list.files(paste0(here::here(), "/data_environment"), 
                    include.dirs = F, recursive=F, full.names = T)
files <- files[stringr::str_detect(files, "5km")]
files

# create empty stack
Env <- raster::stack() 
#unfortunately, "terra" uses too much memory while creating the stack

# load and merge them
for(i in 1:length(files)){
  temp_raster <- raster::raster(files[i])
  
  names(temp_raster) <- gsub("_5.*", "", basename(files[i]))
  
  temp_raster <- raster::mask(temp_raster, grid5k)
  
  Env <- raster::stack(Env, temp_raster)
  
  print(paste0("Stacked file ", names(temp_raster)))
}

Env

# save raster
raster::writeRaster(Env, file=paste0(here::here(), "/results/EnvPredictor_5km.grd"), overwrite=T)
# You may need to create the "results" folder first.

# as dataframe
Env <- terra::rast(paste0(here::here(), "/results/EnvPredictor_5km.grd"))
Env_df <- as.data.frame(Env, xy=TRUE, row.names=FALSE)
save(Env_df, file=paste0(here::here(), "/results/EnvPredictor_5km_df.RData"))
rm(Env_df)

#terra::plot(Env, maxnl=35)
#dev.off()

#- - - - - - - - - - - - - - - - - - - - -
## Scale predictors ####
Env <- raster::stack(paste0(here::here(), "/results/EnvPredictor_5km.grd"))

# scale rasterStack
Env_norm <- raster::scale(Env)
Env_norm <- raster::stack(Env_norm)

# save Env_norm
raster::writeRaster(Env_norm, file=paste0(here::here(), "/results/EnvPredictor_5km_normalized.grd"), overwrite=T)


