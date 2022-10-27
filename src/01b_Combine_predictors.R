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
Env_df <- as.data.frame(Env, xyz=TRUE, row.names=FALSE)
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

# same for dataframe
Env_norm <- terra::rast(paste0(here::here(), "/results/EnvPredictor_2km_normalized.grd"))
Env_norm_df <- as.data.frame(Env_norm, xyz=TRUE, row.names=FALSE)
save(Env_norm_df, file=paste0(here::here(), "/results/EnvPredictor_2km_df_normalized.RData"))
rm(Env, Env_norm, Env_norm_df)

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
Env_df <- as.data.frame(Env, xyz=TRUE, row.names=FALSE)
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

# same for dataframe
Env_norm <- terra::rast(paste0(here::here(), "/results/EnvPredictor_5km_normalized.grd"))
Env_norm_df <- as.data.frame(Env_norm, xyz=TRUE, row.names=FALSE)
save(Env_norm_df, file=paste0(here::here(), "/results/EnvPredictor_5km_df_normalized.RData"))
rm(Env, Env_norm, Env_norm_df)



## NOT NEEDED ####
# #- - - - - - - - - - - - - - - - - - - - -
# ## Create future climate stacks ####
# #- - - - - - - - - - - - - - - - - - - - -
# 
# # define future scenarios
# scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
#                                "mri-esm2-0", "ukesm1-0-ll"), "_",
#                              rep(c("ssp126", "ssp370", "ssp585"),5)))
# 
# # same stack but MAT and MAP (MAT_Seas and MAP_Seas as well) from scenarios
# futureNames <- sort(paste0(rep(c("2011-2040", "2041-2070", "2071-2100"),each=length(scenarioNames)), "_", scenarioNames))
# futureNames
# 
# files <- list.files(paste0(here::here(), "/results/EnvPredictor_5km_normalized.grd"), include.dirs = T, recursive=F, full.names = T)
# files <- files[stringr::str_detect(files, "ssp[:digit:]*_5km_mean.tif$")]
# files
# 
# 
# # load env. stack
# Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_5km_normalized.grd"))
# 
# for(no_future in futureNames){
#   
#   print("====================================")
#   print(paste0("Processing of future scenario ", no_future))
#   
#   temp_files <- files[stringr::str_detect(files, no_future)]
#   
#   temp_Env_norm <- Env_norm
#   
#   if(length(temp_files)!=4) print("Please check... there are not 4 climate variables")
#   
#   for(i in 1:length(temp_files)){
#     
#     # define name of variable
#     temp_name <- stringr::str_extract(temp_files[i], "V[:digit:]{3}_[:alpha:]*_*[:alpha:]*")
#     temp_name <- substr(temp_name, 6, nchar(temp_name))
#     
#     temp_raster <- raster::raster(temp_files[i])
#     names(temp_raster) <- temp_name
#     
#     temp_raster <- raster::mask(temp_raster, grid2k)
#     temp_raster <- raster::scale(temp_raster)
#     
#     temp_Env_norm[[temp_name]] <- temp_raster
#     
#     print(paste0("Stacked file ", names(temp_raster)))
#   }
#   
#   # save Env_norm
#   raster::writeRaster(temp_Env_norm, file=paste0(here::here(), "/results/EnvPredictor_", no_future, "_5km_normalized.grd"), overwrite=T)
#   
#   # same for dataframe
#   temp_Env_df <- as.data.frame(raster::rasterToPoints(temp_Env_norm))
#   save(temp_Env_df, file=paste0(here::here(), "/results/EnvPredictor_", no_future, "_5km_df_normalized.RData"))
#   
# }


