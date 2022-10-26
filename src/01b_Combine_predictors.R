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
library(tidyverse)
library(here)

# change temporary directory
terraOptions(tempdir =  paste0(here::here(), "/Trash"))

# load 2km grid
grid2k <- terra::rast(paste0(here::here(), "/data_environment/grid_2k_0p016.tif"))

#- - - - - - - - - - - - - - - - - - - - - 
## Merge raster files, 2km ####

# get names of the 2km files
files <- list.files(paste0(here::here(), "/data_environment"), 
                    include.dirs = F, recursive=F, full.names = T)

# create empty stack
Env <- raster::stack()

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
Env <- terra::rast("D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_2km.grd")
Env_df <- as.data.frame(Env, xyz=TRUE, row.names=FALSE)
save(Env_df, file="D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_2km_df.RData")

## cut to grid (should not be necessary anymore...)
#Env <- raster::mask(Env, grid1k)
# 
# # trim unnecessary margins
# Env_trimmed <- trim(Env, values = NA)
# Env_trimmed[is.na(Env_trimmed)] <- 0
# 
# Env <- raster::stack(Env)

#raster::plot(Env, maxnl=35)
dev.off()

#- - - - - - - - - - - - - - - - - - - - -
## Scale predictors ####
# scale rasterStack
Env_norm <- raster::scale(Env)
Env_norm <- raster::stack(Env_norm)

# save Env_norm
raster::writeRaster(Env_norm, file="D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_2km_normalized.grd", overwrite=T)

# same for dataframe
Env_norm_df <- as.data.frame(raster::rasterToPoints(Env_norm))
save(Env_norm_df, file="D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_2km_df_normalized.RData")

#- - - - - - - - - - - - - - - - - - - - -
## Create future climate stacks ####
#- - - - - - - - - - - - - - - - - - - - -

# same stack but MAT and MAP (MAT_Seas and MAP_Seas as well) from scenarios
futureNames <- sort(paste0(rep(c("2011-2040", "2041-2070", "2071-2100"),each=length(scenarioNames)), "_", scenarioNames))
futureNames

# define file names
clim_folders <- list.dirs("D:/00_datasets/Climate/V002_MAP", recursive=F)
clim_folders <- c(clim_folders, list.dirs("D:/00_datasets/Climate/V003_MAP_Seas", recursive=F))
clim_folders <- c(clim_folders, list.dirs("D:/00_datasets/Climate/V004_MAT", recursive=F))
clim_folders <- c(clim_folders, list.dirs("D:/00_datasets/Climate/V005_MAT_Seas", recursive=F))
clim_folders

files <- list.files(clim_folders, include.dirs = T, recursive=F, full.names = T)
files <- files[stringr::str_detect(files, "CHELSA")]
files

for(no_future in futureNames){try({
  temp_files <- files[stringr::str_detect(files, no_future)]
  temp_files
  
  for(i in 1:length(temp_files)){try({
    
    # define name of variable
    temp_name <- stringr::str_extract(temp_files[i], "V[:digit:]{3}_[:alpha:]*_*[:alpha:]*")
    temp_name <- substr(temp_name, 6, nchar(temp_name))
    
    makeTo2kmGrid(temp_path = "D:/00_datasets/Climate", raster_grid = grid2k, temp_file = substr(temp_files[i], 24, nchar(temp_files[i])), 
                  file_name=paste0(temp_name, "_", no_future, "_2km_mean.tif"))
  })}
})}


files <- list.files(clim_folders, include.dirs = T, recursive=F, full.names = T)
files <- files[stringr::str_detect(files, "ssp[:digit:]*_2km_mean.tif$")]
files

# check if all exisit
check_files <- paste0(c("V002_MAP/Future/MAP", "V003_MAP_Seas/Future/MAP_Seas", "V004_MAT/Future/MAT", "V005_MAT_Seas/Future/MAT_Seas"), "_", rep(futureNames, 4))
check_files <- paste0("D:/00_datasets/Climate/", check_files, "_2km_mean.tif")
check_files

setdiff(sort(check_files), sort(files))

# load env. stack
Env_norm <- raster::stack("D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_2km_normalized.grd")

for(no_future in futureNames){
  
  print("====================================")
  print(paste0("Processing of future scenario ", no_future))
  
  temp_files <- files[stringr::str_detect(files, no_future)]
  
  temp_Env_norm <- Env_norm
  
  if(length(temp_files)!=4) print("Please check... there are not 4 climate variables")
  
  for(i in 1:length(temp_files)){
    
    # define name of variable
    temp_name <- stringr::str_extract(temp_files[i], "V[:digit:]{3}_[:alpha:]*_*[:alpha:]*")
    temp_name <- substr(temp_name, 6, nchar(temp_name))
    
    temp_raster <- raster::raster(temp_files[i])
    names(temp_raster) <- temp_name
    
    temp_raster <- raster::mask(temp_raster, grid2k)
    temp_raster <- raster::scale(temp_raster)
    
    temp_Env_norm[[temp_name]] <- temp_raster
    
    print(paste0("Stacked file ", names(temp_raster)))
  }
  
  # save Env_norm
  raster::writeRaster(temp_Env_norm, file=paste0("D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_", no_future, "_2km_normalized.grd"), overwrite=T)
  
  # same for dataframe
  temp_Env_df <- as.data.frame(raster::rasterToPoints(temp_Env_norm))
  save(temp_Env_df, file=paste0("D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_", no_future, "_2km_df_normalized.RData"))
  
}



#- - - - - - - - - - - - - - - - - - - - -
## Save all ####
# save raster stack plots
files <- list.files(folders, include.dirs = F, recursive=F, full.names = T)
stack_files <- files[stringr::str_detect(files, "_2km_mean.tif$")]

stack_files

# create empty stack
Env <- raster::stack()

# load and merge them
for(i in 1:length(stack_files)){
  temp_raster <- raster::raster(stack_files[i])
  
  names(temp_raster) <- gsub("_2.*", "", names(temp_raster))
  
  temp_raster <- raster::mask(temp_raster, grid2k)
  
  Env <- raster::stack(Env, temp_raster)
  
  print(paste0("Stacked file ", names(temp_raster)))
}

# save raster
raster::writeRaster(Env, file="D:/00_datasets/EnvPredictor_2km_all.grd", overwrite=T)

#pdf(file="D:/00_datasets/Predictors_Europe_2km.pdf", height=15, width = 18)
raster::plot(Env, maxnl=70)
dev.off()

rm(Env)

