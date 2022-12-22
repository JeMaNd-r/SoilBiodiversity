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
library(terra) #faster as raster but can't handle big stacks (memmory issue)

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


#- - - - - - - - - - - - - - - - - - - - -
## Create future climate stacks ####
#- - - - - - - - - - - - - - - - - - - - -

# define future scenarios
scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr",
                               "mri-esm2-0", "ukesm1-0-ll"), "_",
                             rep(c("ssp126", "ssp370", "ssp585"),5)))

# same stack but MAT and MAP (MAT_Seas and MAP_Seas as well) from scenarios
futureNames <- sort(paste0(rep("2041-2070",each=length(scenarioNames)), "_", scenarioNames))
futureNames

files <- list.files(paste0(here::here(), "/data_environment/future_climate/"), full.names = T)
files <- files[stringr::str_detect(files, "ssp[:digit:]*_5km_mean.tif$")]
files


# load env. stack
Env_clip <- raster::stack(paste0(here::here(), "/results/EnvPredictor_5km_clipped.grd"))

# load scaling parameters
scale_values <- read.csv(file=paste0(here::here(), "/results/EnvPredictor_scaling.csv"))

for(no_future in futureNames){
  
  print("====================================")
  print(paste0("Processing of future scenario ", no_future))
  
  temp_files <- files[stringr::str_detect(files, no_future)]
  
  temp_Env <- Env_clip
  
  if(length(temp_files)!=4) print("Please check... there are not 4 climate variables")
  
  for(i in 1:length(temp_files)){
    
    # define name of variable
    temp_name <- stringr::str_extract(basename(temp_files[i]), "[:graph:]*_2")
    temp_name <- substr(temp_name, 1, nchar(temp_name)-2)
    
    temp_raster <- raster::raster(temp_files[i])
    names(temp_raster) <- temp_name
    
    temp_raster <- raster::mask(temp_raster, Env_clip[[1]])
    
    # scale each predictor by the respective scale values (based on current climate)
    temp_mean <- scale_values[scale_values$predictor==temp_name, "scale_mean"]
    temp_sd <- scale_values[scale_values$predictor==temp_name, "scale_sd"]
      
    temp_raster <- raster::calc(temp_raster, 
                        fun = function(x, na.rm=T){(x - temp_mean) / temp_sd},
                        na.rm=TRUE)
    names(temp_raster) <- temp_name
   
    #temp_raster <- terra::scale(temp_raster) #not used as we want to scale with same values as current climate
    
    temp_Env[[temp_name]] <- temp_raster
    
    print(paste0("Stacked file ", names(temp_raster)))
  }
  
  # same for dataframe
  temp_Env_df <- as.data.frame(temp_Env, xy=TRUE)
  save(temp_Env_df, file=paste0(here::here(), "/results/_FutureEnvironment/EnvPredictor_", no_future, "_5km_df_clipped.RData"))
  
}

