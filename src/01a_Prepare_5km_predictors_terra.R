#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare predictor variables          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(terra)
library(tidyverse)

# change temporary directory
terraOptions(tempdir =  "D:/00_datasets/Trash")

# load 5km grid
grid5k <- terra::rast("D:/00_datasets/Grids/grid_5k_0p041.tif")

# load the predictor table containing the individual file names
pred_tab <- readr::read_csv(file="D:/_students/Romy/SoilBiodiversity/doc/Env_Predictors_table.csv")

# combine ID and Predictor name (to get folder names later on)
pred_tab$Predictor_long <- pred_tab$Predictor
pred_tab[!is.na(pred_tab$ID),]$Predictor_long <- paste0(pred_tab[!is.na(pred_tab$ID),]$ID, "_", pred_tab[!is.na(pred_tab$ID),]$Predictor)

# list the names of the variable's folders that will be included in the analysis
folders <- c(list.dirs("D:/00_datasets/Climate", recursive=F))
folders <- c(folders, list.dirs("D:/00_datasets/LandCover", recursive=F)) 
folders <- c(folders, list.dirs("D:/00_datasets/Location", recursive=F))
folders <- c(folders, list.dirs("D:/00_datasets/Soil", recursive=F))

#folders <- folders[stringr::str_detect(folders,"V[:digit:]{3}_")]
folders

#- - - - - - - - - - - - - - - - - - - - - 
## First, we need to use zonal statistics to project all variables to the 1km
# grid. This didn't work in ArcGIS because of too many unique ID values.

# Note: All raster files need to be in WGS1984 reference system.
#       Use ArcGIS tool "Project Raster" if necessary.

# # prepare for parallal processing
# cluster <- parallel::makeCluster(8)
# doParallel::registerDoParallel(cluster)

# loop through the variable's folders
#foreach(i = 1:length(folders), .packages = c("raster")) %dopar% { try({

makeTo5kmGrid <- function(raster_grid, temp_path, temp_file=NULL, file_name=NULL){
  
  # makeTo5kmGrid takes the temp_file to re-project it into the same extent and resolution as raster_grid.
  
  # raster_grid:  RasterLayer of a grid.
  # temp_path:    Path to the input file and at the same time output directory.
  # temp_file:    Name of the input file that will be re-projected.
  # file_name:    Name of the output file.
  
  print("===========================================")
  
  # define variable name (= folder name)
  temp_pred <- stringr::str_split_fixed(temp_path,"/",4)[4]
  
  # extract file name
  if(is.null(temp_file)) {
    temp_file <- as.character(pred_tab[pred_tab$Predictor_long==temp_pred, "File_name_processed"])
    
    # in case that 1km and 5km have differently processed files:
    if(stringr::str_detect(temp_file, "1km")==TRUE) { 
      temp_file <- stringr::str_replace(temp_file, "1km", "5km")}
  }
  
  print(paste0("Load raster called ", temp_file, "."))
  
  # read tif files with raw predictor variables across Europe
  temp_raster <- terra::rast(paste0(temp_path, "/", temp_file))
  
  print("Average per grid cell(s).")
  
  # calculate average per grid cell
  #temp_raster_mean <- terra::resample(temp_raster, raster_grid, method="average")
  
  # project to same crs
  temp_raster <- terra::project(temp_raster, grid5k, mask=TRUE)
  
  print("Mask raster.")
  
  # mask to same extent (extract out all the cells within the defined extent)
  #temp_raster_mean <- terra::mask(temp_raster_mean, raster_grid)
  
  # extract short predictor's name
  temp_name <- as.character(pred_tab[pred_tab$Predictor_long==temp_pred, "Predictor"])
  
  if(is.null(file_name) ){
    file_name <- paste0(temp_name, "_5km_mean.tif")
  }
  
  print("Save raster.")
  
  # save raster
  terra::writeRaster(temp_raster, file=paste0(temp_path, "/",file_name), overwrite=T)
  
  print(paste0("Raster called ", temp_path, "/", file_name, " saved." ))
  
  rm(temp_raster)
  temp_file <- NULL
  
}

for(i in 1:length(folders)){ try({
  makeTo5kmGrid(temp_path=folders[i], raster_grid=grid5k)
})}

# Note: Error in some variables will be solved individually

#- - - - - - - - - - - - - - - - - - - - - 
## Calculate some missing 5km grids ####

## Check what has been calculated
files <- list.files(folders, include.dirs = F, recursive=F)
files[stringr::str_detect(files, "_5km_mean.tif$")]

# Latitude
temp_raster <- as.data.frame(grid5k, xy=TRUE, row.names = FALSE)
temp_raster$Latitude <- temp_raster$y
temp_raster_mean <- terra::rast(temp_raster %>% dplyr::select(-grid_5k_0p041), type="xyz")
raster::writeRaster(temp_raster_mean, file="D:/00_datasets/Location/V045_Lat/Lat_5km_mean.tif", overwrite=T)

# Clay+Silt
temp_clay <- terra::rast("D:/00_datasets/Soil/Clay/Clay_5km_mean.tif")
temp_silt <- terra::rast("D:/00_datasets/Soil/Silt/Silt_5km_mean.tif")

temp_stack <- c(temp_clay, temp_silt)
temp_raster <- sum(temp_stack)
raster::writeRaster(temp_raster, file="D:/00_datasets/Soil/V062_Clay+Silt/Clay+Silt_5km_mean.tif", overwrite=T)

## Fix land cover proportions
# Agriculture is fine
# Forest

makeTo5kmGrid(temp_path = folders[22], temp_file="Forest_2012_noGrid_reclassified.tif")

makeTo5kmGrid(temp_path = folders[40], temp_file="Forest_Coni_2012_noGrid_reclassified.tif")

makeTo5kmGrid(temp_path = folders[41], temp_file="Forest_Deci_2012_noGrid_reclassified.tif")

# Pasture
makeTo5kmGrid(temp_path = folders[44], temp_file="Pasture_2012_noGrid_reclassified.tif")

# Shrubland
makeTo5kmGrid(temp_path = folders[46], temp_file="Shrubland_2012_noGrid_reclassified.tif")
## land cover done

# Snow (MODIS) for 2000-2009
makeTo5kmGrid(temp_path="D:/00_datasets/Climate/Snow_MODIS", raster_grid = grid5k, temp_file = "Snow_2000-2009_noGrid_mean.tif",
              file_name="Snow_2000-2009_5km_mean.tif")

# Dist_Roads (from shapefile)
# with ArcGIS.... Conversion -> To Raster -> Output cell size = CLC raw file from 2012
# Eucleadian distance...

# SoilTemp
makeTo5kmGrid(temp_path = "D:/00_datasets/Soil/V070_SoilT", raster_grid = grid5k, temp_file = "SBIO1_Annual_Mean_Temperature_5_15cm.tif", 
              file_name="SoilT_5-15cm_5km_mean.tif")

# OCTOP (OC in topsoil from FAO): missing reference system
makeTo5kmGrid(temp_path = "D:/00_datasets/Soil/SOC_OCTOP", raster_grid = grid5k, temp_file = "octop_V121_WGS84.tif", 
              file_name="SOC_OCTOP_5km_mean.tif")

#- - - - - - - - - - - - - - - - - -
## Copy files into one folder

# create folder
#dir.create("D:/_students/Romy/SoilBiodiversity/data_environment")

# select files
files <- list.files(folders, include.dirs = F, recursive=F, full.names = T)
stack_files <- files[stringr::str_detect(files, "_5km_mean.tif$")]

# keep only V0XX variables
stack_files <- stack_files[stringr::str_detect(stack_files, "V[:digit:]{3}")]

# remove Forest (but keep Forest_Deci and Forest_Coni)
stack_files <- stack_files[!stringr::str_detect(stack_files, "Forest_2012_")]
# remove SoilT at 5-15cm and older Snow data
stack_files <- stack_files[!stringr::str_detect(stack_files, "SoilT_5-15cm_")]
stack_files <- stack_files[!stringr::str_detect(stack_files, "Snow_2000-2009")]
stack_files
# remove Latitude
stack_files <- stack_files[!stringr::str_detect(stack_files, "Lat_")]

# copy files in folder
for(i in 1:length(stack_files)){
  temp_raster <- terra::rast(stack_files[i])
  temp_name <- basename(stack_files[i])
  
  terra::writeRaster(temp_raster, paste0("D:/_students/Romy/SoilBiodiversity/data_environment/", temp_name),
                     overwrite=TRUE)
}


#- - - - - - - - - - - - - - - - - - - - -
## Create future climate stacks ####
#- - - - - - - - - - - - - - - - - - - - -

# define future scenarios
scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
                               "mri-esm2-0", "ukesm1-0-ll"), "_",
                             rep(c("ssp126", "ssp370", "ssp585"),5)))

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

for(no_future in futureNames[str_detect(futureNames, "2041-2070")]){try({
  temp_files <- files[stringr::str_detect(files, no_future)]
  temp_files
  
  for(i in 1:length(temp_files)){try({
    
    # define name of variable
    temp_name <- stringr::str_extract(temp_files[i], "V[:digit:]{3}_[:alpha:]*_*[:alpha:]*")
    temp_name <- substr(temp_name, 6, nchar(temp_name))
    
    makeTo5kmGrid(temp_path = "D:/00_datasets/Climate", raster_grid = grid5k, temp_file = substr(temp_files[i], 24, nchar(temp_files[i])), 
                  file_name=paste0(temp_name, "_", no_future, "_5km_mean.tif"))
  })}
})}


files <- list.files(clim_folders, include.dirs = T, recursive=F, full.names = T)
files <- files[stringr::str_detect(files, "ssp[:digit:]*_5km_mean.tif$")]
files

# check if all exisit
check_files <- paste0(c("V002_MAP/Future/MAP", "V003_MAP_Seas/Future/MAP_Seas", "V004_MAT/Future/MAT", "V005_MAT_Seas/Future/MAT_Seas"), "_", rep(futureNames, 4))
check_files <- paste0("D:/00_datasets/Climate/", check_files, "_5km_mean.tif")
check_files

setdiff(sort(check_files), sort(files))

## Move files to one folder

# create folder
dir.create("D:/_students/Romy/SoilBiodiversity/data_environment/future_climate")

# copy files in folder
for(i in 1:length(files)){
  temp_raster <- terra::rast(files[i])
  temp_name <- basename(files[i])
  
  terra::writeRaster(temp_raster, paste0("D:/_students/Romy/SoilBiodiversity/data_environment/future_climate/", temp_name),
                     overwrite=TRUE)
}

