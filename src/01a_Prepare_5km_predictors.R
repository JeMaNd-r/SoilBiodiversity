#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare predictor variables          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# change temporary directory for files
raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

# load 5km grid
grid5k <- raster::raster("D:/00_datasets/Grids/grid_5k_0p041.tif")
raster::crs(grid5k) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

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
  #temp_pred <- stringr::str_extract(temp_path, "V[:digit:]{3}_[:alpha:]*_*[:alpha:]*")
  
  # extract file name
  if(is.null(temp_file)) {
    temp_file <- as.character(pred_tab[pred_tab$Predictor_long==temp_pred, "File_name_processed"])
    
    # in case that 1km and 5km have differencly processed files:
    if(stringr::str_detect(temp_file, "1km")==TRUE) { 
      temp_file <- stringr::str_replace(temp_file, "1km", "5km")}
  }
  
  print(paste0("Load raster called ", temp_file, "."))
  
  # read tif files with raw predictor variables across Europe
  temp_raster <- raster::raster(paste0(temp_path, "/", temp_file))
  
  # add info on coordinate reference system
  #raster::crs(temp_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  ## no cropping needed anymore (?)
  # crop to same extent (extract out all the cells within the defined extent)
  #temp_raster <- raster::crop(temp_raster, raster::extent(raster_grid))
  #raster::extent(temp_raster) <- raster::extent(raster_grid) 
  
  print("Average per grid cell(s).")
  
  # calculate average per grid cell
  temp_raster_mean <- raster::resample(temp_raster, raster_grid)
  #temp_raster_mean <- raster::rasterize(rasterToPolygons(raster_grid), temp_raster, fun=mean)
  
  print("Mask raster.")
  
  # mask to same extent (extract out all the cells within the defined extent)
  temp_raster_mean <- raster::mask(temp_raster_mean, raster_grid)
  #raster::extent(temp_raster) <- raster::extent(raster_grid) 
  
  # extract short predictor's name
  temp_name <- as.character(pred_tab[pred_tab$Predictor_long==temp_pred, "Predictor"])
  
  #temp_raster_mean2 <- crop(temp_raster_mean, raster_grid)
  
  if( is.null(file_name) ){
    file_name <- paste0(temp_name, "_5km_mean.tif")
  }
  
  print("Save raster.")
  
  # save raster
  raster::writeRaster(temp_raster_mean, file=paste0(temp_path, "/",file_name), overwrite=T)
  
  print(paste0("Raster called ", temp_path, "/", file_name, " saved." ))
  
  rm(temp_raster, temp_raster_mean)
  temp_file <- NULL
  
}

#parallel::stopCluster(cluster)

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
temp_raster <- as.data.frame(raster::rasterToPoints(grid5k))
temp_raster$Latitude <- temp_raster$y
temp_raster_mean <- raster::rasterFromXYZ(temp_raster %>% dplyr::select(-grid_5k_0p041))
raster::writeRaster(temp_raster_mean, file="D:/00_datasets/Location/V019_Lat/Lat_5km_mean.tif", overwrite=T)

# Clay+Silt
temp_clay <- raster::raster("D:/00_datasets/Soil/Clay/Clay_5km_mean.tif")
temp_silt <- raster::raster("D:/00_datasets/Soil/Silt/Silt_5km_mean.tif")

temp_stack <- raster::stack(temp_clay, temp_silt)
temp_raster <- sum(temp_stack)
raster::writeRaster(temp_raster, file="D:/00_datasets/Soil/V062_Clay+Silt/Clay+Silt_5km_mean.tif", overwrite=T)

# Snow (MODIS) for 2000-2009
makeTo5kmGrid(temp_path="D:/00_datasets/Climate/Snow_MODIS", raster_grid = grid5k, temp_file = "Snow_2000-2009_noGrid_mean.tif",
           file_name="Snow_2000-2009_5km_mean.tif")

# Dist_Roads (from shapefile)
# with ArcGIS.... Conversion -> To Raster -> Output cell size = CLC raw file from 2012
# Eucleadian distance...

# SoilTemp
makeTo5kmGrid(temp_path = "D:/00_datasets/Soil/V070_SoilT", raster_grid = grid5k, temp_file = "SBIO1_Annual_Mean_Temperature_5_15cm.tif", 
           file_name="SoilT_5-15cm_5km_mean.tif")

# # SOC_OCTOP #not fixed yet...
# temp_raster <- raster::raster("D:/00_datasets/Soil/SOC_OCTOP/octop_V121.asc")
# temp_raster_mean <- raster::resample(temp_raster, grid5k)
# temp_raster_mean

# # future climate
# makeTo5kmGrid(temp_path = "D:/00_datasets/Climate/V002_MAP/Future", raster_grid = grid5k, temp_file = "CHELSA_bio12_2011-2040_gfdl-esm4_ssp126_V.2.1.tif", 
#               file_name="MAP_Future_5km_mean.tif")
# makeTo5kmGrid(temp_path = "D:/00_datasets/Climate/V003_MAP_Seas/Future", raster_grid = grid5k, temp_file = "CHELSA_bio15_2011-2040_gfdl-esm4_ssp126_V.2.1.tif", 
#               file_name="MAP_Seas_Future_5km_mean.tif")
# makeTo5kmGrid(temp_path = "D:/00_datasets/Climate/V004_MAT/Future", raster_grid = grid5k, temp_file = "CHELSA_bio1_2011-2040_gfdl-esm4_ssp126_V.2.1.tif", 
#               file_name="MAT_Future_5km_mean.tif")
# makeTo5kmGrid(temp_path = "D:/00_datasets/Climate/V005_MAT_Seas/Future", raster_grid = grid5k, temp_file = "CHELSA_bio4_2011-2040_gfdl-esm4_ssp126_V.2.1.tif", 
#               file_name="MAT_Seas_Future_5km_mean.tif")

#- - - - - - - - - - - - - - - - - - - - - 
## Mask selected 5km grids ####

# load mask file (one of CORINE land cover files)
temp_mask <- raster::raster("D:/00_datasets/LandCover/V021_Agriculture/Agriculture_5km_mean.tif")

# Aspect: no data for Ukraine
temp_raster <- raster::raster("D:/00_datasets/Location/V041_Aspect/Aspect_5km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/Location/V041_Aspect/Aspect_5km_mean.tif", overwrite=T)

# Dist_Urban: no (CORINE) data for Ukraine
temp_raster <- raster::raster("D:/00_datasets/LandCover/V022_Dist_Urban/Dist_Urban_5km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/LandCover/V022_Dist_Urban/Dist_Urban_5km_mean.tif", overwrite=T)

# Dist_River
temp_raster <- raster::raster("D:/00_datasets/Location/V043_Dist_River/Dist_River_5km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/Location/V043_Dist_River/Dist_River_5km_mean.tif", overwrite=T)

# Dist_Coast
temp_raster <- raster::raster("D:/00_datasets/Location/V042_Dist_Coast/Dist_Coast_5km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/Location/V042_Dist_Coast/Dist_Coast_5km_mean.tif", overwrite=T)

# Dist_Roads
temp_raster <- raster::raster("D:/00_datasets/LandCover/Dist_Roads/Dist_Roads_5km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/LandCover/Dist_Roads/Dist_Roads_5km_mean.tif", overwrite=T)

# Impervious
temp_raster <- raster::raster("D:/00_datasets/LandCover/Impervious/Impervious_5km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/LandCover/Impervious/Impervious_5km_mean.tif", overwrite=T)

rm(temp_raster, temp_mask)

#- - - - - - - - - - - - - - - - - - - - - 
## Merge all raster files, 1km ####
# get names of the 1km files
files <- list.files(folders, include.dirs = F, recursive=F, full.names = T)
stack_files <- files[stringr::str_detect(files, "_5km_mean.tif$")]

# keep only V0XX variables
stack_files <- stack_files[stringr::str_detect(stack_files, "V[:digit:]{3}")]

# remove Forest (but keep Forest_Deci and Forest_Coni)
stack_files <- stack_files[!stringr::str_detect(stack_files, "Forest_2012_")]
# remove SoilT at 5-15cm
stack_files <- stack_files[!stringr::str_detect(stack_files, "SoilT_5-15cm_")]
stack_files <- stack_files[!stringr::str_detect(stack_files, "Snow_2000-2009")]
stack_files
# remove Latitude
stack_files <- stack_files[!stringr::str_detect(stack_files, "Lat_")]

# create empty stack
Env <- raster::stack()

# load and merge them
for(i in 1:length(stack_files)){
  temp_raster <- raster::raster(stack_files[i])
  
  names(temp_raster) <- gsub("_5.*", "", names(temp_raster))
  
  temp_raster <- raster::mask(temp_raster, grid5k)
  
  Env <- raster::stack(Env, temp_raster)
  
  print(paste0("Stacked file ", names(temp_raster)))
}

Env

# save raster
raster::writeRaster(Env, file="D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_5km.grd", overwrite=T)

# as dataframe
Env_df <- as.data.frame(raster::rasterToPoints(Env))
save(Env_df, file="D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_5km_df.RData")

## cut to grid (should not be necessary anymore...)
#Env <- raster::mask(Env, grid5k)
# 
# # trim unnecessary margins
# Env_trimmed <- trim(Env, values = NA)
# Env_trimmed[is.na(Env_trimmed)] <- 0
# 
# Env <- raster::stack(Env)

pdf(file="D:/_students/Romy/SoilBiodiversity/figures/Predictors_Europe_5km.pdf", height=15, width = 18)
raster::plot(Env, maxnl=35)
dev.off()

#- - - - - - - - - - - - - - - - - - - - -
## Scale predictors ####
# scale rasterStack
Env_norm <- raster::scale(Env)
Env_norm <- raster::stack(Env_norm)
  
# save Env_norm
raster::writeRaster(Env_norm, file="D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_5km_normalized.grd", overwrite=T)

# same for dataframe
Env_norm_df <- as.data.frame(raster::rasterToPoints(Env_norm))
save(Env_norm_df, file="D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_5km_df_normalized.RData")

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

# load env. stack
Env_norm <- raster::stack("D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_5km_normalized.grd")

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

      temp_raster <- raster::mask(temp_raster, grid5k)
      temp_raster <- raster::scale(temp_raster)
  
      temp_Env_norm[[temp_name]] <- temp_raster
  
      print(paste0("Stacked file ", names(temp_raster)))
   }

   # save Env_norm
   raster::writeRaster(temp_Env_norm, file=paste0("D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_", no_future, "_5km_normalized.grd"), overwrite=T)

   # same for dataframe
   temp_Env_df <- as.data.frame(raster::rasterToPoints(temp_Env_norm))
   save(temp_Env_df, file=paste0("D:/_students/Romy/SoilBiodiversity/results/EnvPredictor_", no_future, "_5km_df_normalized.RData"))
   
}



#- - - - - - - - - - - - - - - - - - - - -
## Save all ####
# save raster stack plots
files <- list.files(folders, include.dirs = F, recursive=F, full.names = T)
stack_files <- files[stringr::str_detect(files, "_5km_mean.tif$")]

stack_files

# create empty stack
Env <- raster::stack()

# load and merge them
for(i in 1:length(stack_files)){
  temp_raster <- raster::raster(stack_files[i])
  
  names(temp_raster) <- gsub("_5.*", "", names(temp_raster))
  
  temp_raster <- raster::mask(temp_raster, grid5k)
  
  Env <- raster::stack(Env, temp_raster)
  
  print(paste0("Stacked file ", names(temp_raster)))
}

# save raster
raster::writeRaster(Env, file="D:/00_datasets/EnvPredictor_5km_all.grd", overwrite=T)

pdf(file="D:/00_datasets/Predictors_Europe_5km.pdf", height=15, width = 18)
raster::plot(Env, maxnl=70)
dev.off()

rm(Env)
