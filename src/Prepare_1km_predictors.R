#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare predictor variables          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# change temporary directory for files
raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

# load the 1km grid
grid1k <- raster::raster("D:/00_datasets/Grids/grid_1k_0p008.tif")
raster::crs(grid1k) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# load the predictor table containing the individual file names
pred_tab <- readr::read_csv(file="I:/eie/==PERSONAL/RZ SoilBON/SoilBiodiversity/doc/Env_Predictors_table.csv")

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
 
makeToGrid <- function(raster_grid, temp_path, temp_file=NULL, file_name=NULL){
  
  # makeToGrid takes the temp_file to re-project it into the same extent and resolution as raster_grid.
  
  print("===========================================")
  
  # raster_grid:  RasterLayer of a grid.
  # temp_path:    Path to the input file and at the same time output directory.
  # temp_file:    Name of the input file that will be re-projected.
  # file_name:    Name of the output file.
  
  # define variable name (= folder name)
  temp_pred <- stringr::str_split_fixed(temp_path,"/",4)[,4]
  #temp_pred <- stringr::str_extract(temp_path, "V[:digit:]{3}_[:alpha:]*_*[:alpha:]*")
  
  # extract file name
  if(is.null(temp_file)) {
    temp_file <- as.character(pred_tab[pred_tab$Predictor_long==temp_pred, "File_name_processed"])
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
    file_name <- paste0(temp_name, "_1km_mean.tif")
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
  makeToGrid(temp_path=folders[i], raster_grid=grid1k)
})}

# Note: Error in some variables will be solved individually

#- - - - - - - - - - - - - - - - - - - - - 
## Calculate some missing 1km grids ####

## Check what has been calculated
files <- list.files(folders, include.dirs = F, recursive=F)
files[stringr::str_detect(files, "_1km_mean.tif$")]

# SoilTemp
makeToGrid(temp_path = "D:/00_datasets/Soil/SoilT", raster_grid = grid1k, temp_file = "SBIO1_Annual_Mean_Temperature_5_15cm.tif", 
           file_name="SoilT_5-15cm_1km_mean.tif")

# Latitude
temp_raster <- as.data.frame(raster::rasterToPoints(grid1k))
temp_raster$Latitude <- temp_raster$y
temp_raster_mean <- rasterFromXYZ(temp_raster %>% dplyr::select(-grid_1k_0p008))
raster::writeRaster(temp_raster_mean, file="D:/00_datasets/Location/V019_Lat/Lat_1km_mean.tif", overwrite=T)

# Snow (MODIS) for 2000-2009
makeToGrid(temp_path="D:/00_datasets/Climate/Snow", raster_grid = grid1k, temp_file = "Snow_2000-2009_noGrid_mean.tif",
           file_name="Snow_2000-2009_1km_mean.tif")

# Dist_Roads (from shapefile)
# with ArcGIS.... Conversion -> To Raster -> Output cell size = CLC raw file from 2012
# Eucleadian distance...

# Clay+Silt
temp_clay <- raster::raster("D:/00_datasets/Soil/Clay/Clay_1km_mean.tif")
temp_silt <- raster::raster("D:/00_datasets/Soil/Silt/Silt_1km_mean.tif")

temp_stack <- raster::stack(temp_clay, temp_silt)
temp_raster <- sum(temp_stack)
raster::writeRaster(temp_raster, file="D:/00_datasets/Soil/V062_Clay+Silt/Clay+Silt_1km_mean.tif", overwrite=T)

#- - - - - - - - - - - - - - - - - - - - - 
## Mask selected 1km grids ####

# load mask file (one of CORINE land cover files)
temp_mask <- raster::raster("D:/00_datasets/LandCover/V021_Agriculture/Agriculture_1km_mean.tif")

# Aspect: no data for Ukraine
temp_raster <- raster::raster("D:/00_datasets/Location/V015_Aspect/Aspect_1km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/Location/V015_Aspect/Aspect_1km_mean.tif", overwrite=T)

# Dist_Urban: no (CORINE) data for Ukraine
temp_raster <- raster::raster("D:/00_datasets/LandCover/V022_Dist_Urban/Dist_Urban_1km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/LandCover/V022_Dist_Urban/Dist_Urban_1km_mean.tif", overwrite=T)

# Dist_River
temp_raster <- raster::raster("D:/00_datasets/Location/V043_Dist_River/Dist_River_1km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/Location/V043_Dist_River/Dist_River_1km_mean.tif", overwrite=T)

# Dist_Coast
temp_raster <- raster::raster("D:/00_datasets/Location/V042_Dist_Coast/Dist_Coast_1km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/Location/V042_Dist_Coast/Dist_Coast_1km_mean.tif", overwrite=T)

# Dist_Roads
temp_raster <- raster::raster("D:/00_datasets/LandCover/Dist_Roads/Dist_Roads_1km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/LandCover/Dist_Roads/Dist_Roads_1km_mean.tif", overwrite=T)

# Impervious
temp_raster <- raster::raster("D:/00_datasets/LandCover/Impervious/Impervious_1km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/LandCover/Impervious/Impervious_1km_mean.tif", overwrite=T)

rm(temp_raster, temp_mask)

#- - - - - - - - - - - - - - - - - - - - - 
## Merge all raster files, 1km ####
# get names of the 1km files
files <- list.files(folders, include.dirs = F, recursive=F, full.names = T)
stack_files <- files[stringr::str_detect(files, "_1km_mean.tif$")]

# keep only V0XX variables
stack_files <- stack_files[stringr::str_detect(stack_files, "V[:digit:]{3}")]

# remove Forest (but keep Forest_Deci and Forest_Coni)
stack_files <- stack_files[!stringr::str_detect(stack_files, "Forest_2012_")]
# remove SoilT at 5-15cm
stack_files <- stack_files[!stringr::str_detect(stack_files, "SoilT_5-15cm_")]
stack_files <- stack_files[!stringr::str_detect(stack_files, "Snow_2000-2009")]
stack_files

# create empty stacj
Env <- raster::stack()

# load and merge them
for(i in 1:length(stack_files)){
  temp_raster <- raster::raster(stack_files[i])
  
  names(temp_raster) <- gsub("_1.*", "", names(temp_raster))
  
  temp_raster <- raster::mask(temp_raster, grid1k)
  
  Env <- raster::stack(Env, temp_raster)
  
  print(paste0("Stacked file ", names(temp_raster)))
}

# save raster
raster::writeRaster(Env, file="I:/eie/==PERSONAL/RZ SoilBON/SoilBiodiversity/results/EnvPredictor_1km.grd", overwrite=T)

# as dataframe
Env_df <- as.data.frame(raster::rasterToPoints(Env))
save(Env_df, file="I:/eie/==PERSONAL/RZ SoilBON/SoilBiodiversity/results/EnvPredictor_1km_df.RData")

## cut to grid (should not be necessary anymore...)
#Env <- raster::mask(Env, grid1k)
# 
# # trim unnecessary margins
# Env_trimmed <- trim(Env, values = NA)
# Env_trimmed[is.na(Env_trimmed)] <- 0
# 
# Env <- raster::stack(Env)

#pdf(file="I:/eie/==PERSONAL/RZ SoilBON/SoilBiodiversity/figures/Predictors_Europe_1km.pdf", height=15, width = 18)
raster::plot(Env, maxnl=35)
dev.off()

#- - - - - - - - - - - - - - - - - - - - -
## Scale predictors ####
# scale rasterStack
Env_norm <- raster::scale(Env)
Env_norm <- raster::stack(Env_norm)

# save Env_norm
raster::writeRaster(Env_norm, file="I:/eie/==PERSONAL/RZ SoilBON/SoilBiodiversity/results/EnvPredictor_1km_normalized.grd", overwrite=T)

# same for dataframe
Env_norm_df <- as.data.frame(raster::rasterToPoints(Env_norm))
save(Env_norm_df, file="I:/eie/==PERSONAL/RZ SoilBON/SoilBiodiversity/results/EnvPredictor_1km_df_normalized.RData")

#- - - - - - - - - - - - - - - - - - - - -
## Save all ####
# save raster stack plots
files <- list.files(folders, include.dirs = F, recursive=F, full.names = T)
stack_files <- files[stringr::str_detect(files, "_1km_mean.tif$")]

stack_files

# create empty stacj
Env <- raster::stack()

# load and merge them
for(i in 1:length(stack_files)){
  temp_raster <- raster::raster(stack_files[i])
  
  names(temp_raster) <- gsub("_1.*", "", names(temp_raster))
  
  temp_raster <- raster::mask(temp_raster, grid1k)
  
  Env <- raster::stack(Env, temp_raster)
  
  print(paste0("Stacked file ", names(temp_raster)))
}

# save raster
raster::writeRaster(Env, file="D:/00_datasets/EnvPredictor_1km_all.grd", overwrite=T)

#pdf(file="D:/00_datasets/Predictors_Europe_1km.pdf", height=15, width = 18)
raster::plot(Env, maxnl=70)
dev.off()

# #- - - - - - - - - - - - - - - - - - - - - 
# ## OLD CODE for fake analyses ####
# #- - - - - - - - - - - - - - - - - - - - - 
# # get Worldclim raster data layers
# env.bioclim <- raster::getData('worldclim', var='bio', res=10)
# 
# #summary(env.bioclim[[1]])
# plot(env.bioclim[[1]], main="BioClim 1 (raw)")
# 
# ## Correlations
# env.cor <- sdmpredictors::pearson_correlation_matrix(env.bioclim)
# 
# # plot correlations
# corrplot::corrplot.mixed(env.cor)
# 
# pdf(file=paste0("Predictor_correlation_", Taxon_name, ".pdf"))
# sdmpredictors::plot_correlation(env.cor)
# dev.off()
# 
# sdmpredictors::correlation_groups(env.cor)
# 
# # crop to Europe
# for(i in 1:nlayers(env.bioclim)) {
#   env.bioclim[[i]] <- raster::mask(env.bioclim[[i]], 
#                               as(extent(extent_Europe), 'SpatialPolygons'))}
# 
# 
# save also environmental data as raster file, if wanted
# if (checkSave_precitor == TRUE){
#   raster::writeRaster(env.bioclim, file=paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))
# }

# #- - - - - - - - - - - - - - - - - - - - - 
# ## More OLD CODE ####
# #- - - - - - - - - - - - - - - - - - - - - 
# Agriculture <- raster::raster("D:/00_datasets/LandCover/V008_Agriculture/Agriculture_2012_noGrid_WGS84.tif")
# plot(Agriculture)
# 
# r <- raster::raster("D:/00_datasets/00_Trash/baseline_2km_final.tif")
# r <- rgdal::readOGR("D:/00_datasets/00_Trash/baseline.shp")
# 
# extent(Agriculture) <- extent(r)
# 
# r2 <- rasterToPolygons(r)
# 
# raster::zonal(Agriculture, r, fun="mean")
# 
# raster::rasterize(r, Agriculture, fun="mean")
# 
# 
# test <- raster::extract(Agriculture, r, mean)

# #- - - - - - - - - - - - - - - - - - - - - 
# forestVarsPath <- list.files("./DATA/RASTER/ForestVars", pattern=".tif$", full.names = TRUE)
# rsVarsPath <- list.files("./DATA/RASTER/RemoteSensing", pattern=".tif$", full.names = TRUE)
# soilVarsPath <- list.files("./DATA/RASTER/Soil", pattern=".tif$", full.names = TRUE)
# topoVarsPath <- list.files("./DATA/RASTER/Topo", pattern=".tif$", full.names = TRUE)
# 
# allVars <- c(forestVarsPath, rsVarsPath, soilVarsPath, topoVarsPath)
# varNames <- c(bioVarNames, gsub(".tif","",basename(allVars)))

# rstStack <- raster::stack(raster::stack(climVarsPath)[[selVarsClim]],raster::stack(allVars))
# names(rstStack) <- varNames
# rstData <- na.omit(values(rstStack))
#
#saveRDS(rstData, file = "./DATA/fullRasterStack.rds")
