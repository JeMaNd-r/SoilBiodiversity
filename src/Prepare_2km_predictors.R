#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare predictor variables          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# change temporary directory for files
raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

# load 2km grid
grid2k <- raster::raster("D:/00_datasets/Grids/grid_2k_0p016.tif")
raster::crs(grid2k) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# load the predictor table containing the individual file names
pred_tab <- readr::read_csv(file=paste0(here::here(), "/doc/Env_Predictors_table.csv"))

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
 
makeTo2kmGrid <- function(raster_grid, temp_path, temp_file=NULL, file_name=NULL){
  
  # makeTo2kmGrid takes the temp_file to re-project it into the same extent and resolution as raster_grid.
  
  # raster_grid:  RasterLayer of a grid.
  # temp_path:    Path to the input file and at the same time output directory.
  # temp_file:    Name of the input file that will be re-projected.
  # file_name:    Name of the output file.
  
  # define variable name (= folder name)
  temp_pred <- stringr::str_split_fixed(temp_path,"/",4)[4]
  #temp_pred <- stringr::str_extract(temp_path, "V[:digit:]{3}_[:alpha:]*_*[:alpha:]*")
  
  # extract file name
  if(is.null(temp_file)) {
    temp_file <- as.character(pred_tab[pred_tab$Predictor_long==temp_pred, "File_name_processed"])
    
    # in case that 1km and 2km have differencly processed files:
    if(stringr::str_detect(temp_file, "1km")) {
      temp_file <- stringr::str_replace(temp_file, "1km", "2km")}
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
    file_name <- paste0(temp_name, "_2km_mean.tif")
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
  makeTo2kmGrid(temp_path=folders[i], raster_grid=grid2k)
})}

# Note: Error in some variables will be solved individually

#- - - - - - - - - - - - - - - - - - - - - 
## Calculate some missing 2km grids ####

# Latitude
temp_raster <- as.data.frame(raster::rasterToPoints(grid2k))
temp_raster$Latitude <- temp_raster$y
temp_raster_mean <- raster::rasterFromXYZ(temp_raster %>% dplyr::select(-grid_2k_0p016))
raster::writeRaster(temp_raster_mean, file="D:/00_datasets/Location/V019_Lat/Lat_2km_mean.tif", overwrite=T)

#- - - - - - - - - - - - - - - - - - - - - 
## Mask selected 2km grids ####

# load mask file (one of CORINE land cover files)
temp_mask <- raster::raster("D:/00_datasets/LandCover/V008_Agriculture/Agriculture_2012_2km_mean.tif")

# Dist_Coast
temp_raster <- raster::raster("D:/00_datasets/Location/V016_Dist_Coast/Dist_Coast_2km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/Location/V016_Dist_Coast/Dist_Coast_2km_mean.tif", overwrite=T)

# Dist_River
temp_raster <- raster::raster("D:/00_datasets/Location/V017_Dist_River/Dist_River_2km_mean.tif")
temp_raster <- raster::mask(temp_raster, temp_mask)
raster::writeRaster(temp_raster, file="D:/00_datasets/Location/V017_Dist_River/Dist_River_2km_mean.tif", overwrite=T)

#- - - - - - - - - - - - - - - - - - - - - 
## Merge all raster files, 2km ####
# get names of the 2km files
files <- list.files(folders, include.dirs = F, recursive=F, full.names = T)
stack_files <- files[stringr::str_detect(files, "_2km_mean.tif$")]

# remove Forest (but keep Forest_Deci and Forest_Coni)
stack_files <- stack_files[!stringr::str_detect(stack_files, "Forest_2012_")]

# create empty stacj
Env <- raster::stack()

# load and merge them
for(i in 1:length(stack_files)){
  temp_raster <- raster::raster(stack_files[i])
  
  names(temp_raster) <- gsub("_2k.*", "", names(temp_raster))
  
  #temp_raster <- raster::mask(temp_raster, grid2k)
  print(extent(temp_raster)); print(names(temp_raster))
  #Env <- raster::stack(Env, temp_raster)
}

# save raster
raster::writeRaster(Env, file=paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"), overwrite=T)

## cut to grid (should not be necessary anymore...)
#Env <- raster::mask(Env, grid1k)
# 
# # trim unnecessary margins
# Env_trimmed <- trim(Env, values = NA)
# Env_trimmed[is.na(Env_trimmed)] <- 0
# 
# Env <- raster::stack(Env)

# save raster stack plots
par(mfrow=c(5,6))
#pdf(file=paste0(here::here(), "/figures/Predictors_Europe_1km.pdf"), height=15, width = 18)
raster::plot(Env$Aridity)
raster::plot(Env$MAP)
raster::plot(Env$MAP_Seas)
raster::plot(Env$MAT)
raster::plot(Env$MAT_Seas)
raster::plot(Env$Snow)
raster::plot(Env$Agriculture)
raster::plot(Env$Dist_Urb)
raster::plot(Env$Forest_Coni_2012)
raster::plot(Env$Forest_Deci_2012)

raster::plot(Env$Pastures)
raster::plot(Env$Pop_Dens)
raster::plot(Env$Shrubland)
raster::plot(Env$NDVI)
raster::plot(Env$Aspect)
raster::plot(Env$Dist_Coast)
raster::plot(Env$Dist_River)
raster::plot(Env$Elev)
raster::plot(Env$Lat)
raster::plot(Env$Slope)

raster::plot(Env$CEC)
raster::plot(Env$Clay)
raster::plot(Env$Cu)
raster::plot(Env$Hg)
raster::plot(Env$N)
raster::plot(Env$P)
raster::plot(Env$pH)
raster::plot(Env$Silt)
raster::plot(Env$SOC)
raster::plot(Env$Moisture)
dev.off()


