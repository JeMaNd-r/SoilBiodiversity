#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare predictor variables          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load the 1km grid
grid1k <- raster::raster("D:/00_datasets/Grids/grid_1k_0p008.tif")

# load the predictor table containing the individual file names
pred_tab <- readr::read_csv(file=paste0("doc/Env_Predictors_table.csv"))

# list the names of the variable's folders that will be included in the analysis
folders <- c(list.dirs("D:/00_datasets/Climate", recursive=F))
folders <- c(folders, list.dirs("D:/00_datasets/LandCover", recursive=F)) 
folders <- c(folders, list.dirs("D:/00_datasets/Location", recursive=F))
folders <- c(folders, list.dirs("D:/00_datasets/Soil", recursive=F))

folders <- folders[stringr::str_detect(folders,"V[:digit:]{3}_")]

#- - - - - - - - - - - - - - - - - - - - - 
## First, we need to use zonal statistics to project all variables to the 1km
# grid. This didn't work in ArcGIS because of too many unique ID values.

# prepare for parallal processing
cluster <- makeCluster(8)
doParallel::registerDoParallel(cluster)

# loop through the variable's folders
foreach(i = 1:length(folders), .packages = c("raster")) %doPar% {
  
  # define variable name (= folder name)
  temp_pred <- stringr::str_extract(folders[i], "V[:digit:]{3}_[:alpha:]*")
  
  # extract file name
  temp_file <- as.character(pred_tab[pred_tab$ID==stringr::str_extract(temp_pred, "V[:digit:]{3}") &
                          !is.na(pred_tab$ID), "File_name_processed"])
  
  # read tif files with raw predictor variables across Europe
  temp_raster <- raster::raster(paste0(folders[i], "/", temp_file))
  
  # crop to same extent
  temp_raster <- crop(temp_raster, extent(grid1k))
  extent(temp_raster) <- extent(grid1k) 
  
  # calculate average per grid cell
   <- raster::zonal(x=temp_raster, z=grid1k, fun='mean', digits=3, na.rm=TRUE) 

  # create raster based on these averaged values
  temp_1k <- setValues(grid1k)
  
  # extract short predictor's name
  temp_name <- pred_tab[pred_tab$ID==stringr::str_extract(temp_pred, "V[:digit:]{3}") &
             !is.na(pred_tab$ID), "Predictor"]
  
  # save raster
  raster::writeRaster(temp_raster, file=paste0(folders[i], "/", temp_name, "_1km_mean", ".tif"))
}


#- - - - - - - - - - - - - - - - - - - - - 
## Merge all raster files ####
#...


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
