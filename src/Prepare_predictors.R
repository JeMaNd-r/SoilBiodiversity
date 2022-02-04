#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare predictor variables          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#read.csv()

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
# 
test <- raster::extract(Agriculture, r, mean)


#- - - - - - - - - - - - - - - - - - - - - 
# get Worldclim raster data layers
env.bioclim <- raster::getData('worldclim', var='bio', res=10)

#summary(env.bioclim[[1]])
plot(env.bioclim[[1]], main="BioClim 1 (raw)")

## Correlations
env.cor <- sdmpredictors::pearson_correlation_matrix(env.bioclim)

# plot correlations
corrplot::corrplot.mixed(env.cor)

pdf(file=paste0("Predictor_correlation_", Taxon_name, ".pdf"))
sdmpredictors::plot_correlation(env.cor)
dev.off()

sdmpredictors::correlation_groups(env.cor)

# crop to Europe
for(i in 1:nlayers(env.bioclim)) {
  env.bioclim[[i]] <- raster::mask(env.bioclim[[i]], 
                              as(extent(extent_Europe), 'SpatialPolygons'))}


# save also environmental data as raster file, if wanted
if (checkSave_precitor == TRUE){
  raster::writeRaster(env.bioclim, file=paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))
}


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
