#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare protected area raster        #
#          author: Romy Zeiss               #
#            date: 2022-05-04               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(raster)
library(rgdal)
library(sf)
library(exactextractr)

# change temporary directory for files
raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

# load 5km grid
grid5k <- raster::raster("D:/00_datasets/Grids/grid_5k_0p016.tif")
raster::crs(grid5k) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# # load protected area network
# protect_raw_0 <- raster::shapefile("D:/_students/Romy/SoilBiodiversity/data/Shapefiles/WDPA_WDOECM_Dec2021_Public_EU_shp/WDPA_WDOECM_Dec2021_Public_EU_shp_0/WDPA_WDOECM_Dec2021_Public_EU_shp-polygons.shp")
# protect_raw_1 <- raster::shapefile("D:/_students/Romy/SoilBiodiversity/data/Shapefiles/WDPA_WDOECM_Dec2021_Public_EU_shp/WDPA_WDOECM_Dec2021_Public_EU_shp_1/WDPA_WDOECM_Dec2021_Public_EU_shp-polygons.shp")
# protect_raw_2 <- raster::shapefile("D:/_students/Romy/SoilBiodiversity/data/Shapefiles/WDPA_WDOECM_Dec2021_Public_EU_shp/WDPA_WDOECM_Dec2021_Public_EU_shp_2/WDPA_WDOECM_Dec2021_Public_EU_shp-polygons.shp")
# 
# # combine all three shapefiles into 1
# protect_raw <- raster::bind(protect_raw_0, protect_raw_1)
# protect_raw <- raster::bind(protect_raw, protect_raw_2)
# 
# # save combined shapefile
# #rgdal::writeOGR(protect_raw, dsn=getwd(), "WDPA_WDOECM_Dec2021_EU_merged", driver="ESRI Shapefile", overwrite_layer=TRUE)
# 
# rm(protect_raw_0, protect_raw_1, protect_raw_2)

protect_raw <- rgdal::readOGR(paste0(getwd(),"/WDPA_WDOECM_Dec2021_EU_merged.shp"))

# create empty stack to store PA types as layers
protect_stack <- raster::stack()

# split multipolygon by type of protection area
# type of protected area: protect_raw$DESIG_ENG
# IUCN category: IUCN_CAT
for(i in unique(protect_raw$IUCN_CAT)){
	temp_shp <- protect_raw[protect_raw$IUCN_CAT == i,]
	#rgdal::writeOGR(temp_shp, dsn=getwd(), i, driver="ESRI Shapefile", overwrite_layer=TRUE)
	
	# calculate percent of grid cells covered by PA type i
	#temp_cover <- raster::rasterize(temp_shp, grid5k, getCover = TRUE)
	temp_cover <- exactextractr::coverage_fraction(grid5k, sf::st_combine(as(temp_shp, "sf")))[[1]]
	names(temp_cover) <- i

	# replace values higher than 1 by 1
	temp_cover[temp_cover > 1] <- 1
	temp_cover[temp_cover < 0] <- 0

	protect_stack <- raster::stack(protect_stack, temp_cover)

	print(paste0(i, " ready."))
}

protect_stack

# crop to grid extent
protect_stack_crop <- raster::mask(protect_stack, grid5k)
protect_stack_crop <- raster::stack(protect_stack_crop)

# save raster stack
raster::writeRaster(protect_stack_crop, filename="WDPA_WDOECM_IUCNcat.grd")

# transform protected cover into dataframe
protect_df <- as.data.frame(raster::rasterToPoints(protect_stack))
save(protect_df, file="WDPA_WDOECM_IUCNcat_df.RData")

## Crop to Germany
# load germany borders
german_sf <- sf::st_as_sf(maps::map(region = "germany", plot = FALSE, fill = TRUE))

protect_stack_GE <- raster::mask(protect_stack_crop, sf::as_Spatial(german_sf))
protect_stack_GE <- raster::crop(protect_stack_GE, raster::extent(sf::as_Spatial(german_sf)))

raster::writeRaster(protect_stack_GE, filename="I:/eie/==PERSONAL/RZ_SoilBON/2022_MacroCourse/data/WDPA_WDOECM_IUCNcat_Germany.grd", overwrite=T)

pdf(file="I:/eie/==PERSONAL/RZ_SoilBON/2022_MacroCourse/figs/ProtectedAreas_IUCNcategories_Germany.pdf", height=15, width = 18)
raster::plot(protect_stack_GE, maxnl=70)
dev.off()


