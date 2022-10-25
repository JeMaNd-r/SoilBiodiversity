#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Point occurrences to grid            #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)
library(terra)

library(raster)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0("./data_raw/Species_list_", Taxon_name, ".csv"))

# Note: to run the other grid (5km, 10km), you have to comment & un-comment 4 lines

## load grid
r <- terra::rast(paste0(here::here(),"/data_raw/grid_2k_0p016.tif"))

## Load occurrence data
occ <- read.csv(file=paste0(here::here(), "/results/Occurrences_", Taxon_name, ".csv"))
occ$occ <- 1
occ_points <- data.frame(x=0, y=0, year=1800)[0,]

speciesNames$Records <- NA
speciesNames$Records_unique <- NA

#- - - - - - - - - - - - - - - - - - - - - - - 
## For loop through all species ####
for(spID in unique(speciesNames$SpeciesID)){
  
  # ignore errors
  try({
  
  # load occurrences for specific species
  temp_occ <- occ[occ$SpeciesID==spID & !is.na(occ$SpeciesID), c("x", "y", "occ", "SpeciesID", "year")]
  
  ## save number of records for later ####
  speciesNames[speciesNames$SpeciesID==spID,]$Records <- nrow(temp_occ)
  speciesNames[speciesNames$SpeciesID==spID,]$Records_unique <- nrow(temp_occ %>% dplyr::select(-year) %>% unique())
     
  # make occurrences as SpatVector object
  temp_occ <- terra::vect(temp_occ, crs="+proj=longlat +datum=WGS84", geom=c("x", "y"))
  
  # make points to raster
  occ_grid <- terra::rasterize(temp_occ, r, "occ", fun=sum)
  occ_year <- terra::rasterize(temp_occ, r, "year", fun=max)
  names(occ_grid) <- spID
  names(occ_year) <- "year"

  # crop to raster extent
  occ_grid <- terra::mask(occ_grid, r)
  occ_year <- terra::mask(occ_year, r)
  
  # add to point data frame
  temp_points <- as.data.frame(occ_grid, xy=TRUE, row.names = FALSE)
  occ_year <- as.data.frame(occ_year, xy=TRUE, row.names = FALSE)

  temp_points <- dplyr::full_join(temp_points, occ_year)

  occ_points <- dplyr::full_join(occ_points, temp_points, by=c("x", "y", "year"))
  
  print(paste0(spID, " was now tried to add to raster stack."))
  print(paste0("It has ", nrow(occ[occ$SpeciesID==spID,] %>% dplyr::select(latitude, longitude) %>% unique()), " records."))
  print("######################################################")
  
  }) # end of try loop
}

## Calculate number of grid cells with presence ####
speciesNames$NumCells_2km <- 0

for(spID in unique(speciesNames$SpeciesID)){ try({
  temp_records <- nrow(occ_points) - as.numeric(summary(occ_points[,spID])["NA's"])
  speciesNames[speciesNames$SpeciesID==spID,"NumCells_2km"] <- temp_records
  temp_records <- 0
}, silent=T)}
#print("Error messages are fine. They relate to species without records.")

#- - - - - - - - - - - - - - - - - - - - - - - 
## Save ####

# save point data frame
write.csv(occ_points, file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"), row.names=F)

occ_points %>% pivot_longer(cols=(colnames(occ_points %>% dplyr::select(-x, -y, -year)))) %>% filter(!is.na(value))
# 26,351 (1990: 19,993 [OLD]) records in total

# # save individual species as own files
# raster::writeRaster(occ_stack, filename=names(occ_stack), bylayer=TRUE, format="GTiff")

write.csv(speciesNames, file=paste0(here::here(), "/results/Species_list_", Taxon_name, ".csv"), row.names = F)



