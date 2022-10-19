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

library(raster)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0("./results/Species_list_", Taxon_name, ".csv"))

# Note: to run the other grid, you have to comment & un-comment 4 lines

## load grid
r <- raster::raster(paste0(here::here(),"/data/grid_2k_0p016.tif"))

## Load occurrence data
occ <- read.csv(file=paste0(here::here(), "/results/Occurrences_", Taxon_name, ".csv"))
occ$occ <- 1
occ_stack <- raster()
occ_points <- data.frame(x=0, y=0, year=1800)[0,]

speciesNames$Records <- NA

#- - - - - - - - - - - - - - - - - - - - - - - 
## For loop through all species ####
for(sp in unique(speciesNames$SpeciesID)){
  
  # ignore errors
  try({
  
  # load occurrences for specific species
  temp_occ <- occ[occ$SpeciesID==sp & !is.na(occ$SpeciesID), c("x", "y", "occ", "SpeciesID", "year")]
  
  ## save number of records for later ####
  speciesNames[speciesNames$SpeciesID==sp,]$Records <- nrow(temp_occ)
    
  # make occurrences as SpatialPoints object
  coordinates(temp_occ) <- ~x+y
  
  # make points to raster
  occ_grid <- raster::rasterize(temp_occ, r, "occ", fun=sum)
  occ_year <- raster::rasterize(temp_occ, r, "year", fun=max)
  names(occ_grid) <- sp
  names(occ_year) <- "year"

  # crop to raster extent
  occ_grid <- raster::mask(occ_grid, r)
  occ_year <- raster::mask(occ_year, r)
  
  # # add to raster stack
  # occ_stack <- raster::stack(occ_stack, occ_grid)
  
  #plot(occ_grid)
  
  # add to point data frame
  temp_points <- as.data.frame(rasterToPoints(occ_grid))
  occ_year <- as.data.frame(rasterToPoints(occ_year))

  temp_points <- dplyr::full_join(temp_points, occ_year)

  occ_points <- dplyr::full_join(occ_points, temp_points, by=c("x", "y", "year"))
  
  print(paste0(sp, " was now tried to add to raster stack."))
  print(paste0("It has ", nrow(occ[occ$SpeciesID==sp, c("x", "y", "occ", "SpeciesID")]), " records."))
  print("######################################################")
  
  }) # end of try loop
}

# #- - - - - - - - - - - - - - - - - - - - - - - 
# ## Transform raster stack to occurrence table ####
# occ_points <- rasterToPoints(occ_stack)
# occ_points <- as.data.frame(occ_points)

# ## transform in a long data frame
# occ_points_long <- occ_points %>% pivot_longer(cols=3:ncol(occ_points), names_to = "SpeciesID", values_to = "occ") %>% filter(!is.na(occ))
# 
# plot(occ_stack)
# 
# ggplot(occ_points_long, aes(x=x, y=y, col=SpeciesID))+
#   geom_point()
  

## Calculate number of grid cells with presence ####
speciesNames$NumCells_2km <- 0

for(sp in unique(speciesNames$SpeciesID)){ try({
  temp_records <- nrow(occ_points) - as.numeric(summary(occ_points[,sp])["NA's"])
  speciesNames[speciesNames$SpeciesID==sp,"NumCells_2km"] <- temp_records
  temp_records <- 0
}, silent=T)}
print("Error messages are fine. They relate to species without records.")

#- - - - - - - - - - - - - - - - - - - - - - -
# Fix species names Aporr_cali and Aporr_rose
occ_points$Aporr_rose <- rowSums(cbind(occ_points$Aporr_rose, occ_points$Allol_rose), na.rm=T) 
occ_points$Aporr_cali <- rowSums(cbind(occ_points$Aporr_cali, occ_points$Nicod_cali), na.rm=T)
occ_points <- occ_points %>% dplyr::select(-Allol_rose, -Nicod_cali)
occ_points[occ_points==0] <- NA

head(occ_points)

#- - - - - - - - - - - - - - - - - - - - - - - 
## Save ####
# # save raster stack
# raster::writeRaster(occ_stack,file=paste0(here::here(), "/results/OccurrenceGrid_", Taxon_name, ".grd"), format="raster")

# save point data frame
write.csv(occ_points, file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"), row.names=F)

occ_points %>% pivot_longer(cols=(colnames(occ_points %>% dplyr::select(-x, -y, -year)))) %>% filter(!is.na(value))
# 27,086 (1990: 19,993) records in total

# # save individual species as own files
# raster::writeRaster(occ_stack, filename=names(occ_stack), bylayer=TRUE, format="GTiff")

write.csv(speciesNames, file=paste0(here::here(), "/results/Species_list_", Taxon_name, ".csv"), row.names = F)



