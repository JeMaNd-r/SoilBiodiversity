#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Point occurrences to grid            #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# Note: to run the other grid, you have to comment & un-comment 4 lines

## load grid
r <- raster::raster("D:/00_datasets/Grids/grid_2k_0p016.tif")

## Load occurrence data
occ <- read.csv(file=paste0(here::here(), "/results/Occurrences_", Taxon_name, ".csv"))
occ$occ <- 1
occ_stack <- raster()
occ_points <- data.frame(x=0, y=0)[0,]

speciesNames$Records <- NA

#- - - - - - - - - - - - - - - - - - - - - - - 
## For loop through all species ####
for(sp in unique(speciesNames$SpeciesID)){
  
  # ignore errors
  try({
  
  # load occurrences for specific species
  temp_occ <- occ[occ$SpeciesID==sp & !is.na(occ$SpeciesID), c("x", "y", "occ", "SpeciesID")]
  
  ## save number of records for later ####
  speciesNames[speciesNames$SpeciesID==sp,]$Records <- nrow(temp_occ)
    
  # make occurrences as SpatialPoints object
  coordinates(temp_occ) <- ~x+y
  
  # make points to raster
  occ_grid <- raster::rasterize(temp_occ, r, "occ", fun=sum)
  names(occ_grid) <- sp

  # crop to raster extent
  occ_grid <- raster::mask(occ_grid, r)
  
  # # add to raster stack
  # occ_stack <- raster::stack(occ_stack, occ_grid)
  
  #plot(occ_grid)
  
  # add to point data frame
  temp_points <- as.data.frame(rasterToPoints(occ_grid))
  occ_points <- dplyr::full_join(occ_points, temp_points)
  
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
## Save ####
# # save raster stack
# raster::writeRaster(occ_stack,file=paste0(here::here(), "/results/OccurrenceGrid_", Taxon_name, ".grd"), format="raster")

# save point data frame
write.csv(occ_points, file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"), row.names=F)

occ_points %>% pivot_longer(cols=(colnames(occ_points %>% dplyr::select(-x, -y)))) %>% filter(!is.na(value))
# 20,500 records in total

# # save individual species as own files
# raster::writeRaster(occ_stack, filename=names(occ_stack), bylayer=TRUE, format="GTiff")

write.csv(speciesNames, file=paste0(here::here(), "/results/Species_list_", Taxon_name, ".csv"), row.names = F)



