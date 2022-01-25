#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Point occurrences to grid            #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

## create grid (or load grid)
r <- raster::raster(xmn=extent_Europe[1], xmx=extent_Europe[2], 
                    ymn= extent_Europe[3], ymx=extent_Europe[4],
                    res=0.01, 
                    crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

## Load occurrence data
occ <- read.csv(file=paste0(here::here(), "/results/Occurrences_", Taxon_name, ".csv"))
occ$occ <- 1
occ.stack <- raster()
occ.points <- data.frame(x=0, y=0)[0,]

#- - - - - - - - - - - - - - - - - - - - - - - 
## For loop through all Federal State folders ####
for(sp in speciesNames$SpeciesID){
  
  # ignore errors
  try({
          
  # make occurrences as SpatialPoints object
  temp.occ <- occ[occ$SpeciesID==sp, c("x", "y", "occ", "SpeciesID")]
  coordinates(temp.occ) <- ~x+y
  
  # make points to raster
  occ.grid <- raster::rasterize(temp.occ, r, "occ", fun=sum)
  names(occ.grid) <- sp

  # # add to raster stack
  # occ.stack <- raster::stack(occ.stack, occ.grid)
  
  #plot(occ.grid)
  
  # add to point data frame
  temp.points <- as.data.frame(rasterToPoints(occ.grid))
  occ.points <- dplyr::full_join(occ.points, temp.points)
  
  print(paste0(sp, " was now tried to add to raster stack."))
  print(paste0("It has ", nrow(occ[occ$SpeciesID==sp, c("x", "y", "occ", "SpeciesID")]), "records."))
  print("######################################################")
  
  }) # end of try loop
}

# #- - - - - - - - - - - - - - - - - - - - - - - 
# ## Transform raster stack to occurrence table ####
# occ.points <- rasterToPoints(occ.stack)
# occ.points <- as.data.frame(occ.points)

# ## transform in a long data frame
# occ.points.long <- occ.points %>% pivot_longer(cols=3:ncol(occ.points), names_to = "SpeciesID", values_to = "occ") %>% filter(!is.na(occ))
# 
# plot(occ.stack)
# 
# ggplot(occ.points.long, aes(x=x, y=y, col=SpeciesID))+
#   geom_point()
  

#- - - - - - - - - - - - - - - - - - - - - - - 
## Save ####
# # save raster stack
# raster::writeRaster(occ.stack,file=paste0(here::here(), "/results/OccurrenceGrid_", Taxon_name, ".grd"), format="raster")

# save point data frame
write.csv(occ.points, file=paste0(here::here(), "/results/Occurrence_rasterized_", Taxon_name, ".csv"))

# # save individual species as own files
# raster::writeRaster(occ.stack, filename=names(occ.stack), bylayer=TRUE, format="GTiff")





