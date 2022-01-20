#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Point occurrences to grid            #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

occ <- read.csv(file=paste0(here::here(), "/results/Occurrences_", Taxon_name, ".csv"))
occ$occ <- 1

# make occurrences as SpatialPoints object
occ <- occ[occ$SpeciesID==spID, c("x", "y", "occ", "SpeciesID")]
coordinates(occ) <- ~x+y


## create grid (or load grid)
r <- raster::raster(xmn=extent_Europe[1], xmx=extent_Europe[2], 
                    ymn= extent_Europe[3], ymx=extent_Europe[4],
                    res=0.5)

# define reference system: LAEA
crs(r) <- CRS("+init=epsg:4238")

# make points to raster
occ_grid <- raster::rasterize(occ, r, "occ", fun=sum)
plot(occ_grid)

plot(occ)
