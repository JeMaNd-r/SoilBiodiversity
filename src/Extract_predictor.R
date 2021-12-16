#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Extract predictor variables          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# get Worldclim raster data layers
env.bioclim <- getData('worldclim', var='bio', res=10)

#summary(env.bioclim[[1]])
plot(env.bioclim[[1]], main="BioClim 1 (raw)")

#--------------------------------------------------------------#
## Extract climate data at point locations

# load species occurrences
occ <- read.csv(file=paste0(here::here(), "/results/Occurrences_", Taxon_name, ".csv"))

# transform occurrences to raster
coordinates(occ) <- ~ decimalLongitude + decimalLatitude

# extract environmental values
envbd.bioclim <- raster::extract(env.bioclim, occ, df = TRUE)

# #--------------------------------------------------------------#
# ## Outlier detection using bio1(MAT) & bio12(MAP)
# plot(envbd.bioclim$bio1~envbd.bioclim$bio12) # outliers?
# 
# # # histogram first
# # hist(envbd.bioclim$bio1, breaks=20)
# # hist(envbd.bioclim$bio12) #should be normal!
# # hist(log(envbd.bioclim$bio12)) # logarithmic -> better?
# 
# OUTbio12 <- getOutliers(envbd.bioclim$bio12[complete.cases(envbd.bioclim$bio12)],
#                         method="I",
#                         distribution="normal")
# outlierPlot(envbd.bioclim$bio12[complete.cases(envbd.bioclim$bio12)],OUTbio12,mode="qq")
# 
# # lower: OUT$iLeft, upper: OUT$iRight
# # checken, if & which rownumbers are outliers:
# 
# OUTbio12$iLeft
# OUTbio12$iRight
# 
# # combine these vectors
# OUTLIERSbio12 <- c(OUTbio12$iLeft,OUTbio12$iRight)
# 
# # remove that outliers:
# if (length(OUTLIERSbio12) != 0) {         # IF outliers present
#   CHORO.CLIM <- CHOROCLIM[-c(OUTLIERSbio12),] # remove selected
# } else {                                  # else:
#   CHORO.CLIM <- (CHOROCLIM)               # just overwrite
# }

## Save extracted environmental variables

# add coordinates to environment dataframe
envbd.bioclim$Longitude <- coordinates(occ)[,"decimalLongitude"]
envbd.bioclim$Latitude <- coordinates(occ)[,"decimalLatitude"]

# save as csv file
write.csv(envbd.bioclim, file=paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".csv"), row.names = F)

# save also environmental data as raster file
writeRaster(env.bioclim, file=paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))
