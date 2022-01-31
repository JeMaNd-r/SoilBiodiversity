#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#        Vizualize output of SDMs           #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# code to map RF downsampled
library(tmap)
library(scales)

# rescale and stack the predicted maps
maps <- list("default" = pred_def,
             "weighted" = pred_wtd,
             "regression" = pred_reg,
             "downsample" = pred_dws,
             "equal_sample" = pred_eql,
             "shallow_tuned" = pred_shall_tune, 
             "shallow" = pred_rng_shallow, 
             "truth" = virtual_sp) %>% 
  lapply(function(x) calc(x, function(y) rescale(y, c(0,1)))) %>% 
  raster::stack()

# plot the maps
tm_shape(maps) +
  tm_raster(title = "Likelihood", 
            palette = rev(terrain.colors(30)), 
            style = "cont",
            breaks = c(0, 0.25, 0.5, 0.75, 1))





## SSDM or SDM code...
# # define species names (see SDM script)
# sp.occ <- read.csv(file=paste0(here::here(), "/results/Occurrences_wide_", Taxon_name, ".csv"))
# sp.occ <- colSums(mySpeciesOcc[,speciesNames$SpeciesID], na.rm=T) >=5
# sp.names <- speciesNames[sp.occ, "SpeciesID"]
# 
# # replace "_" with point to guarantee matching folder names
# sp.names <- stringr::str_replace(sp.names, "_", ".")
# 
# # set working directory
# setwd(paste0(here::here(), "/results/SDMs/"))
# 
# 
# # define a mask of studied
# alphaMap <- reclassify(subset(myExpl,1), c(-Inf,Inf,0))
# 
# # add all other species map
# for(sp.n in sp.names){
#   # add layer
#   alphaMap <-
#     alphaMap +
#     subset(stack(file.path(sp.n,
#                            "proj_current",
#                            paste("proj_current_",
#                                  sp.n,
#                                  "_TotalConsensus_EMbyTSS_TSSbin.grd", sep=""))), 1)
# }
# # summary of created raster
# alphaMap
# 
# plot(alphaMap, main = expression( paste(alpha, "-diversity based on",
#                                         " TotalConsensus_EMbyTSS_TSSbin outputs")))
