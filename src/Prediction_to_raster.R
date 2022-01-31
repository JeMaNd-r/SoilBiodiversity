#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Predictions as raster             #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load SDM model output
load(file=paste0(here::here(), "/results/", Taxon_name, "/Predicted_SDMs_", spID, ".RData"))

# load environmental variables as raster stack
Env <- stack(paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))


# First, we need to predict all models to the full spatial extent.
# Afterwards, we can transform the point predictions into a raster file.

for(i in 1:length(names(SDMs))){
  temp.model <- names(SDMs)[[i]]
  
  # # if necessary, unlist models
  # if(length(SDMs[[i]])!=3) {
  #   
  # }
  
  # define background dataset (for testing data)
  modelName <- SDMs[[i]][[2]]
  
  # identify and load all relevant data files
  temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                           pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
  lapply(temp.files, load, .GlobalEnv)
  
  # re-loading the species data (we need the x & y column)
  # load background data (pseudo-absences) for each modeling approach
  load(file=paste0(here::here(), "/results/", Taxon_name, "/PA_Env_", Taxon_name, "_", spID, ".RData"))
  data <- bg.list[[modelName]] %>% rename("occ"=1)
  
  
  
  
  # transform SpatialPoints into raster
  rasterize(x, y, field, filename=paste0(here::here(), ""), na.rm=TRUE, ...)
}