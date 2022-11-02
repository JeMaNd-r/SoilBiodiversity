#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Predict SDMs in future climate       #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(raster)

library(biomod2) # also to create pseudo-absences

library(parallel)
library(doParallel)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0("./results/Species_list_", Taxon_name, ".csv"))
speciesSub <- speciesNames %>% filter(NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique() %>% c()
#speciesSub <- speciesNames %>% filter(family == "Lumbricidae" & NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique()
speciesSub <- c(speciesSub$SpeciesID)

# covariates in order of importance (top 10 important)
covarsNames <- c("MAT", "MAP_Seas", "Dist_Coast", "Agriculture", "pH", 
                 "P", "CEC", "Elev", "Clay.Silt", "Pop_Dens")

# define future scenarios
scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
                               "mri-esm2-0", "ukesm1-0-ll"), "_",
                             rep(c("ssp126", "ssp370", "ssp585"),5)))


# Calculate the number of cores
no.cores <-  parallel::detectCores()/2 


#- - - - - - - - - - - - - - - - -
## Predict in future climate ####

setwd(here::here())
setwd(paste0(here::here(), "/results/", Taxon_name))

# load current environmental variables (for projections)
load(paste0(here::here(),"/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df

no.cores <- 3
registerDoParallel(no.cores)
# foreach(spID = speciesSub,
#         .export = c("Env_norm", "Env_norm_df", "form", "Taxon_name"),
#         .packages = c("tidyverse","biomod2")) %dopar% { try({ 
#   

for(spID in speciesSub){ try({ 
            
            # list files in species-specific BIOMOD folder
            temp_files <- list.files(paste0(here::here(), "/results/biomod_files/", stringr::str_replace(spID, "_", ".")), full.names = TRUE)
            
            temp_files
            
            # load model output
            myBiomodModelOut <- temp_files[stringr::str_detect(temp_files,"Modeling.models.out")]
            print(myBiomodModelOut)
            myBiomodModelOut <-get(load(myBiomodModelOut))
            
            # load ensemble model output
            myBiomodEM <- temp_files[stringr::str_detect(temp_files,"Modelingensemble.models.out")]
            myBiomodEM <- get(load(myBiomodEM))
            
            for(no_future in scenarioNames){
              
              ## load env. data with future climate (MAP, MAT, MAP_Seas)
              #load(paste0(here::here(), "/results/_FutureEnvironment/EnvPredictor_2041-2070_", no_future, "_5km_df_normalized.RData")) #temp_Env_df
              
              # one loop per future climate subset, one with both future, each one with only 1 future and 1 current climate
              for(subclim in c("TP", "T", "P")){
                
                if(subclim=="TP"){
                  temp_Env_sub <- temp_Env_df[,c("x", "y", colnames(temp_Env_df)[colnames(temp_Env_df) %in% covarsNames])]
                  
                  temp_MAT <- raster::raster(paste(here::here(), "/data_environment/", temp_name))
                  temp_raster <- raster::scale(temp_raster)
                  
                  temp_Env_norm[[temp_name]] <- temp_raster
                  
                }
                
                if(subclim=="P"){
                  temp_Env_sub <- temp_Env_df[,c("x", "y", colnames(temp_Env_df)[colnames(temp_Env_df) %in% covarsNames])]
                  temp_Env_sub$MAP_Seas <- Env_norm_df$MAP_Seas
                  temp_Env_sub$MAP <- Env_norm_df$MAP
                }
                
                if(subclim=="T"){
                  temp_Env_sub <- temp_Env_df[,c("x", "y", colnames(temp_Env_df)[colnames(temp_Env_df) %in% covarsNames])]
                  temp_Env_sub$MAT <- Env_norm_df$MAT
                }
                
                
                ## NOTE: because biomod output can hardly be stored in list file, we will do calculations based on model output now
                # project single models (also needed for ensemble model)
                myBiomodProj <- biomod2::BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                           new.env = temp_Env_sub[,colnames(temp_Env_sub) %in% covarsNames],        #column/variable names have to perfectly match with training
                                                           proj.name = "modeling",  #name of the new folder being created
                                                           selected.models = "all", #use all models
                                                           binary.meth = NULL,     #binary transformation according to criteria, or no transformation if NULL
                                                           compress = TRUE,         #compression format of objects stored on hard drive
                                                           build.clamping.mask = TRUE, #TRUE: clamping mask will be saved on hard drive different
                                                           do.stack = TRUE,         #save output projections as rasterstack (if not too heavy)
                                                           output.format = ".RData", #what format should projections have: RData, grd or img
                                                           keep.in.memory = TRUE)  #FALSE: only story link to copy to projection file
                
                # project ensemble of all models
                myBiomodEnProj <- biomod2::BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                                                      EM.output = myBiomodEM,
                                                                      #... same arguments as above could be added but are not necessary when loading myBiomodProj
                                                                      selected.models = "all")
                
                #temp_predict_time <- proc.time()[3] - tmp
                
                # save predictions as raster file
                temp_prediction <- myBiomodEnProj@proj@val[,2]
                temp_prediction <- as.numeric(temp_prediction)
                # add names of grid cell (only for those that have no NA in any layer)
                names(temp_prediction) <- rownames(temp_Env_sub)
                temp_prediction <- as.data.frame(temp_prediction)
                temp_prediction$x <- temp_Env_sub$x
                temp_prediction$y <-temp_Env_sub$y
                temp_prediction <- temp_prediction %>% full_join(temp_Env_sub %>% dplyr::select(x,y))
                colnames(temp_prediction)[1] <- "layer"
                temp_prediction$layer <- temp_prediction$layer / 1000
                
                save(temp_prediction, file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_2041-2070_", no_future, "_", subclim, "_biomod_", spID,  ".RData")) 
                rm(temp_prediction, temp_Env_sub, myBiomodEnProj, myBiomodProj)
              }
            }
            
          })}
          
 stopImplicitCluster()
          