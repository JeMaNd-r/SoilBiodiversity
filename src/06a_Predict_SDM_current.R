#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Predict SDMs in current climate      #
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

# Calculate the number of cores
no.cores <-  parallel::detectCores()/2 

#- - - - - - - - - - - - - - - - -
## Predict in current climate at 5km ####

# load environmental variables (for projections)
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_5km_normalized.grd"))
#Env_norm <- stack(Env_norm)

# as dataframe
load(paste0(here::here(),"/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df

registerDoParallel(3)
foreach(spID = speciesSub,
        .export = c("Env_norm", "Env_norm_df", "form"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({   
          
          # list files in species-specific BIOMOD folder
          temp_files <- list.files(paste0(here::here(), "/results/biomod_files/", stringr::str_replace(spID, "_", ".")), full.names = TRUE)
          
          setwd(paste0(here::here(), "/results/", Taxon_name))
          
          # load model output
          myBiomodModelOut <- temp_files[stringr::str_detect(temp_files,"Modeling.models.out")]
          print(myBiomodModelOut)
          myBiomodModelOut <-get(load(myBiomodModelOut))
          
          # load ensemble model output
          myBiomodEM <- temp_files[stringr::str_detect(temp_files,"Modelingensemble.models.out")]
          myBiomodEM <- get(load(myBiomodEM))
          
          tmp <- proc.time()[3]
          ## NOTE: because biomod output can hardly be stored in list file, we will do calculations based on model output now
          # project single models (also needed for ensemble model)
          myBiomodProj <- biomod2::BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                     new.env = Env_norm_df[,colnames(Env_norm_df) %in% covarsNames],        #column/variable names have to perfectly match with training
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
          
          temp_predict_time <- proc.time()[3] - tmp
          
          # extracting the values for ensemble validation
          myEnProjDF <- as.data.frame(get_predictions(myBiomodEM)[,2]) #for weighted probability mean
          
          # see the first few validations
          # note: the prediction scale of biomod is between 0 and 1000
          #head(myEnProjDF)
          
          # Get model evaluation values for later
          myBiomodModelEval <- as.data.frame(biomod2::get_evaluations(myBiomodEM))
          
          # Calculate variable importance across all PA sets, eval runs and algorithms
          # and extract only the one for weighed mean predictions (for later)
          temp_varImp <- biomod2::get_variables_importance(myBiomodEM)[, , 2]
          # average across 10 runs
          temp_varImp <- temp_varImp %>% as.data.frame() %>%
            mutate(mean_vi = as.numeric(rowMeans(temp_varImp, na.rm=T)),
                   Predictor = rownames(temp_varImp))
          colnames(temp_varImp)[colnames(temp_varImp) == "mean_vi"] <- "biomod"
          temp_varImp <- temp_varImp[,c("biomod", "Predictor")]
          
          # save predictions as raster file
          temp_prediction <- myBiomodEnProj@proj@val[,2]
          temp_prediction <- as.numeric(temp_prediction)
          # add names of grid cell (only for those that have no NA in any layer)
          names(temp_prediction) <- rownames(Env_norm_df)
          temp_prediction <- as.data.frame(temp_prediction)
          temp_prediction$x <- Env_norm_df$x
          temp_prediction$y <- Env_norm_df$y
          temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
            rename("layer" = temp_prediction)
          temp_prediction$layer <- temp_prediction$layer / 1000
          
          temp_runs <- 1
          
          biomod_list <- list(time_predict=temp_predict_time, validation=myBiomodModelEval, prediction=temp_prediction, varImp=temp_varImp)
          save(biomod_list, file=paste0("./_SDMs/SDM_biomod_", spID, ".RData"))
          
          rm(biomod_list, temp_predict_time, temp_runs, temp_prediction, temp_varImp, myBiomodEnProj, myBiomodProj, myBiomodModelEval, myEnProjDF, myBiomodModelOut, myBiomodEM)
          
          setwd(here::here())
          
        })}

stopImplicitCluster()
