#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare input data for SDMs          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(biomod2) # also to create pseudo-absences

library(parallel)
library(doParallel)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0("./results/Species_list_", Taxon_name, ".csv"))
speciesSub <- speciesNames %>% filter(NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique() %>% c()
#speciesSub <- speciesNames %>% filter(family == "Lumbricidae" & NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique()
speciesSub <- c(speciesSub$SpeciesID)

#- - - - - - - - - - - - - - - - - - - - -
# note: we will load the datasets before each individual model

# load environmental variables
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km_normalized.grd"))

load(paste0(here::here(), "/results/EnvPredictor_2km_df_normalized.RData")) #Env_norm_df

# Calculate the number of cores
no.cores <-  parallel::detectCores()/2 

#- - - - - - - - - - - - - - - - - - - - -
## Prepare data ####
mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))

registerDoParallel(no.cores)
number_records <- foreach(spID = speciesSub, 
        .combine = rbind,
        .export = c("Env_norm", "mySpeciesOcc"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({
          
          myResp <- mySpeciesOcc[!is.na(mySpeciesOcc[,spID]), c("x","y",spID)]
          
          myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = as.numeric(myResp[,spID]),
                                                        expl.var = Env_norm,
                                                        resp.xy = myResp[,c("x", "y")],
                                                        resp.name = spID,
                                                        PA.nb.rep = 1,
                                                        PA.nb.absences = 10000,
                                                        PA.strategy = "random")
          
          # save data
          save(myBiomodData, file=paste0(here::here(), "/intermediates/BIOMOD_data/BiomodData_", Taxon_name,"_", spID, ".RData"))
          
          #- - - - - - - - - - - - - - - - - - - - -
          ## prepare for Identify top predictors (Maxent)
          # subset covarsNames
          myBiomodData@data.env.var <- myBiomodData@data.env.var[,colnames(myBiomodData@data.env.var) %in% covarsNames]
          
          myData <- cbind(myBiomodData@data.species, myBiomodData@coord, myBiomodData@data.env.var)
          myData$SpeciesID <- spID
          myData <- myData %>% rename("occ" = "myBiomodData@data.species")
          myData[is.na(myData$occ),"occ"] <- 0
          
          print(myData[,c("x", "y","occ", "SpeciesID")])
          
          random.rows <- sample(1:nrow(myData), nrow(myData))
          
          validation <- myData[random.rows[1:round(0.2*nrow(myData))], 
                               c("x","y", "SpeciesID", "occ", covarsNames[covarsNames %in% colnames(myData)])]
          
          training <- myData[random.rows[round(0.2*nrow(myData)):nrow(myData)],]
          
          # subset uncorrelated covariates
          training <- training[, c("occ", covarsNames[covarsNames %in% colnames(myData)])]
          
          # save all datasets
          save(training, file=paste0(here::here(), "/intermediates/MaxEnt_data/MaxentData_train_", Taxon_name,"_", spID, ".RData"))
          save(validation, file=paste0(here::here(), "/intermediates/MaxEnt_data/MaxentData_valid_", Taxon_name,"_", spID, ".RData"))
          
          rm(myBiomodData, myResp, myRespCoord, spID, na.id, myData, training, validation)
        })}

stopImplicitCluster()

number_records


