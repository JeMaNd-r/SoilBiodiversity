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

# covariates
corMatPearson <- as.matrix(read.csv(file=paste0(here::here(), "/results/corMatPearson_predictors.csv")))
dimnames(corMatPearson)[[1]] <- dimnames(corMatPearson)[[2]]
# based on Valavi et al. 2021: Pearson 0.8
env_exclude <- caret::findCorrelation(corMatPearson, cutoff = 0.8, names=TRUE)
covarsNames <- dimnames(corMatPearson)[[1]][!(dimnames(corMatPearson)[[1]] %in% env_exclude)]
covarsNames <- covarsNames[covarsNames != "x" & covarsNames != "y"]
# exclude based on VIF
env_vif <- read.csv(file=paste0(here::here(), "/results/VIF_predictors.csv"))
env_exclude <- env_vif %>% filter(is.na(VIF)) %>% dplyr::select(Variables) %>% as.character()
covarsNames <- covarsNames[!(covarsNames %in% env_exclude)]
# excluded:
print("=== We excluded the following variables based on VIF and Pearson correlation: ===")
setdiff(env_vif$Variables, covarsNames)

# final predictor variables
print("=== And we kept the following, final predictor variables: ===")
covarsNames

#- - - - - - - - - - - - - - - - - - - - -
# note: we will load the datasets before each individual model

# load environmental variables
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km_normalized.grd"))

# Calculate the number of cores
no.cores <-  parallel::detectCores()/2 

#- - - - - - - - - - - - - - - - - - - - -
## Prepare data ####
mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))

registerDoParallel(3)
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
          
          rm(myBiomodData, myResp, spID, myData)
        })}

stopImplicitCluster()

number_records


