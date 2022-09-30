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
covarsNames <- c("MAT", "Dist_Coast", "MAP_Seas", "CEC", "Elev",
                 "P", "Pop_Dens", "Agriculture", "pH", "Clay.Silt")

#- - - - - - - - - - - - - - - - - - - - -
# note: we will load the datasets before each individual model

# load environmental variables (for projections)
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km_normalized.grd"))
#Env_norm <- stack(Env_norm)

# Calculate the number of cores
no.cores <-  parallel::detectCores()/2 

#- - - - - - - - - - - - - - - - - - - - -
## Prepare data ####
mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))

registerDoParallel(no.cores)
foreach(spID = speciesSub, 
        .export = c("Env_norm", "form", "mySpeciesOcc"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({
          
          myResp <- as.numeric(mySpeciesOcc[,spID])
          
          # get NAs id
          na.id <- which(is.na(myResp))
          # remove NAs to enforce PA sampling to be done on explanatory rasters
          myResp <- myResp[-na.id]
          
          myRespCoord <- mySpeciesOcc[-na.id,c('x','y')]
          
          myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                                        expl.var = Env_norm,
                                                        resp.xy = myRespCoord,
                                                        resp.name = spID,
                                                        PA.nb.rep = 1,
                                                        PA.nb.absences = 10000,
                                                        PA.strategy = "random")
          
          # save data
          save(myBiomodData, file=paste0(here::here(), "/results/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData"))
          
          rm(myBiomodData, myResp, myRespCoord, spID, na.id)
        })}

stopImplicitCluster()
