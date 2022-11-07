#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Extract model performance            #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(raster)

library(biomod2) # also to create pseudo-absences

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0("./results/Species_list_", Taxon_name, ".csv"))
speciesSub <- speciesNames %>% filter(NumCells_2km_biomod >=100) %>% dplyr::select(SpeciesID) %>% unique() %>% c()
#speciesSub <- speciesNames %>% filter(family == "Lumbricidae" & NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique()
speciesSub <- c(speciesSub$SpeciesID)

# covariates in order of importance (top 10 important)
covarsNames <- c("MAT", "MAP_Seas", "Dist_Coast", "Agriculture", "pH", 
                 "P", "CEC", "Elev", "Clay.Silt", "Pop_Dens")

# define future scenarios
scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
                               "mri-esm2-0", "ukesm1-0-ll"), "_",
                             rep(c("ssp126", "ssp370", "ssp585"),5)))


#- - - - - - - - - - - - - - - - - - - - -
## Model performance ####
#- - - - - - - - - - - - - - - - - - - - -

data_eval <- data.frame("SpeciesID"="species", "Kappa"=0, "TSS"=0, "ROC"=0)

# for loop through all species
for(spID in speciesSub){ try({
  
  ## Load probability maps 
  load(file=paste0(here::here(), "/results/_SDMs/SDM_biomod_", spID, ".RData")) #biomod_list
  temp_validation <- biomod_list$validation[,str_detect( colnames(biomod_list$validation), "Testing.data")]
  temp_validation <- temp_validation[,str_detect( colnames(temp_validation), "EMcaByTSS")]
  
  data_eval <- rbind(data_eval, c(spID, temp_validation))
  
  print(paste0(spID, " successfully loaded."))
}, silent=T)}

data_eval <- data_eval %>% filter(SpeciesID!="species") 
data_eval$Kappa <- as.numeric(data_eval$Kappa)
data_eval$ROC <- as.numeric(data_eval$ROC)
data_eval$TSS <- as.numeric(data_eval$TSS)

data_eval; str(data_eval)
summary(data_eval)

write.csv(data_eval, paste0(here::here(), "/results/Model_evaluation_", Taxon_name, ".csv"), row.names = F)


