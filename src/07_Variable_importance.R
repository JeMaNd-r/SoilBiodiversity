#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Variable importance               #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(biomod2) # also to create pseudo-absences

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

#- - - - - - - - - - - - - - - - - - - - -
## Variable importance ####
#- - - - - - - - - - - - - - - - - - - - -

## Calculate variable importance (VI) ####

# create result data frame
var_imp <- data.frame("Predictor"= c("Aridity", "MAP", "MAP_Seas", "MAT", 
                                     "MAT_Seas", "Snow", "Agriculture", "Dist_Urban",
                                     "Forest_Coni", "Forest_Deci", "NDVI", 
                                     "Pastures", "Pop_Dens", "Shrubland", "Aspect",
                                     "Dist_Coast", "Dist_River", "Elev", 
                                     "Slope", "CEC", "Clay.Silt", "Cu", "Hg",
                                     "Moisture", "N", "P", "pH", "SOC", "SoilT"),
                      "biomod"=NA, "Species"=NA)

for(spID in unique(speciesNames[speciesNames$NumCells_5km >= 10,]$SpeciesID)){ try({
  
  print("=====================================")
  print(spID)
  
  load(file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_biomod_", spID, ".RData")) #biomod_list
  temp_vi <- biomod_list$varImp
  
  # round variable importance to 3 decimals
  temp_vi[,"biomod"] <- round(temp_vi[,"biomod"], 3)
  
  # add species name
  temp_vi$Species <- spID
  
  var_imp <- rbind(var_imp, temp_vi)
  
  rm(temp_vi, biomod_list)
  
  #var_imp
  
})}

var_imp <- var_imp[!is.na(var_imp$Species),] %>% unique()
rownames(var_imp) <- 1:nrow(var_imp)

var_imp
str(var_imp)

## Save ####
write_csv(var_imp, file=paste0(here::here(), "/results/Variable_importance_biomod_", Taxon_name, ".csv"))



