# - - - - - - - - - - - - - - - - - - #
#                                     #
#   Assign ecogroups of Bottinelli    #
#       author: Romy Zeiss            #
#                                     #
# - - - - - - - - - - - - - - - - - - #

# reference: Bottinelli et al. 2020, Geoderma

library(tidyverse)

## Load data
# read classification of Bottinelli
Bottinelli <- read.csv(file="Bottinelli_2020_Earthworm_classification.csv")

# read earthworm species list
ew_list <- read.csv(file="../doc/Species_list_Crassiclitellata.csv")


## Prepare data for merge
# add SpeciesID column to Bottinelli
Bottinelli$SpeciesID <- paste0(substr(Bottinelli$Genus, 1,4), "_", 
                               substr(Bottinelli$Species, 1,4))

# check how different the SpeciesID's are
setdiff(Bottinelli$SpeciesID, ew_list$SpeciesID)

## merge Buoche species and genus name
for(i in 1:nrow(Bottinelli)){
  # check if Subgenus_Bouche is NA
  if(is.na(Bottinelli[i,]$Subgenus_Bouche)){
    # check if Subspecies_Bouche is NA
    if(is.na(Bottinelli[i,]$Subspecies_Bouche)){
    Bottinelli[i,]$Species_Bouche <-paste0(Bottinelli[i,]$Genus, " ", Bottinelli[i,]$Species)
    }else{
      Bottinelli[i,]$Species_Bouche <-paste0(Bottinelli[i,]$Genus, " ", Bottinelli[i,]$Subspecies_Bouche)
    }
  }else{
    # check if Subspecies_Bouche is NA
    if(is.na(Bottinelli[i,]$Subspecies_Bouche)){
      Bottinelli[i,]$Species_Bouche <-paste0(Bottinelli[i,]$Subgenus_Bouche, " ", Bottinelli[i,]$Species)
    }else{
    Bottinelli[i,]$Species_Bouche <-paste0(Bottinelli[i,]$Subgenus_Bouche, " ", Bottinelli[i,]$Subspecies_Bouche)  
    }
  }
}

# try with [genus] [species] names
Bottinelli$Species <- paste0(Bottinelli$Genus, " ", Bottinelli$Species)

setdiff(Bottinelli$Species, ew_list$Species)
setdiff(Bottinelli$Species_Bouche, ew_list$Species)


# keep only unique rows based on some criteria
for(i in unique(Bottinelli$Species)){
  data <- Bottinelli[Bottinelli$Species==i,]
  temp.species <- NA
  temp.subspecies <- NA
  
  if(nrow(data)==1){next()}
  temp.species <- !stringr::str_detect(data$Subspecies_Bouche, data$Species)
  
  if(summary(temp.species==TRUE)[2]<=2){
    temp.subspecies <- stringr::str_detect(data$Variety_Bouche, "int")
  }
  
  data %>% filter(temp.species)
  
  aggregate(data, by=list(data$Family, data$Genus, data$Species), mean)
  
}


## Add ecological groups to ew_list
ew_list2 <- dplyr::full_join(ew_list, Bottinelli)
