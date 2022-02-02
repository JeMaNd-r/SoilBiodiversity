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
Bottinelli$Species_Bouche <- NA
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
    Bottinelli[i,]$Species_Bouche <- paste0(Bottinelli[i,]$Subgenus_Bouche, " ", Bottinelli[i,]$Subspecies_Bouche)  
    }
  }
  
  # remove double space in some species' names
  Bottinelli[i,]$Species_Bouche <- stringr::str_squish(Bottinelli[i,]$Species_Bouche)
}


# try with [genus] [species] names
Bottinelli$Species <- paste0(Bottinelli$Genus, " ", Bottinelli$Species)

setdiff(Bottinelli$Species, ew_list$Species)
setdiff(Bottinelli$Species_Bouche, ew_list$Species)

# keep only unique rows based on some criteria
ew_data <- Bottinelli[0,c("Species", "Percent_Epigeic", "Percent_Anecic",
                        "Percent_Endogeic")]
for(i in unique(Bottinelli$Species)){
  data <- Bottinelli[Bottinelli$Species==i,]
  temp.species <- i
  temp.subspecies <- NA
  
  # take mean if having multiple Bouche names (= rows) for that species
  if(nrow(data)!=1){
    temp.species <- !stringr::str_detect(data$Subspecies_Bouche, data$Species)
  }  
  # average percentages
  temp.data <- data %>% select(Species, Percent_Epigeic, Percent_Anecic, Percent_Endogeic) %>% 
    group_by(Species) %>% summarise_all(mean)

  # take most frequently used Botinelli group for this species
  temp.data$Ecogroup_Bottinelli <- names(table(data$Ecogroup_Bottinelli)[
    table(data$Ecogroup_Bottinelli)==max(table(data$Ecogroup_Bottinelli))])[1]
  
  ew_data <- rbind(ew_data, temp.data)
}

ew_data

## Add ecological groups to ew_list
ew_list2 <- dplyr::left_join(ew_list, ew_data)

## Save expanded ew_list
write.csv(ew_list2, file="../data/Species_list_Crassiclitellata.csv", row.names = F)
