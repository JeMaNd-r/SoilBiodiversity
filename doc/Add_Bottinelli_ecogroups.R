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

# read earthworm file from Phillips et a. 2019
Phillips <- read.csv(file="../data/SppOccData_sWorm_v2.csv")

# read Edaphobase data
Edapho <- read.csv(file="../data/Edaphobase_download_24-Feb-2021_Lumbricidae_Europe.csv")

#- - - - - - - - - - - - - - - - - - - - -
## Add ecological groups assigned in Phillips et al. 2019 ####

Phillips <- Phillips %>% 
  dplyr::select(SpeciesBinomial, Ecological_group, Family, Genus) %>% 
  unique() %>%
  rename("Ecogroup_Phillips" = Ecological_group,
         "Species" = SpeciesBinomial) %>%
  filter(!is.na(Species))

Phillips[Phillips$Ecogroup_Phillips=="Unknown", "Ecogroup_Phillips"] <- NA

ew_list2 <- dplyr::full_join(ew_list, Phillips)

#- - - - - - - - - - - - - - - - - - - - -
## Add species from Edaphobase ####

# extract the first 2 words of the species' names
Edapho$Species <- word(Edapho$Valid.taxon, 1, 2)
Edapho$Species <- sub(",*", "\\1", Edapho$Species)
Edapho$Species <- gsub(",", "", Edapho$Species)

Edapho <- Edapho %>% dplyr::select(Species) %>% unique() %>%
  filter(!str_detect(Species, " [:upper:]"))

Edapho$Genus <- word(Edapho$Species, 1, 1)

ew_list2 <- dplyr::full_join(ew_list2, Edapho %>% filter(!Species %in% ew_list2$Species))

# - - - - - - - - - - - - - - - - - -
## Add missing species ID ####
ew_list2$SpeciesID <- paste0(substr(word(ew_list2$Species, 1), 1, 4), "_", 
                               substr(word(ew_list2$Species, 2), 1, 4))

# - - - - - - - - - - - - - - - - - -
## Add Bottinelli to species list ####

## Prepare data for merge
# add SpeciesID column to Bottinelli
Bottinelli$SpeciesID <- paste0(substr(Bottinelli$Genus, 1,4), "_", 
                               substr(Bottinelli$Species, 1,4))

# check how different the SpeciesID's are
setdiff(Bottinelli$SpeciesID, ew_list2$SpeciesID)

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

setdiff(Bottinelli$Species, ew_list2$Species)
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
ew_list2 <- dplyr::left_join(ew_list2, ew_data)

#- - - - - - - - - - - - - - - - - - - - -
## Create summarizing eco-group column ####
ew_list2$Ecogroup <- NA

for(i in 1:nrow(ew_list2)){
  temp_row <- ew_list2[i,]
  
  temp_group <- temp_row$Ecogroup_Bottinelli
  
  # if there is no Ecogroup_Bottinelli, take Phillips
  if(is.na(temp_group)) temp_group <- temp_row$Ecogroup_Phillips
  
  ew_list2[i,"Ecogroup"] <- temp_group
}


#- - - - - - - - - - - - - - - - - - - - -
## Add family if missing ####

Phylogeny <- ew_list2 %>% 
  dplyr::select(Group_name, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  unique() %>%
  filter(!is.na(Family), !is.na(Kingdom))

ew_list2 <- dplyr::inner_join(ew_list2, Phylogeny)

#- - - - - - - - - - - - - - - - - - - - -
## Save expanded ew_list ####
ew_list2 <- ew_list2 %>% arrange("Order", "Family", "Genus", "Species")

write.csv(ew_list2, file="../data/Species_list_Crassiclitellata.csv", row.names = F)

