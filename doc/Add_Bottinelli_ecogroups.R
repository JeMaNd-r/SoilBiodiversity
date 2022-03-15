# - - - - - - - - - - - - - - - - - - #
#                                     #
#   Assign ecogroups of Bottinelli    #
#       author: Romy Zeiss            #
#                                     #
# - - - - - - - - - - - - - - - - - - #

# reference: Bottinelli et al. 2020, Geoderma

library(tidyverse)
library(rgbif)

## Load data
# read classification of Bottinelli
Bottinelli <- read.csv(file="Bottinelli_2020_Earthworm_classification.csv")

# read earthworm species list
#ew_list <- read.csv(file="../doc/Species_list_Crassiclitellata.csv")

# read raw GBIF data
load(file=paste0("../results/RawOccurrences_", Taxon_name, ".RData")) #dat
GBIF <- dat$data; rm(dat)

# read earthworm file from Phillips et a. 2019
Phillips <- read.csv(file="../data/SppOccData_sWorm_v2.csv")

# read Edaphobase data
Edapho <- read.csv(file="../data/Edaphobase_download_24-Feb-2021_Lumbricidae_Europe.csv")


#- - - - - - - - - - - - - - - - - - - - -
## Empty earthworm list ####

ew_list <- tibble::tibble("Group_name" = "Earthworms", "Kingdom" = NA, "Phylum" = NA, 
                      "Class" = NA, "Order" = NA, "Family" = "family", "Genus" = "genus",
                      "Species" = "species", "SpeciesID" = NA, "Fct_group" = NA, "Comment" = NA,
                      "Species_synonym" = NA)[0,]
ew_list

#- - - - - - - - - - - - - - - - - - - - -
## Load species and ecological groups assigned in Phillips et al. 2019 ####

Phillips <- Phillips %>% 
  dplyr::select(SpeciesBinomial, Ecological_group, Family, Genus) %>% 
  unique() %>%
  rename("Ecogroup_Phillips" = Ecological_group,
         "Species" = SpeciesBinomial) %>%
  filter(!is.na(Species))

Phillips[Phillips$Ecogroup_Phillips=="Unknown", "Ecogroup_Phillips"] <- NA

ew_list <- ew_list %>% full_join(Phillips)
ew_list

#- - - - - - - - - - - - - - - - - - - - -
## Add species from GBIF ####

GBIF <- GBIF %>%
  full_join(GBIF %>% dplyr::select(species, #acceptedScientificName, taxonomicStatus,
                                   kingdom, phylum, order, family, genus) %>% unique(),
            by=c("Species" = "species", "Kingdom" = "kingdom", "Phylum" = "phylum",
                 "Order" = "order", "Family" = "family", "Genus" = "genus"))

ew_list <- ew_list %>% full_join(GBIF %>% dplyr::select(species) %>% unique(), 
                                 by = c("Species" = "species"))
ew_list

#- - - - - - - - - - - - - - - - - - - - -
## Add species from Edaphobase ####

# extract the first 2 words of the species' names
Edapho$Species <- word(Edapho$Valid.taxon, 1, 2)
Edapho$Species <- sub(",*", "\\1", Edapho$Species)
Edapho$Species <- gsub(",", "", Edapho$Species)

Edapho <- Edapho %>% dplyr::select(Species) %>% unique() %>%
  filter(!str_detect(Species, " [:upper:]"))

Edapho$Genus <- word(Edapho$Species, 1, 1)

ew_list <- dplyr::full_join(ew_list, Edapho %>% filter(!Species %in% ew_list$Species), by=c("Species")) %>% unique()
ew_list

# - - - - - - - - - - - - - - - - - -
## Add missing species ID ####
ew_list$SpeciesID <- paste0(substr(word(ew_list$Species, 1), 1, 4), "_", 
                               substr(word(ew_list$Species, 2), 1, 4))

#- - - - - - - - - - - - - - - - - - - - -
## Define SpeciesID while accounting for synonyms ####

# delete species with Identifier name as species name (e.g., Lumbricus Linneous)
ew_list <- ew_list %>% filter(!str_detect(Species, " [:upper:]"))

#tax_key <- rgbif::name_suggest(q=Taxon_name, rank=Taxon_rank)

# fix synonyms using ITIS (Integrated Taxonomic Information System) in taxize package
temp <- taxize::synonyms(ew_list$Species, db="itis")
synonym_ids <- grep(pattern = "acc_name", temp) #is this the optimal solution?
#accepted_names <- unlist(cbind(lapply(temp[synonym_ids], '[', "acc_name"), lapply(temp[synonym_ids], '[', "syn_name")), use.names = FALSE)
all_names <- lapply(temp[synonym_ids], cbind)
#ew_list$Species[synonym_ids] <- accepted_names
ew_list <- ew_list %>% full_join(all_names %>% dplyr::select(acc_name, syn_name),
                                 by=c("Species" = "acc_name")) 

# fix misspellings
species2 <- unique(species)
temp <- taxize::gnr_resolve(species2, best_match_only = TRUE, canonical = TRUE)
temp

species2 <- temp$matched_name2
# here We will need to recover repeated species in an eficient way, as the are dropped.

# keep accepted names only
taxize::itis_acceptname(get_tsn(species2))
vapply(x, itis_acceptname, "")

out <- list()
for(i in 1:length(species2)){
  out[[i]] <- itis_acceptname(get_tsn(species2[i]))
}
#All accepted, wich is not what I want.

#this provides nicer output and can be used to drop unknown species, AND keep synonims.
taxas <- tax_name(query = species2, get = "species", verbose = TRUE)
#fails because not all has species. in a for loop will work.
out <- list()
for(i in 1:length(species2)){
  out[[i]] <- tax_name(species2[i], get = "species")
}
out2 <- plyr::ldply(out, data.frame)
species2[-which(is.na(out2$species))]

#note, using genus do not work, because all has genus now.
taxas <- tax_name(query = species2, get = "genus", verbose = TRUE)


# - - - - - - - - - - - - - - - - - -
## Add Bottinelli to species list ####

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
  temp.data <- data %>% dplyr::select(Species, Percent_Epigeic, Percent_Anecic, Percent_Endogeic) %>% 
    group_by(Species) %>% summarise_all(mean)

  # take most frequently used Botinelli group for this species
  temp.data$Ecogroup_Bottinelli <- names(table(data$Ecogroup_Bottinelli)[
    table(data$Ecogroup_Bottinelli)==max(table(data$Ecogroup_Bottinelli))])[1]
  
  ew_data <- rbind(ew_data, temp.data)
}

ew_data

## Add ecological groups to ew_list
ew_list <- dplyr::left_join(ew_list, ew_data)

#- - - - - - - - - - - - - - - - - - - - -
## Create summarizing eco-group column ####
ew_list$Ecogroup <- NA

for(i in 1:nrow(ew_list)){
  temp_row <- ew_list[i,]
  
  temp_group <- temp_row$Ecogroup_Bottinelli
  
  # if there is no Ecogroup_Bottinelli, take Phillips
  if(is.na(temp_group)) temp_group <- temp_row$Ecogroup_Phillips
  
  ew_list[i,"Ecogroup"] <- temp_group
}


#- - - - - - - - - - - - - - - - - - - - -
## Add family if missing ####

Phylogeny <- ew_list %>% 
  dplyr::select(Group_name, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  unique() %>%
  filter(!is.na(Family), !is.na(Kingdom))

ew_list <- dplyr::inner_join(ew_list, Phylogeny)

#- - - - - - - - - - - - - - - - - - - - -
## Save expanded ew_list ####
ew_list <- ew_list %>% arrange("Order", "Family", "Genus", "Species")

write.csv(ew_list, file="../data/Species_list_Crassiclitellata.csv", row.names = F)

