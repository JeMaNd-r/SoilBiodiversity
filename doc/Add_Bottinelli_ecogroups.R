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

ew_list <- tibble::tibble("Group_name" = "Earthworms", 
                      "Species" = "species", "SpeciesID" = NA, "Fct_group" = NA, "Comment" = NA)[0,]
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

ew_list <- ew_list %>% full_join(Phillips[,c("Species", "Ecogroup_Phillips")] %>% unique())
ew_list

#- - - - - - - - - - - - - - - - - - - - -
## Add species from GBIF ####

ew_list <- ew_list %>%
  full_join(GBIF %>% dplyr::select(species, #acceptedScientificName, taxonomicStatus
                                   ) %>% unique(),
            by=c("Species" = "species"))
ew_list

#- - - - - - - - - - - - - - - - - - - - -
## Add species from Edaphobase ####

# extract the first 2 words of the species' names
Edapho$Species <- word(Edapho$Valid.taxon, 1, 2)
Edapho$Species <- sub(",*", "\\1", Edapho$Species)
Edapho$Species <- gsub(",", "", Edapho$Species)

Edapho <- Edapho %>% dplyr::select(Species) %>% unique() %>%
  filter(!str_detect(Species, " [:upper:]"))

#Edapho$Genus <- word(Edapho$Species, 1, 1)

ew_list <- dplyr::full_join(ew_list, Edapho %>% filter(!Species %in% ew_list$Species), 
                            by=c("Species") %>% unique())
ew_list

#- - - - - - - - - - - - - - - - - - - - -
## Define SpeciesID while accounting for synonyms ####

# delete species with Identifier name as species name (e.g., Lumbricus Linneous)
ew_list <- ew_list %>% filter(!str_detect(Species, " [:upper:]"))

#tax_key <- rgbif::name_suggest(q=Taxon_name, rank=Taxon_rank)

# # fix synonyms using ITIS (Integrated Taxonomic Information System) in taxize package
# temp <- taxize::synonyms(ew_list$Species, db="itis")
# synonym_ids <- grep(pattern = "acc_name", temp) #is this the optimal solution?
# all_names <- lapply(temp[synonym_ids], cbind)
# all_names <- do.call(rbind, all_names)
# 
# all_names$Species <- rownames(all_names)
# all_names$Acc_name <- word(all_names$acc_name, 1, 2)
# all_names$Synonyms <- word(all_names$syn_name, 1, 2)
# 
# ew_list <- ew_list %>% full_join(all_names[all_names$Synonyms %in% ew_list$Species, c("Acc_name", "Synonyms", "Species")])

# fix misspellings
speciesA <- unique(ew_list$Species)
temp <- taxize::gnr_resolve(speciesA, best_match_only = TRUE, canonical = TRUE)
ew_list <- ew_list %>% full_join(temp[,c("user_supplied_name", "matched_name2")], by = c("Species" = "user_supplied_name")) %>%
  rename("Acc_name" = matched_name2)
 
# # alternative approach to only keep accepted names
# temp <- taxize::itis_acceptname(taxize::get_tsn(speciesA))
#  
# #this provides nicer output and can be used to drop unknown species, AND keep synonims.
# taxas <- taxize::tax_name(query = speciesA, get = "species", verbose = TRUE)
# 
# ew_list <- ew_list %>% full_join(taxas %>% mutate("tax_name" = species) %>% dplyr::select(query, tax_name), 
#                                  by = c("Species" = "query"))

# - - - - - - - - - - - - - - - - - -
## Add missing species ID ####
ew_list$SpeciesID <- paste0(substr(word(ew_list$Acc_name, 1), 1, 5), "_", 
                            substr(word(ew_list$Acc_name, 2), 1, 4))
# replace "..._NA" with speciesID created by original species name
for(i in 1:nrow(ew_list)){
  if(stringr::str_detect(ew_list[i,]$SpeciesID, "_[:upper:]")) ew_list[i,]$SpeciesID <- paste0(substr(word(ew_list[i,]$Species, 1), 1, 5), "_", 
                                                                         substr(word(ew_list[i,]$Species, 2), 1, 4))
}

length(unique(ew_list$SpeciesID))

# - - - - - - - - - - - - - - - - - -
## Add Bottinelli to species list ####
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


# fix misspellings in Bottinelli
speciesB <- unique(Bottinelli$Species_Bouche)
temp <- taxize::gnr_resolve(speciesB, best_match_only = TRUE, canonical = TRUE)
Bottinelli <- Bottinelli %>% full_join(temp[,c("user_supplied_name", "matched_name2")], by = c("Species_Bouche" = "user_supplied_name")) %>%
  rename("Acc_name" = matched_name2)

# add SpeciesID column to Bottinelli
Bottinelli$SpeciesID <- paste0(substr(word(Bottinelli$Acc_name, 1), 1, 5), "_", 
                               substr(word(Bottinelli$Acc_name, 2), 1, 4))

# replace "..._NA" with speciesID created by original species name
for(i in 1:nrow(Bottinelli)){
  if(stringr::str_detect(Bottinelli[i,]$SpeciesID, "_[:upper:]")) Bottinelli[i,]$SpeciesID <- paste0(substr(word(Bottinelli[i,]$Species_Bouche, 1), 1, 4), "_", 
                                                                                                     substr(word(Bottinelli[i,]$Species_Bouche, 2), 1, 4))
}

# check how different the SpeciesID's are
setdiff(Bottinelli$SpeciesID, ew_list$SpeciesID)
setdiff(ew_list$SpeciesID, Bottinelli$SpeciesID)

# try with [genus] [species] names
Bottinelli$Species_Bottinelli <- paste0(Bottinelli$Genus, " ", Bottinelli$Species)

# keep only unique rows based on some criteria
ew_data <- Bottinelli[0,c("SpeciesID", "Percent_Epigeic", "Percent_Anecic",
                        "Percent_Endogeic")]
for(i in unique(Bottinelli$SpeciesID)){
  data <- Bottinelli[Bottinelli$SpeciesID==i,]
  temp.species <- i
  temp.subspecies <- NA
  
  # take mean if having multiple Bouche names (= rows) for that species
  if(nrow(data)!=1){
    temp.species <- !stringr::str_detect(data$Subspecies_Bouche, data$Species)
  }  
  # average percentages
  temp.data <- data %>% dplyr::select(SpeciesID, Percent_Epigeic, Percent_Anecic, Percent_Endogeic) %>% 
    group_by(SpeciesID) %>% summarise_all(mean)

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

phylogeny <- taxize::classification(ew_list$Acc_name, db="itis", rows=NA) #rows=NA means all rows
phylogeny_list <- lapply(1:length(phylogeny), function(x) {
  if(!is.na(phylogeny[[x]])) pivot_wider(phylogeny[[x]] %>% dplyr::select(-id),names_from = rank, values_from=name)})
phylogeny_df <- bind_rows(phylogeny_list)

ew_list <- ew_list %>% full_join(phylogeny_df, by=c("Acc_name" = "species"))

ew_list[is.na(ew_list$genus),"genus"] <- word(ew_list[is.na(ew_list$genus),]$Acc_name, 1, 1)

ew_list %>% inner_join(ew_list %>% dplyr::select(kingdom, phylum, order, family, genus) %>% unique(), by=c("kingdom", "phylum", "order", "family", "genus"))


ew_list$Group_name <- "Earthworms"

# #... TODO: use taxize to get information :)
# Phylogeny <- ew_list %>% 
#   dplyr::select(Group_name, Kingdom, Phylum, Class, Order, Family, Genus) %>%
#   unique() %>%
#   filter(!is.na(Family), !is.na(Kingdom))
# 
# ew_list <- dplyr::full_join(ew_list, Phylogeny)

#- - - - - - - - - - - - - - - - - - - - -
## Save expanded ew_list ####
ew_list <- ew_list[order(ew_list$Acc_name),]
ew_list

write.csv(ew_list, file="../data/Species_list_Crassiclitellata.csv", row.names = F)

# BY HAND:
# add remaining (missing) classification based on already present genus, or with internet.
# remove duplicates
# add wrong speciesID (Leptogaster in Lennogaster)
# make Firzingeria depressa classification to Firzingeria (not otehr genus)
# result: 351 unique species names ...


