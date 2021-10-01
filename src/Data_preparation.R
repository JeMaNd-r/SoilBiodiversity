#- - - - - - - - - - - - - - - - -#
#   Data preparation for SDMs     #
#                                 #
#     author: Romy Zeiss          #
#       date: 18.08.2021          #
#- - - - - - - - - - - - - - - - -#

#library(here)  #instead of setwd()
library(tidyverse)

setwd("I:/eie/==PERSONAL/RZ SoilBON/SoilBiodiversity")

# - - - - - - - - - - - - - - - - - - -
## Data from sWorm ####
sworm_occ <- readr::read_csv("data/SppOccData_sWorm_v2.csv")
sworm_site <- readr::read_csv("data/SiteData_sWorm_v2.csv")

sworm <- dplyr::full_join(sworm_occ, sworm_site)

sworm <- sworm %>% filter(Abundance>0)

sworm$datasource <- "sWorm"

rm(sworm_occ, sworm_site)

# - - - - - - - - - - - - - - - - - - -
## Data from GBIF ####
gbif.wd <- "I:/eie/==PERSONAL/Macroecology/Students/Jessica/Grid/GBIF_Data"
species.folders <- list.dirs(path=gbif.wd, recursive=F)[-46]

gbif <- tibble::tibble(countryCode="DE", decimalLatitude=1.1, decimalLongitude=1.1, 
                       family="Family", genus="Genus", specificEpithet="Species")[0,]

for(i in 1:length(species.folders)){
  temp.wd <- species.folders[i]

  temp.data <- read.delim(paste0(temp.wd, "/occurrence.txt")) 
  temp.data <- temp.data[,c("countryCode", "decimalLatitude", 
                            "decimalLongitude", "family", "genus", "specificEpithet")]
  
  gbif <- dplyr::full_join(gbif, temp.data)
}

gbif$species <- paste(gbif$genus, gbif$specificEpithet)
gbif$datasource <- "GBIF"

gbif
#write.csv(gbif, here::here("data", "GBIF_Earthworms_combined.csv"))

rm(temp.data, gbif.wd, temp.wd, species.folders, i)

# - - - - - - - - - - - - - - - - - - -
## Data from Edaphobase ####

edapho <- readr::read_csv("data/Edaphobase_download_24-Feb-2021_Lumbricidae_Europe.csv")

#!!! Manually: We added the missing "Valid taxon" for Helodrilus sp.

# extract the first 2 words of the species' names
edapho$species <- word(edapho$'Valid taxon', 1, 2)
edapho$species <- sub(",*", "\\1", edapho$species)

edapho$datasource <- "Edaphobase"

# have a look at the data
edapho[,c("species", "Latitude", "Longitude", "datasource")]
  
# - - - - - - - - - - - - - - - - - - -
## Merge all together ####

data <- dplyr::full_join(sworm, gbif, 
                         by=c("SpeciesBinomial"="species", 
                              "Latitude_decimal_degrees"="decimalLatitude", 
                              "Longitude_decimal_degrees"="decimalLongitude", 
                              "datasource"))

data <- dplyr::full_join(data, edapho, 
                 by=c("SpeciesBinomial"="species", 
                      "Latitude_decimal_degrees"="Latitude", 
                      "Longitude_decimal_degrees"="Longitude",
                      "datasource"))

data <- tibble::tibble(species=data$SpeciesBinomial, 
                       latitude=data$Latitude_decimal_degrees, 
                       longitude=data$Longitude_decimal_degrees, 
                       datasource=data$datasource)

data
data <- data[complete.cases(data$longitude, data$latitude),]

data$OBJECTID <- 1:nrow(data) 

## Save data ####
write.csv(data, "data/Earthworm_occurrence_GBIF-sWorm-Edaphobase.csv",
          row.names = F)

#!!! In ArcGIS: Select only points that fall into German shapefile.
#               Save them as "Earthworm_occurrence_Germany.txt"

## Count occurrences per species & datasource ####
count.data <- data %>% group_by(datasource, species) %>% count()
View(count.data)
