#- - - - - - - - - - - - - - - - -#
#   Data preparation for SDMs     #
#                                 #
#     author: Romy Zeiss          #
#       date: 18.08.2021          #
#- - - - - - - - - - - - - - - - -#

# - - - - - - - - - - - - - - - - - - -
## Data from sWorm ####
sworm_occ <- readr::read_csv(file=paste0(here::here(), "/data/SppOccData_sWorm_v2.csv"))
sworm_site <- readr::read_csv(file=paste0(here::here(), "/data/SiteData_sWorm_v2.csv"))

sworm <- dplyr::full_join(sworm_occ, sworm_site)

sworm <- sworm %>% filter(Abundance>0)

sworm$datasource <- "sWorm"

rm(sworm_occ, sworm_site)

# - - - - - - - - - - - - - - - - - - -
## Data from GBIF ####

load(file=paste0(here::here(), "/results/Occurrences_GBIF_Crassiclitellata.RData")) #dat

gbif <- dat$data[,c("species", "decimalLatitude", "decimalLongitude")] %>%
  mutate(datasource = "GBIF")

rm(dat)

# - - - - - - - - - - - - - - - - - - -
## Data from Edaphobase ####

edapho <- readr::read_csv(file=paste0(here::here(), "/data/Edaphobase_download_24-Feb-2021_Lumbricidae_Europe.csv"))

#!!! Manually: We already added one missing "Valid taxon" for Helodrilus sp.

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

# - - - - - - - - - - - - - - - - - - -
## Save data ####
write.csv(data, file=paste0(here::here(), "/data/Earthworm_occurrence_GBIF-sWorm-Edaphobase.csv"),
          row.names = F)

#!!! In ArcGIS: Select only points that fall into German shapefile.
#               Save them as "Earthworm_occurrence_Germany.txt"

# - - - - - - - - - - - - - - - - - - -
## Count occurrences per species & datasource ####
count.data <- data %>% group_by(datasource, species) %>% count() %>%
  pivot_wider(names_from = datasource, values_from = n)
count.data$TotalOcc <- rowSums(count.data[,2:4])

count.data

ggplot(data=count.data, aes(x=TotalOcc, y=species))+
  geom_bar(stat="identity")+
  theme(axis.text = element_text(size=4))

# - - - - - - - - - - - - - - - - - - -
## CoordinateCleaner ####
# flag problems with coordinates
dat_cl <- data.frame(data)
flags <- CoordinateCleaner::clean_coordinates(x = dat_cl, lon = "longitude", lat = "latitude",
                                              species = "species", tests = c("capitals", "centroids", "equal", "gbif", "zeros", "seas"), #normally: test "countries"
                                              country_ref = rnaturalearth::ne_countries("small"), 
                                              country_refcol = "iso_a3")
sum(flags$.summary) #those not flagged! = 75388

# remove flagged records from the clean data (i.e., only keep non-flagged ones)
dat_cl <- dat_cl[flags$.summary, ]

## Plot flagged records
world.inp <- map_data("world")

ggplot() + 
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") + 
  xlim(min(data$longitude, na.rm = T), max(data$longitude, na.rm = T)) + 
  ylim(min(data$latitude, na.rm = T), max(data$latitude, na.rm = T)) + 
  geom_point(data = data, aes(x = longitude, y = latitude), colour = "darkred", size = 1, show.legend = T) +
  geom_point(data = dat_cl, aes(x = longitude, y = latitude), colour = "darkgreen", size = 1, show.legend = T) + 
  coord_fixed() + 
  scale_color_manual(name='CoordinateCleaner',
                     values=c('RawRecords = red'='darkred', 'CleanRecords = green'='darkgreen'))+ 
  theme_bw() + theme(axis.title = element_blank())

# - - - - - - - - - - - - - - - - - - -
## Some structuring ####
# add species ID
dat_cl <- dplyr::right_join(dat_cl, speciesNames[,c("Species", "SpeciesID")], 
                            by=c("species" = "Species"))

# create x and y column
dat_cl$x <- dat_cl$longitude
dat_cl$y <- dat_cl$latitude

# remove NA in coordinates
dat_cl <- dat_cl[complete.cases(dat_cl$x),]
dat_cl <- dat_cl[complete.cases(dat_cl$y),]

# - - - - - - - - - - - - - - - - - - -
## Save clean data ####
write.csv(dat_cl, file=paste0(here::here(), "/data/Earthworm_occurrence_GBIF-sWorm-Edaphobase.csv"),
          row.names = F)

# load number of records during cleaning process
df.cleaning <- read.csv(file=paste0(here::here(), "/results/NoRecords_cleaning_", Taxon_name, ".csv"))

df.cleaning <- df.cleaning %>% add_row(CleaningStep=c("mergedRecords","coordinateCleaner"), 
                                       NumberRecords=c(nrow(data), nrow(dat_cl)))

df.cleaning

# save updated number of records during cleaning process
write.csv(df.cleaning, file=paste0(here::here(), "/results/NoRecords_cleaning_", Taxon_name, ".csv"), row.names = F)

# - - - - - - - - - - - - - - - - - - -
## OLD CODE ####
# to add grids per species (prepared in ArcGIS) with other occurrence data

# gbif.wd <- "I:/eie/==PERSONAL/Macroecology/Students/Jessica/Grid/GBIF_Data"
# species.folders <- list.dirs(path=gbif.wd, recursive=F)[-46]
# 
# gbif <- tibble::tibble(countryCode="DE", decimalLatitude=1.1, decimalLongitude=1.1, 
#                        family="Family", genus="Genus", specificEpithet="Species")[0,]
# 
# for(i in 1:length(species.folders)){
#   temp.wd <- species.folders[i]
#   
#   temp.data <- read.delim(paste0(temp.wd, "/occurrence.txt")) 
#   temp.data <- temp.data[,c("countryCode", "decimalLatitude", 
#                             "decimalLongitude", "family", "genus", "specificEpithet")]
#   
#   gbif <- dplyr::full_join(gbif, temp.data)
# }
# 
# gbif$species <- paste(gbif$genus, gbif$specificEpithet)
# gbif$datasource <- "GBIF"
# 
# gbif
# #write.csv(gbif, here::here("data", "GBIF_Earthworms_combined.csv"))
# 
# rm(temp.data, gbif.wd, temp.wd, species.folders, i)

