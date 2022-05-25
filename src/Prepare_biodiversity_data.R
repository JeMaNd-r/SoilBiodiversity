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

sworm <- sworm %>% filter(Abundance>0, !is.na(SpeciesBinomial))

sworm$datasource <- "sWorm"

rm(sworm_occ, sworm_site)

# - - - - - - - - - - - - - - - - - - -
## Data from GBIF ####

dat <- read.csv(file=paste0(here::here(), "/results/Occurrences_GBIF_Crassiclitellata.csv")) #dat

gbif <- dat[,c("species", "decimalLatitude", "decimalLongitude")] %>%
  mutate(datasource = "GBIF")

rm(dat)

# - - - - - - - - - - - - - - - - - - -
## Data from Edaphobase ####

edapho <- readr::read_csv(file=paste0(here::here(), "/data/Edaphobase_download_24-Feb-2021_Lumbricidae_Europe.csv"))

#!!! Manually: We already added one missing "Valid taxon" for Helodrilus sp.

# extract the first 2 words of the species' names
edapho$species <- word(edapho$"Valid taxon", 1, 2)
edapho$species <- sub(",*", "\\1", edapho$species)
edapho$species <- gsub(",", "", edapho$species)

# exclude species with only Genus name (followed by identifier's name)
edapho <- edapho %>% filter(!str_detect(species, " [:upper:]"))

# remove observations before 1990
edapho$date <- as.Date(edapho$"Observation date (Sampling event)", format="%m/%d/%y") 
edapho <- edapho %>% filter(date >= "1990-01-01")

edapho$datasource <- "Edaphobase"

# have a look at the data
edapho[,c("species", "Latitude", "Longitude", "datasource")] #n=13.473

# - - - - - - - - - - - - - - - - - - -
## Data from SoilReCon project (Portugal) ####
  
recon <- readr::read_csv(file=paste0(here::here(), "/data/SoilReCon_earthworms_clean.csv"))

# remove non-species level and NA species
recon <- recon[recon$Species!="Octolasion sp." & !is.na(recon$Species),]

# filter relevant columns
recon <- recon %>% dplyr::select(Species, POINT_X, POINT_Y) %>%
  mutate("datasource"="SoilReCon")

# - - - - - - - - - - - - - - - - - - -
## Data from Jérome Matthieu (Jema) ####

jema <-  read.csv(file=paste0(here::here(),"/data/worm_spd_europe_jerome.csv"), sep=";")
# nrow=27800

jema$Species <- jema$species_name

# rename "A. caliginosa trapezoides"
jema[jema$Species=="A. caliginosa trapezoides","Species"] <- "Aporrectodea caliginosa"

# make upper letter
jema$Species <- stringr::str_to_sentence(jema$Species)

# keep only species, not subspecies
jema$Species <- word(jema$Species, 1, 2)

# remove - from species names
jema$Species <- gsub("[[:punct:]]", "",jema$Species)

# remove sp species
jema <- jema %>% filter(jema$Species %in% jema$Species[!stringr::str_detect(jema$Species, "sp[p]?[:blank:]")])

jema$lat <- as.double(stringr::str_replace(jema$lat, ",", "."))
jema$lon <- as.double(stringr::str_replace(jema$lon, ",", "."))

# remove data before 1990
jema <- jema %>% filter(obs_year >=1990)

jema <- jema %>% rename("datasource"=source) %>% dplyr::select(Species, lat, lon, datasource)
# nrow=18166

# - - - - - - - - - - - - - - - - - - -
## Merge all together ####

data <- dplyr::full_join(sworm, gbif, 
                         by=c("SpeciesBinomial"="species", 
                              "Latitude_decimal_degrees"="decimalLatitude", 
                              "Longitude_decimal_degrees"="decimalLongitude", 
                              "datasource"))

data <- data %>% dplyr::full_join(edapho, 
                 by=c("SpeciesBinomial"="species", 
                      "Latitude_decimal_degrees"="Latitude", 
                      "Longitude_decimal_degrees"="Longitude",
                      "datasource"))

data <- data %>% dplyr::full_join(recon, 
                         by=c("SpeciesBinomial"="Species", 
                              "Latitude_decimal_degrees"="POINT_Y", 
                              "Longitude_decimal_degrees"="POINT_X",
                              "datasource"))

data <- data %>% dplyr::full_join(jema, 
                         by=c("SpeciesBinomial"="Species", 
                              "Latitude_decimal_degrees"="lat", 
                              "Longitude_decimal_degrees"="lon",
                              "datasource"))

data <- tibble::tibble(species=data$SpeciesBinomial, 
                       latitude=data$Latitude_decimal_degrees, 
                       longitude=data$Longitude_decimal_degrees, 
                       datasource=data$datasource)

data #nrow = 101,555
data_merged <- nrow(data)
data <- data[complete.cases(data$longitude, data$latitude, data$species),] #nrow=83,625
data

# remove species with only sp in name
data <- data %>% filter(!stringr::str_detect(data$species, "[[:blank:]]sp"))
# nrow=83,602

data$OBJECTID <- 1:nrow(data) 
data_filter <- nrow(data)

# - - - - - - - - - - - - - - - - - - -
## Save data ####
write.csv(data, file=paste0(here::here(), "/data/Earthworm_occurrence_GBIF-sWorm-Edapho-SoilReCon-JM.csv"),
          row.names = F)

# - - - - - - - - - - - - - - - - - - -
## CoordinateCleaner ####
# flag problems with coordinates
dat_cl <- data.frame(data)
flags <- CoordinateCleaner::clean_coordinates(x = dat_cl, lon = "longitude", lat = "latitude",
                                              species = "species", tests = c("capitals", "centroids", "equal", "gbif", "zeros", "seas"), #normally: test "countries"
                                              country_ref = rnaturalearth::ne_countries("small"), 
                                              country_refcol = "iso_a3")
sum(flags$.summary) #those not flagged! = 80930 out of 83,602

# remove flagged records from the clean data (i.e., only keep non-flagged ones)
dat_cl <- dat_cl[flags$.summary, ]

## Plot flagged records
world.inp <- map_data("world")

pdf(paste0(here::here(), "/figures/CoordinateCleaner_flagged_records_Crassiclitellata.pdf"), width=20)
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
dev.off()

# - - - - - - - - - - - - - - - - - - -
## Some structuring ####
# add species ID
dat_cl <- dplyr::full_join(dat_cl, speciesNames[,c("Species", "SpeciesID")], 
                            by=c("species" = "Species"))

# create x and y column
dat_cl <- dat_cl %>% mutate("x"=longitude, "y"=latitude)

# remove NA in coordinates
dat_cl <- dat_cl[complete.cases(dat_cl$x),]
dat_cl <- dat_cl[complete.cases(dat_cl$y),]

# - - - - - - - - - - - - - - - - - - -
## Save clean data ####
write.csv(dat_cl, file=paste0(here::here(), "/results/Occurrences_", Taxon_name, ".csv"),
          row.names = F)

# load number of records during cleaning process
df_cleaning <- read.csv(file=paste0(here::here(), "/results/NoRecords_cleaning_", Taxon_name, ".csv"))

df_cleaning <- df_cleaning %>% add_row(CleaningStep=c("mergedRecords", "filteredRecords", "coordinateCleaner"), 
                                       NumberRecords=c(data_merged, data_filter, nrow(dat_cl)))

df_cleaning

# save updated number of records during cleaning process
write.csv(df_cleaning, file=paste0(here::here(), "/results/NoRecords_cleaning_", Taxon_name, ".csv"), row.names = F)

