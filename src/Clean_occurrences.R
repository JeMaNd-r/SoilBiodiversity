#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#        Clean the occurrence data          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

load(file=paste0(here::here(), "/results/RawOccurrences_", Taxon_name, ".R"))

## Cleaning coordinates based on meta-data ####

# create a table to see how many records get removed.
df.cleaning <- tibble(CleaningStep="RawData", NumberRecords=nrow(dat$data))

# remove records with low coordinate precision (lower than 1km)
#hist(dat$data$coordinateUncertaintyInMeters, breaks = 30)
dat_cl <- dat$data %>% filter(coordinateUncertaintyInMeters <= 1000 | is.na(coordinateUncertaintyInMeters))

df.cleaning <- df.cleaning %>% add_row(CleaningStep="coordinateUncertainty", NumberRecords=nrow(dat_cl))

# remove unsuitable data sources, especially fossils, and keep only those listed here
#table(dat_cl$basisOfRecord) 
dat_cl <- filter(dat_cl, basisOfRecord == "LIVING_SPECIMEN" | basisOfRecord == "HUMAN_OBSERVATION" | 
                   basisOfRecord == "PRESERVED_SPECIMEN" | is.na(basisOfRecord))

df.cleaning <- df.cleaning %>% add_row(CleaningStep="basisOfRecord", NumberRecords=nrow(dat_cl))

# Individual count: remove records with less than 1 observation
#table(dat_cl$individualCount)
dat_cl <- dat_cl %>% filter(individualCount > 0 | is.na(individualCount)) %>% 
  filter(individualCount < 99 | is.na(individualCount))  # high counts are not a problem

df.cleaning <- df.cleaning %>% add_row(CleaningStep="individualCount_min1", NumberRecords=nrow(dat_cl))

# Age of records
print("Occurrence records per year"); print(table(dat_cl$year)); print("Records before 1990 will be removed.")
dat_cl <- dat_cl %>% filter(year >= 1990) 

df.cleaning <- df.cleaning %>% add_row(CleaningStep="year_1990", NumberRecords=nrow(dat_cl))

# taxonomic problems
print("Please indicate if the listed families look good to you.")
print(table(dat_cl$family))  #that looks good

#table(dat_cl$taxonRank)  # We will only include records identified to species level
dat_cl <- dat_cl %>% filter(taxonRank == "SPECIES" | is.na(taxonRank))

df.cleaning <- df.cleaning %>% add_row(CleaningStep="taxonRank_Species", NumberRecords=nrow(dat_cl))

# check how many excluded
print("Records during the cleaning process:"); print(df.cleaning)
print("Records removed [%]:"); print(round((nrow(dat$data) - nrow(dat_cl)) / nrow(dat$data) * 100, 0))

# flag problems
dat_cl <- data.frame(dat_cl)
flags <- CoordinateCleaner::clean_coordinates(x = dat_cl, lon = "decimalLongitude", lat = "decimalLatitude", countries = "countryCode", 
                                              species = "species", tests = c("capitals", "centroids", "equal", "gbif", "zeros", "seas"), #normally: test "countries"
                                              country_ref = rnaturalearth::ne_countries("small"), 
                                              country_refcol = "iso_a3", seas_ref = buffland)
sum(flags$.summary) #those not flagged! = 287841

## Plot flagged records
world.inp <- map_data("world")

print({
  ggplot() + 
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") + 
    xlim(min(dat$data$decimalLongitude, na.rm = T), max(dat$data$decimalLongitude, na.rm = T)) + 
    ylim(min(dat$data$decimalLatitude, na.rm = T), max(dat$data$decimalLatitude, na.rm = T)) + 
    geom_point(data = dat$data, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkred", size = 1, show.legend = T) + 
    geom_point(data = dat_cl, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkgreen", size = 1, show.legend = T) + 
    coord_fixed() + 
    scale_color_manual(name='CoordinateCleaner',
                       values=c('RawRecords = red'='darkred', 'CleanRecords = green'='darkgreen'))+ 
    theme_bw() + theme(axis.title = element_blank())
})
# remove flagged records from the clean data (i.e., only keep non-flagged ones)
dat_cl <- dat_cl[flags$.summary, ]

# filter data columns
dat_cl <- dat_cl %>% select(key, acceptedScientificName, decimalLatitude, decimalLongitude, issues, occurrenceStatus, 
                            speciesKey, order, family, species, genus, specificEpithet, coordinateUncertaintyInMeters,
                            year, month, references, license, geodeticDatum, countryCode, collectionCode, individualCount,
                            samplingProtocol, habitat)

# save only if we want to do so
if (checkSave_cleanData==T){

  # save number of records during cleaning process
  write.csv(df.cleaning, file=paste0(here::here(), "/results/NoRecords_cleaning_", Taxon_name, ".csv"), row.names = F)

  # save cleaned occurrence records
  write.csv(dat_cl, file=paste0(here::here(), "/results/Occurrences_", Taxon_name, ".csv"), row.names = F)
}

# remove temporal R objects
rm(dat, dat_cl, df.cleaning, flags, tax_key, world.inp)

