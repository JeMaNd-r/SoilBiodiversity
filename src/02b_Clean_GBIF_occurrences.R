#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#        Clean the occurrence data          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

## Load data ####

# use data downloaded via R script
#load(file=paste0(here::here(), "/results/RawOccurrences_", Taxon_name, ".RData"))
#dat <- dat$data

# OR use data with DOI downloaded directly from GBIF
dat <- readr::read_delim(file=paste0(here::here(), "/data/GBIF_", Taxon_name, "/occurrence.txt"))

#- - - - - - - - - - - - - - - - - - - - - -
## Cleaning coordinates based on meta-data ####

# create a table to see how many records get removed.
df_cleaning <- tibble::tibble(CleaningStep="GBIF_RawData", NumberRecords=nrow(dat))

# remove records outside of Europe
dat_cl <- dat %>% filter(extent_Europe[1] <= decimalLongitude &  decimalLongitude <= extent_Europe[2]) %>% 
  filter(extent_Europe[3] <= decimalLatitude &  decimalLatitude <= extent_Europe[4])

df_cleaning <- df_cleaning %>% add_row(CleaningStep="GBIF_Europe", NumberRecords=nrow(dat_cl))

# remove records from iNaturalist
dat_cl <- dat_cl %>% filter(publisher!="iNaturalist.org")

df_cleaning <- df_cleaning %>% add_row(CleaningStep="GBIF_iNat", NumberRecords=nrow(dat_cl))

# remove records with low coordinate precision (lower than 1km)
#hist(dat$coordinateUncertaintyInMeters, breaks = 30)
dat_cl <- dat_cl %>% filter(coordinateUncertaintyInMeters <= 1000 | is.na(coordinateUncertaintyInMeters))

df_cleaning <- df_cleaning %>% add_row(CleaningStep="GBIF_coordinateUncertainty", NumberRecords=nrow(dat_cl))

# remove records with issues in species identification or lat/long
unique(unlist(str_split(unique(dat_cl$issue), ";")))

dat_cl <- dat_cl %>% filter(!stringr::str_detect(issue, "RECORDED_DATE_INVALID"))
df_cleaning <- df_cleaning %>% add_row(CleaningStep="GBIF_issueDate", NumberRecords=nrow(dat_cl))

dat_cl <- dat_cl %>% filter(!stringr::str_detect(issue, "GEODETIC_DATUM_INVALID"))
df_cleaning <- df_cleaning %>% add_row(CleaningStep="GBIF_issueGeoDate", NumberRecords=nrow(dat_cl))

# note: coordinate_rounded are rounded to 1m precision

# remove unsuitable data sources, especially fossils, and keep only those listed here
#table(dat_cl$basisOfRecord) 
dat_cl <- filter(dat_cl, basisOfRecord == "LIVING_SPECIMEN" | basisOfRecord == "HUMAN_OBSERVATION" | 
                   basisOfRecord == "PRESERVED_SPECIMEN" | is.na(basisOfRecord))

df_cleaning <- df_cleaning %>% add_row(CleaningStep="GBIF_basisOfRecord", NumberRecords=nrow(dat_cl))

# Individual count: remove records with less than 1 observation
#table(dat_cl$individualCount)
dat_cl <- dat_cl %>% filter(individualCount > 0 | is.na(individualCount)) %>% 
  filter(individualCount < 99 | is.na(individualCount))  # high counts are not a problem

df_cleaning <- df_cleaning %>% add_row(CleaningStep="GBIF_individualCount_min1", NumberRecords=nrow(dat_cl))

# Age of records
print("Occurrence records per year"); print(table(dat_cl$year)); print("Records before 1970 will be removed.")
dat_cl <- dat_cl %>% filter(year >= 1970) 

df_cleaning <- df_cleaning %>% add_row(CleaningStep="GBIF_year1970", NumberRecords=nrow(dat_cl))

# taxonomic problems
print("Please check if the listed families look good to you.")
print(table(dat_cl$family))  #that looks good

#table(dat_cl$taxonRank)  # We will only include records identified to species level
dat_cl <- dat_cl %>% filter(taxonRank == "SPECIES" | is.na(taxonRank))

df_cleaning <- df_cleaning %>% add_row(CleaningStep="GBIF_taxonRank_Species", NumberRecords=nrow(dat_cl))

# check how many excluded
print("Records during the cleaning process:"); print(df_cleaning)
print("Records removed [%]:"); print(round((nrow(dat) - nrow(dat_cl)) / nrow(dat) * 100, 0))

print("Records removed from Europe [%]:")
print(round((as.numeric(df_cleaning[2,2]) - nrow(dat_cl)) / as.numeric(df_cleaning[2,2]) * 100, 0))

## will be done in a later step 
# # flag problems
# dat_cl <- data.frame(dat_cl)
# flags <- CoordinateCleaner::clean_coordinates(x = dat_cl, lon = "decimalLongitude", lat = "decimalLatitude", countries = "countryCode", 
#                                               species = "species", tests = c("capitals", "centroids", "equal", "gbif", "zeros", "seas"), #normally: test "countries"
#                                               country_ref = rnaturalearth::ne_countries("small"), 
#                                               country_refcol = "iso_a3", seas_ref = buffland)
# sum(flags$.summary) #those not flagged! = 287841
# 
# # remove flagged records from the clean data (i.e., only keep non-flagged ones)
# dat_cl <- dat_cl[flags$.summary, ]

## Plot flagged records
world.inp <- map_data("world")

pdf(paste0(here::here(), "/figures/CoordinateCleaner_flagged_records_Crassiclitellata_1970.pdf"), width=18, height=6)
  ggplot() +
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(min(dat$decimalLongitude, na.rm = T), max(dat$decimalLongitude, na.rm = T)) +
    ylim(min(dat$decimalLatitude, na.rm = T), max(dat$decimalLatitude, na.rm = T)) +
    geom_point(data = dat, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkred", size = 1, show.legend = T) +
    geom_point(data = dat_cl, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkgreen", size = 1, show.legend = T) +
    coord_fixed() +
    scale_color_manual(name='ManualCleaning',
                       values=c('RawRecords = red'='darkred', 'CleanRecords = green'='darkgreen'))+
    theme_bw() + theme(axis.title = element_blank())
dev.off()


# filter data columns
dat_cl <- dat_cl %>% dplyr::select(acceptedScientificName, decimalLatitude, decimalLongitude, issue, occurrenceStatus, 
                            speciesKey, order, family, species, genus, specificEpithet, coordinateUncertaintyInMeters,
                            year, month, references, license, verbatimCoordinateSystem, countryCode, collectionCode, individualCount,
                            samplingProtocol, habitat)



# remove commata
dat_cl$issue <- gsub(","," ",dat_cl$issue)
dat_cl$acceptedScientificName <- gsub(","," ",dat_cl$acceptedScientificName)

# transform into wide format
dat_wide <- dat_cl
dat_wide$presence <- 1
dat_wide <- pivot_wider(dat_wide, id_cols=c(decimalLongitude, decimalLatitude), 
                            names_from = species, values_from = presence,
                            values_fn = max)
dat_wide <- as.data.frame(dat_wide)

# save only if we want to do so
#if (checkSave_cleanData==T){

  # save number of records during cleaning process
  write.csv(df_cleaning, file=paste0(here::here(), "/results/NoRecords_cleaning_", Taxon_name, ".csv"), row.names = F)

  # save cleaned occurrence records
  write.csv(dat_cl, file=paste0(here::here(), "/results/Occurrences_GBIF_", Taxon_name, ".csv"), row.names = F)

  # save cleaned occurrence records in wide format
  write.csv(dat_wide, file=paste0(here::here(), "/results/Occurrences_GBIF_wide_", Taxon_name, ".csv"), row.names = F)
  
#}

# remove temporal R objects
rm(dat, dat_cl, df_cleaning, flags, tax_key, world.inp)

