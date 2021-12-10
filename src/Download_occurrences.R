#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#        Download occurrence data           #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#taxize::get_ids("Lumbricus terrestris")

# get style of taxon name(s) in GBIF database
tax_key <- rgbif::name_suggest(q=Taxon_name, rank=Taxon_rank)
print(tax_key); print("Is this correct?")

# # define a list with the European country codes
# country_code <- c('AL', 'AD', 'AM', 'AT', 'BY', 'BE', 'BA', 'BG', 'CH', 'CY', 'CZ', 'DE',
#                   'DK', 'EE', 'ES', 'FO', 'FI', 'FR', 'GB', 'GE', 'GI', 'GR', 'HU', 'HR',
#                   'IE', 'IS', 'IT', 'LI', 'LT', 'LU', 'LV', 'MC', 'MK', 'MT', 'NO', 'NL', 'PL',
#                   'PT', 'RO', 'RS', 'RU', 'SE', 'SI', 'SK', 'SM', 'TR', 'UA', 'VA')


# count occurrences and save in result table
# we only want occurrence records that have coordinates (georeferenced), were collected after 1990 (from) 
# and in one of the European countries (country)
dat <- rgbif::occ_search(taxonKey = tax_key$data$key, geometry = search_polygon, hasCoordinate = T, limit=1000)

print("Number of occurrences found:"); print(dat$meta$count)
print("Number of occurrences downloaded"); print(nrow(dat$data)) 
print("Note: If you have 100,000 records here, you have to manually download the remaining data.")
print("The limit per search request is set to 100,000.")


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


# save cleaned occurrence records
#write.csv(data_cl, file=paste0("Occurrences_", Taxon_name, ""))

