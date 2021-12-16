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

# do so only if you want to download new data
if (checkDownload_rawData==TRUE){
  dat <- rgbif::occ_search(taxonKey = tax_key$data$key, geometry = Search_polygon, hasCoordinate = T, limit=No_records)

  print("Number of occurrences found:"); print(dat$meta$count)
  print("Number of occurrences downloaded"); print(nrow(dat$data)) 
  print("Note: If you have 100,000 records here, you have to manually download the remaining data.")
  print("The limit per search request is set to 100,000.")

  save(dat, file=paste0("RawOccurrences_", Taxon_name, ".RData"))

  rm(dat)
}