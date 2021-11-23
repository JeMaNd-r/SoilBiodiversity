#- - - - - - - - - - - - - - - - -#
#  Data preparation for ArcGIS    #
#                                 #
#     author: Romy Zeiss          #
#       date: 22.11.2021          #
#- - - - - - - - - - - - - - - - -#

#library(here)  #instead of setwd()
library(tidyverse)
library(rgbif)
library(sp)
library(countrycode)
library(CoordinateCleaner)

setwd("I:/eie/==PERSONAL/RZ SoilBON/SoilBiodiversity")

#- - - - - - - - - - - - - - - - -
## Download data from GBIF ####

# Use the name_suggest function to get the gbif taxon key
tax_key <- name_suggest(q = "Lumbricidae", rank = "Order")

# Sometimes groups have multiple taxon keys, in this case three, so we will
# check how many records are available for them
lapply(tax_key$key, "occ_count")

# Here the first one is relevant, check for your group!
tax_key <- tax_key$key[1]

# count occurrences
occ_count(tax_key, country = "DE")

# define a study region
study_a <- "POLYGON((-35 -4.5, -38.5 -4.5, -38.5 -7, -35 -7, -35 -4.5))"

# define region using a shapefile
amz <- readOGR('inst', layer = 'Amazonia')
#rgeos::writeWKT(amz)
# Or, best use the extent of the shape, since it is simple:
ex <- raster::extent(amz)
ex <- as(ex, "SpatialPolygons")
ex <- rgeos::writeWKT(ex)

# get the occurrences
dat_ne <- occ_search(taxonKey = tax_key, return = "data", hasCoordinate = T, 
                     geometry = study_a, limit = 1000)

## Data cleaning ####
#... see script course Macroecology
clean_coordinates(x = dat,
                  lon = "decimalLongitude",
                  lat = "decimalLatitude",
                  countries = "countryCode",
                  species = "species",
                  tests = c("capitals", "centroids", "equal","gbif", "institutions",
                            "zeros", "countries"))

## Vizualize difference between raw and cleaned data ####
world.inp <- map_data("world")

ggplot() + geom_map(data = world.inp, map = world.inp, 
                    aes(x = long, y = lat, map_id = region), 
                    fill = "grey80") + 
  xlim(min(dat$decimalLongitude, na.rm = T),
       max(dat$decimalLongitude, na.rm = T)) + 
  ylim(min(dat$decimalLatitude, na.rm = T), 
       max(dat$decimalLatitude, na.rm = T)) + 
  geom_point(data = dat, aes(x = decimalLongitude,
                             y = decimalLatitude), 
             colour = "darkred", size = 1) + 
  geom_point(data = dat_cl, aes(x = decimalLongitude, 
                                y = decimalLatitude), 
             colour = "darkgreen", size = 1) + 
  coord_fixed() + theme_bw() + theme(axis.title = element_blank())


ggplot() + geom_map(data = world.inp, map = world.inp, 
                    aes(x = long, y = lat, map_id = region), 
                    fill = "grey80") + 
  xlim(min(dat$decimalLongitude, na.rm = T), 
       max(dat$decimalLongitude, na.rm = T)) + 
  ylim(min(dat$decimalLatitude, na.rm = T), 
       max(dat$decimalLatitude, na.rm = T)) + 
  geom_point(data = dat_cl, aes(x = decimalLongitude,
                                y = decimalLatitude, 
                                colour = dataset), size = 1) + 
  coord_fixed() + theme_bw() + 
  theme(axis.title = element_blank())


