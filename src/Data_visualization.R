#- - - - - - - - - - - - - - - - -#
#   General data visualization    #
#                                 #
#     author: Romy Zeiss          #
#       date: 18.08.2021          #
#- - - - - - - - - - - - - - - - -#

library(here)  #instead of setwd()
library(tidyverse)
library(rworldmap)
library(sp)
library(rgdal)
#library(maps)       # Provides functions that let us plot the maps
#library(mapdata)    # Contains the hi-resolution points that mark out the countries.
#library(ggmap)

## Load occurrence data ####
data <- readr::read_csv(here::here("data", "Earthworm_occurrence_Germany.txt"))

## Plot point occurrences ####
ggplot(data=data, aes(x="longitude", y="latitude"))+
  geom_point()

# newmap <- getMap(resolution = "low")
# setMapExtents(mapRegion = "Germany")
# plot(newmap, xlim = c(5.9, 15), ylim = c(47.3, 54.9), asp = 1)
# points(data$longitude, data$latitude, col=2, pch=18, cex=0.1)

# colors representing the data sources
plot(newmap, xlim = c(5.9, 15), ylim = c(47.3, 54.9), asp = 1)
for(i in 1:length(unique(data$datasource))){
  temp.data <- data[data$datasource==unique(data$datasource)[i],]
  points(temp.data$longitude, temp.data$latitude, col=i, pch=18, cex=0.1)
}

mi_counties <- ggplot2::map_data("world", "Germany") %>% 
  select(lon = long, lat, group)
head(mi_counties)

ggplot(mi_counties, aes(lon, lat)) +
  geom_polygon(fill = "white", colour = "grey50") +
  geom_point(data = data, mapping = aes(x = longitude, y = latitude), colour = "red")+
  coord_quickmap()

# plot cropped to Germany
coordinates(data) <- c("longitude","latitude")
proj4string(data) <- CRS("+proj=longlat")  

shp.german <- readOGR(dsn = here::here("data/Shapefile_Germany"), layer = "GISPORTAL_GISOWNER01_GERMANY_INTBORDER_15") #--> read in the shapefile one is using 

plot(shp.german)
points(data, pch = 20,cex = 0.5, col= "#B000B5")


