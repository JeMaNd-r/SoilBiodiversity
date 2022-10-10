#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare results for Faktencheck      #
#          author: Romy Zeiss               #
#            date: 2022-10-10               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(raster)

library(ggplot2) 

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0("./results/Species_list_", Taxon_name, ".csv"))
speciesSub <- speciesNames %>% filter(NumCells_2km >=100) %>% dplyr::select(SpeciesID) %>% unique() %>% c()
#speciesSub <- speciesNames %>% filter(family == "Lumbricidae" & NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique()
speciesSub <- c(speciesSub$SpeciesID)


## Load grid system
grid10 <- raster::raster("data/grid_10k_Germany.tif")
grid10

grid10_info <- read.csv("data/Faktencheck_MTB.csv")
grid10_info <- grid10_info %>% rename("ID"="ï..ID")
head(grid10_info)

#- - - - - - - - - - - - - - - - - - - - -
## Transform data into the 10 km grid ####
#- - - - - - - - - - - - - - - - - - - - -

# load uncertainty extent for all maps
load(file=paste0(here::here(), "/results/_Maps/SDM_Uncertainty_extent_", Taxon_name, ".RData")) #extent_df

# load species distributions
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
head(species_stack)

earthworm_stack <- extent_df %>% inner_join(species_stack)

# transform richness to 10km grid system
earthworm_stack_rich <- raster::rasterize(earthworm_stack %>% dplyr::select(x,y), grid10, earthworm_stack$Richness, fun=mean)
earthworm_stack_rich
raster::plot(earthworm_stack_rich)

# add grid ID
earthworm_stack_rich <- raster::stack(grid10, earthworm_stack_rich)
earthworm_stack_rich

names(earthworm_stack_rich) <- c("ID", "Artensumme")
earthworm_stack_rich

# transform raster to data frame
earthworm_rich <- as.data.frame(rasterToPoints(earthworm_stack_rich))
head(earthworm_rich)

# remove grid cells outside of grid (without ID)
earthworm_rich <- earthworm_rich %>% filter(!is.na(ID))
head(earthworm_rich)

# add grid cell information
earthworm_rich <- grid10_info %>% full_join(earthworm_rich %>% dplyr::select(ID, Artensumme))
head(earthworm_rich)

#- - - - - - - - - - - - - - - - - - - - -
## Explore data ####
#- - - - - - - - - - - - - - - - - - - - -

# plot richness
world.inp <- map_data("world")
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(5, 17) +
  ylim(46, 56) +
  
  geom_point(data=earthworm_rich %>% filter(Artensumme>0), 
            aes(x=X, y=Y, col=Artensumme), shape=15)+
  ggtitle("Species richness (number of species)")+
  scale_color_viridis_c()+
  geom_point(data=earthworm_rich %>% filter(Artensumme==0), aes(x=X, y=Y), col="grey60", shape=15)+
  geom_point(data=earthworm_rich %>% filter(is.na(Artensumme)), aes(x=X, y=Y), col="white", shape=15)+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(), legend.text = element_text(size=10),
        legend.position = c(0.1,0.9), legend.direction = "horizontal")

# amphibian data
amphib_rich <- read.delim("data/Faktencheck_Rept_Amph_MTB_lat_long.txt", sep=";")
head(amphib_rich)

ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(5, 17) +
  ylim(46, 56) +
  
  geom_point(data=amphib_rich %>% filter(Artensumme>0), 
             aes(x=X, y=Y, col=Artensumme), shape=15)+
  ggtitle("Species richness (number of species)")+
  scale_color_viridis_c()+
  geom_point(data=amphib_rich %>% filter(Artensumme==0), aes(x=X, y=Y), col="grey60", shape=15)+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(), legend.text = element_text(size=10),
        legend.position = c(0.1,0.9), legend.direction = "horizontal")


#- - - - - - - - - - - - - - - - - - - - -
## Save richness at 10km in Germany ####
#- - - - - - - - - - - - - - - - - - - - -

write.table(earthworm_rich, "results/Faktencheck_Regenwurm_Artensumme_RZ_202210.txt",
            row.names = F)


# save data at 5km²
write.table(earthworm_stack %>% dplyr::select(x,y,Richness), 
            "results/Faktencheck_Regenwurm_Richness_5km_RZ_202210.txt",
            row.names = F)


# # - - - - - - - - - - - - - - - - - - - - -
# # # MAYBE: Transform species data into the 10 km grid ####
# # - - - - - - - - - - - - - - - - - - - - -
# 
# # transform species occurrence to grid
# earthworm_stack_species <- raster::rasterize(earthworm_stack %>% dplyr::select(x,y), grid10, earthworm_stack %>% dplyr::select(Allol_chlo_current:Satch_mamm_current), fun=max)
# earthworm_stack_species <- raster::mask(earthworm_stack_species, grid10)
# earthworm_stack_species
# raster::plot(earthworm_stack_species[[1]])
# 
# # transform raster to data frame
# earthworm_species <- as.data.frame(rasterToPoints(earthworm_stack_species))
# head(earthworm_species)
# 
# # rename columns
# colnames(earthworm_species) <- substr(colnames(earthworm_species), 1, 10)
# 
# #- - - - - - - - - - - - - - - - - - - - -
# ## Add true occurrences with their sampling date ####
# #- - - - - - - - - - - - - - - - - - - - -
# 
# # load data with sampling year information
# occ_points <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_10km_", Taxon_name, ".csv"))
# str(occ_points)
# 
# # add year information
# earthworm_species <- earthworm_species %>% full_join(occ_points)


