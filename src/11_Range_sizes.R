#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#    Calculate range sizes of species       #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(raster)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0("./results/Species_list_", Taxon_name, ".csv"))
speciesSub <- speciesNames %>% filter(NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique() %>% c()
#speciesSub <- speciesNames %>% filter(family == "Lumbricidae" & NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique()
speciesSub <- c(speciesSub$SpeciesID)

# covariates in order of importance (top 10 important)
covarsNames <- c("MAT", "MAP_Seas", "Dist_Coast", "Agriculture", "pH", 
                 "P", "CEC", "Elev", "Clay.Silt", "Pop_Dens")

# define future scenarios
scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
                               "mri-esm2-0", "ukesm1-0-ll"), "_",
                             rep(c("ssp126", "ssp370", "ssp585"),5)))

# load future species stack
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_future_species_", Taxon_name, ".RData")) #future_stack

## Calculate species ranges (area) ####
range_df <- data.frame("scenario"=colnames(future_stack), 
                       "cells"=colSums(future_stack, na.rm=T), 
                       "area_km2"=colSums(future_stack, na.rm=T)*5) %>%
  filter(scenario!="x" & scenario!="y") %>%
  tidyr::separate(scenario, c("SpeciesID", "scenario"), "[.]")
rownames(range_df) <- NULL
head(range_df)

range_sum <- range_df[range_df$SpeciesID %in% range_df$SpeciesID[stringr::str_detect(range_df$SpeciesID, "_future")],] %>% group_by(SpeciesID) %>% dplyr::select(-scenario) %>% 
  summarise(across(everything(), list(min = min, max = max, mean = mean, sd = sd)))

# add current ranges
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
range_df_current <- data.frame("scenario"=colnames(species_stack), 
                               "cells"=colSums(species_stack, na.rm=T), 
                               "area_km2"=colSums(species_stack, na.rm=T)*5) %>%
  filter(scenario!="x" & scenario!="y")
rownames(range_df_current) <- NULL
head(range_df_current)  

range_sum <- range_sum %>% cbind(range_df_current %>% filter(scenario!="Richness"))
range_sum

range_sum$area_km2_change <- range_sum$area_km2_mean - range_sum$area_km2
range_sum$area_km2_change_p <- range_sum$area_km2_change / range_sum$area_km2
range_sum$area_km2_p_sd <- range_sum$area_km2_sd / range_sum$area_km2

range_sum %>% arrange(area_km2_change_p)
mean(range_sum$area_km2_change_p); sd(range_sum$area_km2_change_p)

write.csv(range_sum, file=paste0(here::here(), "/results/Range_shift_", Taxon_name, ".csv"), row.names=F)
