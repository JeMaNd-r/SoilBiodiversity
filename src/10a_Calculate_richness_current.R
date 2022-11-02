#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Calculate species richness under     #
#             current climate               #
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

# load environmental data 5km
load(paste0(here::here(),"/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df

#- - - - - - - - - - - - - - - - - - - - - -
## Create maps and calculate richness ####
#- - - - - - - - - - - - - - - - - - - - - -
# create empty data frame
species_stack <- Env_norm_df %>% dplyr::select(x, y)

# for loop through all species
for(spID in speciesSub){ try({
  
  ## Load probability maps 
  load(file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_biomod_", spID, ".RData")) #biomod_list
  best_pred <- biomod_list$prediction
  
  print(paste0(spID, " successfully loaded."))
  
  ## Transform to binary maps ####
  
  # extract threshold to define presence/absence: TSS [row 2]
  temp_thresh <- biomod_list$validation[2,str_detect(colnames(biomod_list$validation), "EMcaByTSS_mergedAlgo_mergedRun_mergedData.Cutoff")]/1000
  if(is.na(temp_thresh)) temp_tresh <- 0.9
  
  # change to binary
  best_pred[best_pred$layer>=temp_thresh & !is.na(best_pred$layer), "layer"] <- 1
  best_pred[best_pred$layer<temp_thresh & !is.na(best_pred$layer), "layer"] <- 0
  
  best_pred[,paste0(spID,"_current")] <- best_pred$layer
  best_pred <- best_pred[,c("x","y",paste0(spID,"_current"))]
  
  # save binary
  save(best_pred, file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_bestPrediction_binary_", Taxon_name, "_", spID, ".RData"))
  
  print(paste0("Saved binary prediction of ", spID))
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Stack species binary maps ####
  
  # add species dataframe to stacked dataframe
  species_stack <- species_stack %>% full_join(best_pred, by=c("x","y"))
  
  print(paste0("Added binary prediction of ", spID, " to the species stack"))
  
  rm(temp_thresh, best_pred)
}, silent=T)}  

head(species_stack)


#- - - - - - - - - - - - - - - - - - - - - -
## Calculate richness ####
species_stack$Richness <- rowSums(species_stack %>% dplyr::select(-x, -y), na.rm=F)


#- - - - - - - - - - - - - - - - - - - - - -
## Save species stack ####
save(species_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData"))

#load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack

#- - - - - - - - - - - - - - - - - - - - - -
## Beta diversity ####
# vegan::vegdist() cannot handle such big vectors
# alternatively, we calculate clusters
temp_matrix <- species_stack[species_stack$Richness>0 & !is.na(species_stack$Richness),colnames(species_stack)[stringr::str_detect("_current", string = colnames(species_stack))]]
temp_matrix <- temp_matrix[complete.cases(temp_matrix),]

temp_clust <- kmeans(temp_matrix, 10) 
temp_clust <- cbind(species_stack[species_stack$Richness>0 & !is.na(species_stack$Richness),],
                    "Cluster"=temp_clust$cluster)

species_stack <- species_stack %>% full_join(temp_clust)


## View clusters
world.inp <- map_data("world")

png(file=paste0(here::here(), "/figures/Clusters_K10_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(-23, 60) +
  ylim(31, 75) +
  
  geom_tile(data=species_stack %>% filter(Richness>0), aes(x=x, y=y, fill=as.factor(Cluster)))+
  ggtitle("Clusters (k=10)")+
  scale_fill_viridis_d()+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()


