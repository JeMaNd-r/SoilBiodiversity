#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Extract uncertainty               #
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
speciesSub <- speciesNames %>% filter(NumCells_2km_biomod >=100) %>% dplyr::select(SpeciesID) %>% unique() %>% c()
#speciesSub <- speciesNames %>% filter(family == "Lumbricidae" & NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique()
speciesSub <- c(speciesSub$SpeciesID)

# load environmental space as data frame
load(paste0(here::here(),"/results/EnvPredictor_5km_df_clipped.RData")) #Env_clip_df

uncertain_df <- Env_clip_df %>% dplyr::select(x, y)

for(spID in speciesSub){try({
  
  print(paste0("Species: ", spID))

  # list files in species-specific BIOMOD folder
  temp_files <- list.files(paste0(here::here(), "/results/biomod_files/", stringr::str_replace(spID, "_", "."), "/proj_modeling"), full.names = TRUE)
          
  #setwd(paste0(here::here(), "/results/", Taxon_name))

  myBiomodEnProj <- get(load(temp_files[stringr::str_detect(temp_files,"ensemble.RData")]))
       
  # save predictions as raster file
  temp_prediction <- myBiomodEnProj[,1] #column with CoV
  temp_prediction <- as.numeric(temp_prediction)
  # add names of grid cell (only for those that have no NA in any layer)
  names(temp_prediction) <- rownames(Env_clip_df)
  temp_prediction <- as.data.frame(temp_prediction)
  temp_prediction$x <- Env_clip_df$x
  temp_prediction$y <- Env_clip_df$y
  temp_prediction <- temp_prediction %>% full_join(Env_clip_df %>% dplyr::select(x,y)) %>%
     rename("layer" = temp_prediction)
  temp_prediction$layer <- temp_prediction$layer / 1000
  temp_prediction[,spID] <- temp_prediction$layer
          
  # add layer to stack
  uncertain_df <- full_join(uncertain_df, temp_prediction %>% dplyr::select(x,y, spID))
 
})}

uncertain_df$Mean <- rowMeans(uncertain_df %>% dplyr::select(-x, -y), na.rm=T)

# calculate sd of predictions
uncertain_df$SD <- apply(uncertain_df %>% dplyr::select(-x, -y, -Mean), 1, sd, na.rm = TRUE)

head(uncertain_df)

# save species' uncertainty map
save(uncertain_df, file=paste0(here::here(), "/results/_Maps/SDM_Uncertainty_", Taxon_name, ".RData"))
load(file=paste0(here::here(), "/results/_Maps/SDM_Uncertainty_", Taxon_name, ".RData")) #uncertain_df


# extract area with uncertainty lower than threshold
summary(uncertain_df$Mean)

extent_df <- uncertain_df %>% filter(Mean<0.1 & !is.na(Mean)) %>% dplyr::select(x,y)
save(extent_df, file=paste0(here::here(), "/results/_Maps/SDM_Uncertainty_extent_", Taxon_name, ".RData"))


#- - - - - - - - - - - - - - - - - - - - -
## Clamping mask (extrapolated areas outside of calibration data range) ####
#- - - - - - - - - - - - - - - - - - - - -

clamping_df <- Env_clip_df %>% dplyr::select(x, y)

for(spID in speciesSub){try({
  
  print(paste0("Species: ", spID))
  
  # list files in species-specific BIOMOD folder
  temp_files <- list.files(paste0(here::here(), "/results/biomod_files/", stringr::str_replace(spID, "_", "."), "/proj_modeling"), full.names = TRUE)
  
  myBiomodEnProj <- get(load(temp_files[stringr::str_detect(temp_files,"ClampingMask.RData")]))
  myBiomodEnProj <- data.frame(myBiomodEnProj)
  colnames(myBiomodEnProj) <- spID
  
  # add layer to stack
  clamping_df <- cbind(clamping_df, myBiomodEnProj)
  
})}

head(clamping_df) #number of variables outside the range of calibration data

clamping_df$Mean <- rowMeans(clamping_df %>% dplyr::select(-x, -y), na.rm=T)
clamping_df$Sum <- rowSums(clamping_df %>% dplyr::select(-x, -y), na.rm=T)

# calculate sd of predictions
clamping_df$SD <- apply(clamping_df %>% dplyr::select(-x, -y, -Mean, -Sum), 1, sd, na.rm = TRUE)
clamping_df$Max <- apply(clamping_df %>% dplyr::select(-x, -y, -Mean, -Sum, -SD), 1, max, na.rm = TRUE)

head(clamping_df)

save(clamping_df, file=paste0(here::here(), "/results/_Maps/SDM_ClampingMask_current_", Taxon_name, ".RData"))

#- - - - - - - - - - - - - - - - - - - - -
# same for future projections
clamping_df <- Env_clip_df %>% dplyr::select(x, y)

for(spID in speciesSub){try({
  
  print(paste0("Species: ", spID))
  
  for(no_future in scenarioNames){
    
    # list files in species-specific BIOMOD folder
    temp_files <- list.files(paste0(here::here(), "/results/biomod_files/", stringr::str_replace(spID, "_", "."), "/proj_modeling_future_", no_future, "_TP"), full.names = TRUE)
    
    myBiomodEnProj <- get(load(temp_files[stringr::str_detect(temp_files,"ClampingMask.RData")]))
    myBiomodEnProj <- data.frame(myBiomodEnProj)
    colnames(myBiomodEnProj) <- paste0(spID, "_", no_future)
    
    # add layer to stack
    clamping_df <- cbind(clamping_df, myBiomodEnProj)
    
  }
})}

head(clamping_df) #number of variables outside the range of calibration data

clamping_df$Mean <- rowMeans(clamping_df %>% dplyr::select(-x, -y), na.rm=T)
clamping_df$Sum <- rowSums(clamping_df %>% dplyr::select(-x, -y), na.rm=T)

# calculate sd of predictions
clamping_df$SD <- apply(clamping_df %>% dplyr::select(-x, -y, -Mean, -Sum), 1, sd, na.rm = TRUE)
clamping_df$Max <- apply(clamping_df %>% dplyr::select(-x, -y, -Mean, -Sum, -SD), 1, max, na.rm = TRUE)

head(clamping_df)

save(clamping_df, file=paste0(here::here(), "/results/_Maps/SDM_ClampingMask_future_", Taxon_name, ".RData"))

