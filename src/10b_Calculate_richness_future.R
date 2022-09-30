#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Calculate species richness under     #
#             future climate                #
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
covarsNames <- c("MAT", "Dist_Coast", "MAP_Seas", "CEC", "Elev",
                 "P", "Pop_Dens", "Agriculture", "pH", "Clay.Silt")

# define future scenarios
scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
                               "mri-esm2-0", "ukesm1-0-ll"), "_",
                             rep(c("ssp126", "ssp370", "ssp585"),5)))

# load environmental data 5km
load(paste0(here::here(),"/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df


#- - - - - - - - - - - - - - - - - - - - - -
## Create future maps and calculate richness ####
for(no_future in scenarioNames){
  
  # one loop per future climate subset, one with both future, each one with only 1 future and 1 current climate
  for(subclim in c("TP", "T", "P")){
    
    #- - - - - - - - - - - - - - - - - - - - - -
    # create empty data frame
    species_stack <- Env_norm_df %>% dplyr::select(x, y)
    
    # for loop through all species
    for(spID in speciesSub){ try({
      
      ## Load probability maps 
      load(file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_2041-2070_", no_future, "_", subclim, "_biomod_", spID,  ".RData")) #biomod_list
      best_pred <- temp_prediction
      
      # load model information 
      load(file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_biomod_", spID, ".RData")) #biomod_list
      
      print(paste0(spID, " successfully loaded."))
      
      ## Transform to binary maps ####
      
      # extract threshold to define presence/absence
      temp_thresh <- biomod_list$validation[2,str_detect(colnames(biomod_list$validation), "EMcaByTSS_mergedAlgo_mergedRun_mergedData.Cutoff")]/1000
      if(is.na(temp_thresh)) temp_tresh <- 0.9
      
      # change to binary
      best_pred[best_pred$layer>=temp_thresh & !is.na(best_pred$layer), "layer"] <- 1
      best_pred[best_pred$layer<temp_thresh & !is.na(best_pred$layer), "layer"] <- 0
      
      best_pred[,paste0(spID,"_future")] <- best_pred$layer
      best_pred <- best_pred[,c("x","y",paste0(spID,"_future"))]
      
      # save binary
      save(best_pred, file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_bestPrediction_binary_2041-2070_", no_future, "_", subclim, "_biomod_", spID,  ".RData"))
      
      print(paste0("Saved binary prediction of ", spID))
      
      #- - - - - - - - - - - - - - - - - - - - - -
      ## Stack species binary maps ####
      
      # add species dataframe to stacked dataframe
      species_stack <- species_stack %>% full_join(best_pred, by=c("x","y"))
      
      print(paste0("Added binary prediction of ", spID, " to the species stack"))
      
      rm(temp_thresh, best_pred, temp_prediction)
    }, silent=T)}  
    
    head(species_stack)
    
    
    #- - - - - - - - - - - - - - - - - - - - - -
    ## Calculate richness ####
    species_stack$Richness <- rowSums(species_stack %>% dplyr::select(-x, -y), na.rm=F)
    
    #- - - - - - - - - - - - - - - - - - - - - -
    ## Save species stack ####
    save(species_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", "2041-2070_", no_future, "_", subclim, ".RData"))
    
  }
}


#- - - - - - - - - - - - - - - - - - - - -
## Species-specific mean future predictions ####

future_stack <- get(load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", "2041-2070_", scenarioNames[1], "_TP.RData"))) #species_stack

for(no_future in scenarioNames[2:length(scenarioNames)]){
  
  load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", "2041-2070_", no_future, "_TP.RData")) #species_stack
  
  # add layer to stack
  future_stack <- full_join(future_stack, species_stack, suffix=c("", paste0(".", no_future)), by=c("x", "y"))
  
}
colnames(future_stack)[3:22] <- paste0(colnames(future_stack)[3:22], ".", scenarioNames[1])
colnames(future_stack)

# calculate average future prediction per species
for(spID in unique(speciesNames[speciesNames$NumCells_2km>=100,]$SpeciesID)){ try({
  future_stack[,as.character(paste0(spID, ".future_mean"))] <- rowMeans(future_stack[,stringr::str_detect(colnames(future_stack), spID)], na.rm=T)
  future_stack[,as.character(paste0(spID, ".future_max"))] <- matrixStats::rowMaxs(as.matrix(future_stack[,stringr::str_detect(colnames(future_stack), spID)]), na.rm=T)
  future_stack[,as.character(paste0(spID, ".future_min"))] <- matrixStats::rowMins(as.matrix(future_stack[,stringr::str_detect(colnames(future_stack), spID)]), na.rm=T)
})}
colnames(future_stack)

save(future_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_future_species_", Taxon_name, ".RData"))


#- - - - - - - - - - - - - - - - - - - - -
## Average future predictions ####
world.inp <- map_data("world")
average_stack <- Env_norm_df %>% dplyr::select(x, y)

for(no_future in scenarioNames){
  
  # only plot subclim scenario TP
  subclim <- "TP"
  
  load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", "2041-2070_", no_future, "_", subclim, ".RData")) #species_stack
  species_stack[,as.character(no_future)] <- species_stack$Richness
  species_stack <- species_stack[,c("x","y",no_future)]
  
  # add layer to stack
  average_stack <- full_join(average_stack, species_stack)
  
}

average_stack$Mean <- rowMeans(average_stack %>% dplyr::select(-x, -y), na.rm=T)

# calculate average per SSP
for(temp_ssp in c("ssp126", "ssp370", "ssp585")){
  temp_cols <- colnames(average_stack)[stringr::str_detect(colnames(average_stack), temp_ssp)]
  average_stack[,as.character(paste0(temp_ssp, "_mean"))] <- rowMeans(average_stack[,temp_cols])
  average_stack[,as.character(paste0(temp_ssp, "_max"))] <- matrixStats::rowMaxs(as.matrix(average_stack[,temp_cols]))
  average_stack[,as.character(paste0(temp_ssp, "_min"))] <- matrixStats::rowMins(as.matrix(average_stack[,temp_cols]))
  average_stack[,as.character(paste0(temp_ssp, "_sd"))] <- matrixStats::rowSds(as.matrix(average_stack[,temp_cols]))
}
colnames(average_stack)

# Calculate percent change in distribution
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
average_stack <- average_stack %>% dplyr::rename(FutureRichness=Mean) %>%
  full_join(species_stack %>% dplyr::select(x,y,Richness))

average_stack$Change <- (average_stack$FutureRichness - average_stack$Richness)
average_stack$Change_f <- cut(average_stack$Change, 
                              breaks=c(-15, -10, -5, 0, 5, 10),
                              labels=c("[-15,-10]", "[-10,-5]", "[-5,0]", "[0,5]", "[5,10]"))

average_stack$Change_ssp126 <- (average_stack$ssp126_mean - average_stack$Richness)
average_stack$Change_f_ssp126 <- cut(average_stack$Change_ssp126, 
                                     breaks=c(-15, -10, -5, 0, 5, 10, 15),
                                     labels=c("[-15,-10]", "[-10,-5]", "[-5,0]", "[0,5]", "[5,10]", "[10,15]"))
average_stack$Change_ssp370 <- (average_stack$ssp370_mean - average_stack$Richness)
average_stack$Change_f_ssp370 <- cut(average_stack$Change_ssp370, 
                                     breaks=c(-15, -10, -5, 0, 5, 10, 15),
                                     labels=c("[-15,-10]", "[-10,-5]", "[-5,0]", "[0,5]", "[5,10]", "[10,15]"))
average_stack$Change_ssp585 <- (average_stack$ssp585_mean - average_stack$Richness)
average_stack$Change_f_ssp585 <- cut(average_stack$Change_ssp585, 
                                     breaks=c(-15, -10, -5, 0, 5, 10, 15),
                                     labels=c("[-15,-10]", "[-10,-5]", "[-5,0]", "[0,5]", "[5,10]", "[10,15]"))

save(average_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_future_richness_change_", Taxon_name, ".RData"))


## Agreement: Estimate number of scenarios that follow same trends ####
# (cf. Delgado-Baquerizo et al. 2020, Soil borne pathogens, Fig. 4b)
average_stack$ssp126_gain <- 0; average_stack[average_stack$Change_ssp126>0 & !is.na(average_stack$Change_ssp126),]$ssp126_gain <- 1
average_stack$ssp370_gain <- 0; average_stack[average_stack$Change_ssp370>0 & !is.na(average_stack$Change_ssp370),]$ssp370_gain <- 1
average_stack$ssp585_gain <- 0; average_stack[average_stack$Change_ssp585>0 & !is.na(average_stack$Change_ssp585),]$ssp585_gain <- 1

average_stack$ssp126_loss <- 0; try(average_stack[average_stack$Change_ssp126<0 & !is.na(average_stack$Change_ssp126),]$ssp126_loss <- -1)
average_stack$ssp370_loss <- 0; average_stack[average_stack$Change_ssp370<0 & !is.na(average_stack$Change_ssp370),]$ssp370_loss <- -1
average_stack$ssp585_loss <- 0; average_stack[average_stack$Change_ssp585<0 & !is.na(average_stack$Change_ssp585),]$ssp585_loss <- -1

average_stack$ssp126_unchanged <- 0; try(average_stack[average_stack$Change_ssp126==0 & !is.na(average_stack$Change_ssp126),]$ssp126_unchanged <- 1)
average_stack$ssp370_unchanged <- 0; average_stack[average_stack$Change_ssp370==0 & !is.na(average_stack$Change_ssp370),]$ssp370_unchanged <- 1
average_stack$ssp585_unchanged <- 0; average_stack[average_stack$Change_ssp585==0 & !is.na(average_stack$Change_ssp585),]$ssp585_unchanged <- 1

average_stack$Gain <- rowSums(average_stack[,c("ssp126_gain", "ssp370_gain", "ssp585_gain")])
average_stack$Loss <- rowSums(average_stack[,c("ssp126_loss", "ssp370_loss", "ssp585_loss")])
average_stack$Unchanged <- rowSums(average_stack[,c("ssp126_unchanged", "ssp370_unchanged", "ssp585_unchanged")])

average_stack %>% filter(Loss==-3) %>% count()
average_stack %>% filter(Gain==3) %>% count()
average_stack %>% filter(Unchanged==3) %>% count()

average_stack$No_change <- "mixed"
average_stack[average_stack$Gain==1 & average_stack$Unchanged==2,]$No_change <- "1"
average_stack[average_stack$Loss==-1 & average_stack$Unchanged==2,]$No_change <- "-1"
average_stack[average_stack$Gain==2 & average_stack$Unchanged==1,]$No_change <- "2"
average_stack[average_stack$Loss==-2 & average_stack$Unchanged==1,]$No_change <- "-2"
average_stack[average_stack$Gain==3,]$No_change <- "3"
average_stack[average_stack$Loss==-3,]$No_change <- "-3"
average_stack[average_stack$Unchanged==3,]$No_change <- "no changes"

average_stack %>% filter(No_change==-3) %>% count()
average_stack %>% filter(No_change==3) %>% count()
average_stack %>% filter(No_change==0) %>% count()

average_stack$No_change <- factor(average_stack$No_change, levels = c("3", "2", "1", "no changes", "mixed", "-1", "-2", "-3"))

save(average_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_future_richness_agreement_", Taxon_name, ".RData"))
