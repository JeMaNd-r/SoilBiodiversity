#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Calculate species richness           #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# change temporary directory for files
raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - - -
## Create empty data frame for all species ####
# load environmental space as data frame
load(paste0(here::here(),"/results/EnvPredictor_2km_df_normalized.RData")) #Env_norm_df

# create empty data frame
species_stack <- Env_norm_df %>% dplyr::select(x, y)

# load model evaluation
mod_eval_all <- read.csv(file=paste0(here::here(), "/results/ModelEvaluation_", Taxon_name, ".csv"))

# for loop through all species
for(spID in unique(speciesNames$SpeciesID)){ try({

  # read model performance evaluation table (for threshold of best model)
  mod_eval <- mod_eval_all[mod_eval_all$species==spID,]
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Load probability maps ####
  load(file=paste0(here::here(), "/results/_Maps/SDM_bestPrediction_", Taxon_name, "_", spID, ".RData")) #best_pred
  
  # extract name of best model
  temp_model <- mod_eval[mod_eval$best==1 & !is.na(mod_eval$best), "model"]
  
  print(paste0(spID, " successfully loaded. Best model is ", temp_model))
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Transform to binary maps ####
  
  # extract threshold to define presence/absence
  temp_thresh <- mod_eval[mod_eval$model==temp_model,"thres.maxTSS"]
  if(is.na(temp_thresh)) temp_tresh <- 0.9
  
  # change to binary
  best_pred[best_pred$layer>=temp_thresh & !is.na(best_pred$layer), "layer"] <- 1
  best_pred[best_pred$layer<temp_thresh & !is.na(best_pred$layer), "layer"] <- 0
  
  best_pred[,paste0(spID,"_", temp_model)] <- best_pred$layer
  best_pred <- best_pred[,c("x","y",paste0(spID,"_", temp_model))]
  
  # save binary
  save(best_pred, file=paste0(here::here(), "/results/_Maps/SDM_bestPrediction_binary_", Taxon_name, "_", spID, ".RData"))
  
  print(paste0("Saved binary prediction of ", spID))
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Stack species binary maps ####
  
  # add species dataframe to stacked dataframe
  species_stack <- species_stack %>% full_join(best_pred, by=c("x","y"))
  
  print(paste0("Added binary prediction of ", spID, " to the species stack"))
  
  rm(temp_thresh, best_pred, temp_model)
}, silent=T)}  

head(species_stack)


#- - - - - - - - - - - - - - - - - - - - - -
## Calculate richness ####
species_stack$Richness <- rowSums(species_stack %>% dplyr::select(-x, -y), na.rm=T)

#- - - - - - - - - - - - - - - - - - - - - -
## Save species stack ####
save(species_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData"))


#- - - - - - - - - - - - - - - - - - - - - -
## View individual binary maps and species stack ####
# species richness
png(file=paste0(here::here(), "/figures/SpeciesRichness_", Taxon_name, ".png"),width=1000, height=1000)
ggplot(data=species_stack, aes(x=x, y=y, fill=Richness))+
  geom_tile()+
  ggtitle("Species richness (number of species)")+
  scale_fill_viridis_c()+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()

while (!is.null(dev.list()))  dev.off()


# map binary species distributions
plots <- lapply(3:(ncol(species_stack)-1), function(s) {try({
  print(s-2)
  ggplot(data=species_stack, aes(x=x, y=y, fill=as.factor(species_stack[,s])))+
    geom_tile()+
    ggtitle(colnames(species_stack)[s])+
    scale_fill_manual(values=c("1"="#fde725","0"="#440154","NA"="lightgrey"))+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.4))


require(gridExtra)
#pdf(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_", Taxon_name, ".pdf"))
png(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_", Taxon_name, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots)
dev.off()

while (!is.null(dev.list()))  dev.off()

