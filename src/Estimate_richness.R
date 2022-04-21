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

# for loop through all species
for(sp in unique(speciesNames$SpeciesID)){ try({

  # read model performance evaluation table (for threshold of best model)
  mod_eval <- read.csv(file=paste0(here::here(), "/results/ModelEvaluation_", Taxon_name, "_", sp, ".csv"))
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Load probability maps ####
  load(file=paste0(here::here(), "/results/_Maps/SDM_bestPrediction_", Taxon_name, "_", sp, ".RData")) #best_pred
  
  # extract name of best model
  temp_model <- mod_eval[mod_eval$tss == max(mod_eval$tss, na.rm=T) & !is.na(mod_eval$model), "model"]
  
  print(paste0(sp, " successfully loaded. Best model is ", temp_model))
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Transform to binary maps ####
  
  # extract threshold to define presence/absence
  temp_thresh <- mod_eval[mod_eval$model==temp_model,"thres.maxTSS"]
  
  # change to binary
  best_pred[best_pred$layer>=temp_thresh & !is.na(best_pred$layer), "layer"] <- 1
  best_pred[best_pred$layer<temp_thresh & !is.na(best_pred$layer), "layer"] <- 0
  
  best_pred[,sp] <- best_pred$layer
  best_pred <- best_pred[,c("x","y",sp)]
  
  # save binary
  save(best_pred, file=paste0(here::here(), "/results/_Maps/SDM_bestPrediction_binary_", Taxon_name, "_", sp, ".RData"))
  
  print(paste0("Saved binary prediction of ", sp))
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Stack species binary maps ####
  
  # add species dataframe to stacked dataframe
  species_stack <- species_stack %>% full_join(best_pred, by=c("x","y"))
  
  print(paste0("Added binary prediction of ", sp, " to the species stack"))
  
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
#png(file=paste0(here::here(), "/figures/SpeciesRichness_", Taxon_name, ".png"),width=1000, height=1000)
ggplot(data=species_stack, aes(x=x, y=y, fill=Richness))+
  geom_tile()+
  ggtitle("Species richness (number of species)")+
  scale_fill_viridis_c()+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()


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
})})

require(gridExtra)
#pdf(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_", Taxon_name, ".pdf"))
#png(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_", Taxon_name, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots)
dev.off()


