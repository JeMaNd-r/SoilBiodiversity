#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Uncertainty analysis              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load environmental space as data frame
load(paste0(here::here(),"/results/EnvPredictor_2km_df_normalized.RData")) #Env_norm_df


for(spID in unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID)){ try({
print(paste0("Species: ", spID))

# load species-specific models
load(file=paste0(here::here(), "/results/_Maps/SDM_Predictions_", Taxon_name, "_", spID, ".RData")) #model_list

# create empty data frame
uncertain_df <- Env_norm_df %>% dplyr::select(x, y)

# loop through all models
for(m in 1:length(model_list)){ try({
  # load model prediction
  temp_model <- model_list[[m]]
  
  # rename layer column
  temp_model[,names(model_list)[m]] <- temp_model$layer
  temp_model <- temp_model[,c("x","y",names(model_list)[m])]
  
  # add species dataframe to stacked dataframe
  uncertain_df <- uncertain_df %>% full_join(temp_model, by=c("x","y"))
  
})}

head(uncertain_df)

# calculate sd of predictions
uncertain_df$SD <- apply(uncertain_df %>% dplyr::select(-x, -y), 1, sd, na.rm = TRUE)

head(uncertain_df)

# save species' uncertainty map
save(uncertain_df, file=paste0(here::here(), "/results/_Maps/SDM_Uncertainty_", Taxon_name, "_", spID, ".RData"))

# view uncertainty in map 
# png(file=paste0(here::here(), "/figures/Uncertainty_", Taxon_name, ".png"), width=1000, height=1000)
# ggplot(data=uncertain_df, aes(x=x, y=y, fill=SD))+
#   geom_tile()+
#   ggtitle("SD between model predictions")+
#   scale_fill_viridis_c()+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))
# dev.off()

})}

#- - - - - - - - - - - - - - - - - - - - - -
## Stack species uncertainty maps ####

# create empty data frame
species_stack <- Env_norm_df %>% dplyr::select(x, y)

# for loop through all species
for(spID in unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID)){ try({

  load(file=paste0(here::here(), "/results/_Maps/SDM_Uncertainty_", Taxon_name, "_", spID, ".RData")) #uncertain_df
  uncertain_df[,paste0("SD.", spID)] <- uncertain_df$SD

  # add species dataframe to stacked dataframe
  species_stack <- species_stack %>% full_join(uncertain_df %>% dplyr::select(x,y,paste0("SD.", spID)), by=c("x","y"))
  
  print(paste0("Added binary prediction of ", spID, " to the species stack"))
  
  rm(uncertain_df)
}, silent=T)}

head(species_stack)

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate mean uncertainty ####
species_stack$SD_mean <- rowMeans(species_stack %>% dplyr::select(-x, -y), na.rm=T)

# save
save(species_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_uncertainty_", Taxon_name, ".RData"))

#- - - - - - - - - - - - - - - - - - - - - -
## Plotting ####
# mean uncertainty
png(file=paste0(here::here(), "/figures/Uncertainty_", Taxon_name, ".png"),width=1000, height=1000)
ggplot(data=species_stack, aes(x=x, y=y, fill=SD_mean))+
  geom_tile()+
  ggtitle("Mean standard deviation between species-spcific algorithms")+
  scale_fill_viridis_c(option="magma")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()

while (!is.null(dev.list()))  dev.off()


