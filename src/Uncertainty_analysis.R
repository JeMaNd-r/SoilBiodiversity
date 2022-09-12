#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Uncertainty analysis              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load environmental space as data frame
load(paste0(here::here(),"/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df

uncertain_df <- Env_norm_df %>% dplyr::select(x, y)

for(spID in unique(speciesNames$SpeciesID)){try({
  
  print(paste0("Species: ", spID))

  # list files in species-specific BIOMOD folder
  temp_files <- list.files(paste0(here::here(), "/results/", Taxon_name, "/", stringr::str_replace(spID, "_", "."), "/proj_modeling"), full.names = TRUE)
          
  #setwd(paste0(here::here(), "/results/", Taxon_name))

  myBiomodEnProj <- get(load(temp_files[stringr::str_detect(temp_files,"ensemble.RData")]))
       
  # save predictions as raster file
  temp_prediction <- myBiomodEnProj[,1]
  temp_prediction <- as.numeric(temp_prediction)
  # add names of grid cell (only for those that have no NA in any layer)
  names(temp_prediction) <- rownames(Env_norm_df)
  temp_prediction <- as.data.frame(temp_prediction)
  temp_prediction$x <- Env_norm_df$x
  temp_prediction$y <- Env_norm_df$y
  temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
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

# view uncertainty in map 
world.inp <- map_data("world")

png(file=paste0(here::here(), "/figures/Uncertainty_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-23, 40) +
    ylim(31, 75) +

  geom_tile(data=uncertain_df %>% filter(Mean!=0), aes(x=x, y=y, fill=Mean))+
  ggtitle("Coefficient of variation averaged across SDMs")+
  scale_fill_viridis_c(option="E")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()

# extract area with uncertainty lower than threshold
summary(uncertain_df$Mean)

extent_df <- uncertain_df %>% filter(Mean<0.1 & !is.na(Mean)) %>% dplyr::select(x,y)
save(extent_df, file=paste0(here::here(), "/results/_Maps/SDM_Uncertainty_extent_", Taxon_name, ".RData"))

temp_thresh <- 0.10
png(file=paste0(here::here(), "/figures/Uncertainty_", temp_thresh, "_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(-23, 40) +
  ylim(31, 75) +
  
  geom_tile(data=uncertain_df %>% filter(Mean<temp_thresh), aes(x=x, y=y, fill=Mean))+
  ggtitle("Coefficient of variation averaged across SDMs")+
  scale_fill_viridis_c(option="E")+
  geom_tile(data=uncertain_df %>% filter(Mean>=temp_thresh), aes(x=x, y=y), fill="linen")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()


#- - - - - - - - - - - - - - - - - - - - - -
## Map species uncertainty maps ####

plots <- lapply(3:(ncol(uncertain_df)-2), function(s) {try({
  print(s-2)
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-23, 60) +
    ylim(31, 75) +
    
    geom_tile(data=uncertain_df[!is.na(uncertain_df[,s]) & uncertain_df[,s]>0,], 
              aes(x=x, y=y, fill=uncertain_df[!is.na(uncertain_df[,s]),s]))+
    ggtitle(colnames(uncertain_df)[s])+
    scale_fill_viridis_c(option="E")+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.4))
  })
})


require(gridExtra)
#pdf(file=paste0(here::here(), "/figures/Uncertainty_allSpecies_", Taxon_name, ".pdf"))
png(file=paste0(here::here(), "/figures/Uncertainty_allSpecies_", Taxon_name, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots)
dev.off()

