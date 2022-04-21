#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Uncertainty analysis              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load species-specific models
load(file=paste0(here::here(), "/results/_Maps/SDM_Predictions_", Taxon_name, "_", spID, ".RData")) #model_list

# load environmental space as data frame
load(paste0(here::here(),"/results/EnvPredictor_2km_df_normalized.RData")) #Env_norm_df

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

# view uncertainty in map 
#png(file=paste0(here::here(), "/figures/Uncertainty_", Taxon_name, ".png"), width=1000, height=1000)
ggplot(data=uncertain_df, aes(x=x, y=y, fill=SD))+
  geom_tile()+
  ggtitle("SD between model predictions")+
  scale_fill_viridis_c()+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()


