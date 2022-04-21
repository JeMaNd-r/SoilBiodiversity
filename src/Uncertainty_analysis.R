#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Uncertainty analysis              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load species_stack
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) 

# calculate sd of predictions
species_stack$SD <- apply(species_stack %>% dplyr::select(-x, -y, -Richness),1, sd, na.rm = TRUE)

head(species_stack)

# view uncertainty in map 
#png(file=paste0(here::here(), "/figures/Uncertainty_", Taxon_name, ".png"), width=1000, height=1000)
ggplot(data=species_stack, aes(x=x, y=y, fill=SD))+
  geom_tile()+
  ggtitle("SD between model predictions")+
  scale_fill_viridis_c()+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()
