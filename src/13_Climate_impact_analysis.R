#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#        Climate impact analysis            #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(raster)

library(biomod2) # also to create pseudo-absences

library(ggplot2) # for plotting the curves
library(ggpubr)

library(parallel)
library(doParallel)

# plotting
library(gridExtra)

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


#- - - - - - - - - - - - - - - - - - - - -
## Climate impact analysis ####
#- - - - - - - - - - - - - - - - - - - - -
world.inp <- map_data("world")

# load current environmental variables (for projections)
load(paste0(here::here(),"/results/EnvPredictor_5km_df_clipped.RData")) #Env_clip_df

for(subclim in c("T", "P", "TP")){
  
  average_stack <- Env_clip_df %>% dplyr::select(x, y)
  
  for(no_future in scenarioNames){
    
    load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", "2041-2070_", no_future, "_", subclim, ".RData")) #species_stack
    species_stack[,as.character(no_future)] <- species_stack$Richness
    species_stack <- species_stack[,c("x","y",no_future)]
    
    # add layer to stack
    average_stack <- full_join(average_stack, species_stack)
    
  }
  average_stack$Mean <- rowMeans(average_stack %>% dplyr::select(-x, -y), na.rm=T)
  
  save(average_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_", subclim, "_future_richness_", Taxon_name, ".RData"))
  
}

# merge into one dataframe
full_stack <- Env_clip_df %>% dplyr::select(x, y)

for(subclim in c("T", "P", "TP")){
  
  load(file=paste0(here::here(), "/results/_Maps/SDM_stack_", subclim, "_future_richness_", Taxon_name, ".RData")) #average_stack
  
  average_stack$climate <- subclim
  
  full_stack <- full_join(full_stack, average_stack)
  
}

head(full_stack)

full_stack <- full_stack %>% full_join(Env_clip_df, by=c("x", "y"))

# ANOVA 
lm1 <- lm(data=full_stack, Mean~climate)
anova(lm1)

library(emmeans)
em1 <- emmeans::emmeans(lm1, "climate", data=full_stack)
pairs(em1, adjust="tukey")

#png(paste0(here::here(), "/figures/Emmeans_lm1_climate_", Taxon_name, ".png"))
plot(em1, comparison=T)
#dev.off()

ggplot(data=full_stack, aes(x=climate, y=Mean))+
  geom_violin(width=1.4, alpha=0.7)+
  # geom_boxplot(width=0.1, color="black", fill="white", alpha=1)+
  stat_summary(fun = "mean",geom = "point",color = "black", size=3.5, show.legend = FALSE)+
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                position =  position_dodge(width = 0.9),
                width=0.1) +#geom_jitter(alpha=0.6, width=0.2)+
  theme_bw()

lm2 <- lm(data=full_stack, Mean~MAT+Dist_Coast+MAP_Seas+CEC+Elev+P+Pop_Dens+Agriculture+pH+Clay.Silt+climate)
summary(lm2)
anova(lm2)

em2 <- emmeans::emmeans(lm2, "climate", data=full_stack)
pairs(em2, adjust="tukey")

#png(paste0(here::here(), "/figures/Emmeans_lm2_climate_", Taxon_name, ".png"))
plot(em2, comparison=T)
#dev.off()

lm_varImp <- data.frame("F_value"=anova(lm2)[,"F value"])
lm_varImp$Predictor <- rownames(anova(lm2))
lm_varImp <- lm_varImp %>% filter(Predictor != "Residuals")
lm_varImp$F_abs <- abs(lm_varImp$F_value)
lm_varImp$Direction <- factor(sign(lm_varImp$F_value), 1:(-1), c("positive", "neutral", "negative"))

# load predictor table to get classification of variables
pred_tab <- readr::read_csv(file=paste0(here::here(), "/data_environment/METADATA_Predictors.csv"))

# transform to long format and add variable categories
lm_varImp <- lm_varImp%>%
  left_join(pred_tab %>% dplyr::select(Predictor, Category), by="Predictor")

# add category for clay.silt
lm_varImp[lm_varImp$Predictor=="Clay.Silt","Category"] <- "Soil"

plotTopVI <- lm_varImp %>% dplyr::select(F_abs, Predictor, Category, Direction) %>% arrange(desc(F_abs)) %>%
  ggplot(aes(x=reorder(Predictor, F_abs), y=F_abs, fill=Category)) + 
  geom_segment(aes(x=reorder(Predictor, F_abs), xend=reorder(Predictor, F_abs), y=0, yend=F_abs, lty=Direction), color="black") +
  geom_point(aes(color=Category), size=4, alpha=1) +
  coord_flip() +
  xlab("Predictors")+ylab("F value (SR)")+
  theme_bw()+theme(aspect.ratio=1/1)
plotTopVI

png(paste0(here::here(), "/figures/VariableImportance_biomod_top10_lm_climate_", Taxon_name, ".png")); plotTopVI; dev.off()

emmip(lm2, ~ climate, CIs = TRUE)

# save model summary
sink(paste0(here::here(), "/results/Summary_lm2_Crassiclitellata_climate.txt"))
print(summary(lm1))
print(anova(lm1))
print(em1)
print("###################################################")
print(summary(lm2))
print(anova(lm2))
print(em2)
sink()

# ## Paired T test
# 
# a <- ggplot(data=full_stack %>%
#               pivot_wider(id_cols=c(x,y), names_from = climate, values_from = Mean), aes(x=TP, y=P))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1, color="red")+
#   theme_bw()
# 
# b <- ggplot(data=full_stack %>%
#               pivot_wider(id_cols=c(x,y), names_from = climate, values_from = Mean), aes(x=TP, y=T))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1, color="red")+
#   theme_bw()
# 
# require(grid)
# pdf(file=paste0(here::here(), "/figures/SpeciesRichness_cert0.1_", "2041-2070_TP_xy_", Taxon_name, ".pdf"),width=10, height=10)
# gridExtra::grid.arrange(a,b)
# dev.off()
# 
# # plot future mean distribution
# png(file=paste0(here::here(), "/figures/SpeciesRichness_cert0.1_", "2041-2070_TP_", Taxon_name, ".png"),width=1000, height=1000)
# print(ggplot()+
#         geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
#         xlim(-10, 30) +
#         ylim(35, 70) +
#         
#         geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0)), 
#                   aes(x=x, y=y, fill=FutureRichness))+
#         ggtitle(paste0("Future species richness (number of species)"))+
#         scale_fill_viridis_c()+
#         geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(Richness==0)), aes(x=x, y=y), fill="grey60")+
#         theme_bw()+
#         theme(axis.title = element_blank(), legend.title = element_blank(),
#               legend.position = c(0.1,0.4)))
# dev.off()

