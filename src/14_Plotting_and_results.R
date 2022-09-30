#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Plotting and other results           #
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
covarsNames <- c("MAT", "Dist_Coast", "MAP_Seas", "CEC", "Elev",
                 "P", "Pop_Dens", "Agriculture", "pH", "Clay.Silt")

# define future scenarios
scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
                               "mri-esm2-0", "ukesm1-0-ll"), "_",
                             rep(c("ssp126", "ssp370", "ssp585"),5)))


#- - - - - - - - - - - - - - - - - - - - -
## Variable importance ####
#- - - - - - - - - - - - - - - - - - - - -

var_imp <- read.csv(file=paste0(here::here(), "/results/Variable_importance_biomod_", Taxon_name, ".csv"))
var_imp

# load predictor table to get classification of variables
pred_tab <- readr::read_csv(file=paste0(here::here(), "/doc/Env_Predictors_table.csv"))

# transform to long format and add variable categories
var_imp <- var_imp %>%
  left_join(pred_tab %>% dplyr::select(Predictor, Category, Long_name) %>%
              mutate("Long_name"=str_replace_all(Long_name, "_", " ")), by="Predictor")

# add category for clay.silt
var_imp[var_imp$Predictor=="Clay.Silt","Category"] <- "Soil"
var_imp[var_imp$Predictor=="Clay.Silt","Long_name"] <- "Clay and silt content"


# plot VIF
plotVarImp <- ggplot(data=var_imp, aes(x=biomod, y=reorder(Long_name, biomod), fill=Category))+
  geom_boxplot(cex=0.2, outlier.size=1.5, show.legend = F)+
  guides(fill=(guide_legend(override.aes=list(alpha = 1))))+
  geom_point(alpha = 0, shape=21, show.legend = T)+
  geom_jitter(height=0.2, alpha=0.3, show.legend = F)+
  ylab("")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 25), axis.title = element_blank(), axis.text.x = element_text(size=15),
        legend.position=c(0.75, 0.15), legend.text = element_text(size=20), legend.title = element_blank(),
        panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank())
plotVarImp

png(paste0(here::here(), "/figures/VariableImportance_biomod_", Taxon_name, ".png"), height=600, width=700); plotVarImp; dev.off()

# plot barplot with top 10
plotTopVI <- var_imp %>% dplyr::select(biomod, Predictor, Category) %>%
  group_by(Predictor, Category) %>% summarize_all(mean, na.rm=T) %>% arrange(desc(biomod)) %>%
  ggplot(aes(x=reorder(Predictor, biomod), y=biomod, fill=Category)) + 
  geom_segment(aes(x=reorder(Predictor, biomod), xend=reorder(Predictor, biomod), y=0, yend=biomod), color="black") +
  geom_point(aes(color=Category), size=4, alpha=1) +
  coord_flip() +
  xlab("Predictors")+ylab("Mean variable importance")+
  theme_bw()+theme(aspect.ratio=1/1)
plotTopVI

png(paste0(here::here(), "/figures/VariableImportance_biomod_top10_", Taxon_name, ".png")); plotTopVI; dev.off()

# mean varImp
var_imp %>% group_by(Predictor) %>% dplyr::select(-Species, -Category) %>% summarize_all(mean)

# plot varImp of each species
var_imp$Predictor <- factor(var_imp$Predictor, levels=c("MAP_Seas", "MAT",
                                                        "Dist_Coast", "Elev",
                                                        "Agriculture", "Pop_Dens",
                                                        "CEC", "Clay.Silt", "P", "pH"))
plotAllVI <- ggplot(var_imp, aes(fill=Predictor, alpha=Predictor, y=biomod, x=reorder(Species, biomod))) + 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  xlab("Species")+
  scale_y_continuous(expand = c(0, 0))+
  scale_alpha_manual(values=c("MAP_Seas"=0.75, "MAT"=0.5, "Dist_Coast"=0.75, "Elev"=0.5,
                              "Agriculture"=0.75, "Pop_Dens"=0.5, "CEC"=0.75,"Clay.Silt"=0.55, "P"=0.35, "pH"=0.15))+
  scale_fill_manual(values=c("MAP_Seas"="#F8766D", "MAT"="#F8766D", "Dist_Coast"="#00BFC4", "Elev"="#00BFC4",
                             "Agriculture"="#7CAE00", "Pop_Dens"="#7CAE00", "CEC"="#C77CFF","Clay.Silt"="#C77CFF", "P"="#C77CFF", "pH"="#C77CFF"))+
  theme_bw()+
  theme(legend.position = "bottom")

png(paste0(here::here(), "/figures/VariableImportance_biomod_species_", Taxon_name, ".png"), height=800, width=600); plotAllVI; dev.off()


#- - - - - - - - - - - - - - - - - - - - -
## Calculate variable importance for richness (lm) ####
#- - - - - - - - - - - - - - - - - - - - -

data_stack <- average_stack %>% full_join(Env_norm_df)

lm1 <- lm(data=data_stack, Richness~MAT+Dist_Coast+MAP_Seas+CEC+Elev+P+Pop_Dens+Agriculture+pH+Clay.Silt)
summary(lm1)

lm_varImp <- data.frame("t_value"=summary(lm1)[["coefficients"]][,"t value"])
lm_varImp$Predictor <- rownames(lm_varImp)
lm_varImp <- lm_varImp %>% filter(Predictor != "(Intercept)")
lm_varImp$t_abs <- abs(lm_varImp$t_value)
lm_varImp$Direction <- factor(sign(lm_varImp$t_value), 1:(-1), c("positive", "neutral", "negative"))

# transform to long format and add variable categories
lm_varImp <- lm_varImp%>%
  left_join(pred_tab %>% dplyr::select(Predictor, Category), by="Predictor")

# add category for clay.silt
lm_varImp[lm_varImp$Predictor=="Clay.Silt","Category"] <- "Soil"

plotTopVI <- lm_varImp %>% dplyr::select(t_abs, Predictor, Category, Direction) %>% arrange(desc(t_abs)) %>%
  ggplot(aes(x=reorder(Predictor, t_abs), y=t_abs, fill=Category)) + 
  geom_segment(aes(x=reorder(Predictor, t_abs), xend=reorder(Predictor, t_abs), y=0, yend=t_abs, lty=Direction), color="black") +
  geom_point(aes(color=Category), size=4, alpha=1) +
  coord_flip() +
  xlab("Predictors")+ylab("Variable importance (SR)")+
  theme_bw()+theme(aspect.ratio=1/1)
plotTopVI

png(paste0(here::here(), "/figures/VariableImportance_biomod_top10_lm_", Taxon_name, ".png")); plotTopVI; dev.off()

# save model summary
sink(paste0(here::here(), "/results/Summary_lm1_Crassiclitellata_varImp.txt"))
print(summary(lm1))
sink()



#- - - - - - - - - - - - - - - - - - - - - -
## Load UNCERTAINTY (relevant for all plots) ####
#- - - - - - - - - - - - - - - - - - - - - -

# load uncertainty extent for all maps
load(file=paste0(here::here(), "/results/_Maps/SDM_Uncertainty_extent_", Taxon_name, ".RData")) #extent_df


#- - - - - - - - - - - - - - - - - - - - - -
## Uncertainty ####
#- - - - - - - - - - - - - - - - - - - - - -

# view uncertainty in map 
world.inp <- map_data("world")

png(file=paste0(here::here(), "/figures/Uncertainty_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(-10, 30) +
  ylim(35, 70) +
  
  geom_tile(data=uncertain_df %>% filter(Mean!=0), aes(x=x, y=y, fill=Mean))+
  ggtitle("Coefficient of variation averaged across SDMs")+
  scale_fill_viridis_c(option="E")+
  theme_bw()+  
  
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position ="bottom",legend.direction = "horizontal",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()


temp_thresh <- 0.10
png(file=paste0(here::here(), "/figures/Uncertainty_", temp_thresh, "_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(-10, 30) +
  ylim(35, 70) +
  
  geom_tile(data=uncertain_df %>% filter(Mean<temp_thresh), aes(x=x, y=y, fill=Mean))+
  ggtitle("Coefficient of variation averaged across SDMs")+
  scale_fill_viridis_c(option="E")+
  geom_tile(data=uncertain_df %>% filter(Mean>=temp_thresh), aes(x=x, y=y), fill="linen")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4),
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()


#- - - - - - - - - - - - - - - - - - - - - -
## Map species uncertainty maps ####

plots <- lapply(3:(ncol(uncertain_df)-2), function(s) {try({
  print(s-2)
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-10, 30) +
    ylim(35, 70) +
    
    geom_tile(data=uncertain_df[!is.na(uncertain_df[,s]) & uncertain_df[,s]>0,], 
              aes(x=x, y=y, fill=uncertain_df[!is.na(uncertain_df[,s]),s]))+
    ggtitle(colnames(uncertain_df)[s])+
    scale_fill_viridis_c(option="E")+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.4),
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
})
})


require(gridExtra)
#pdf(file=paste0(here::here(), "/figures/Uncertainty_allSpecies_", Taxon_name, ".pdf"))
png(file=paste0(here::here(), "/figures/Uncertainty_allSpecies_", Taxon_name, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots)
dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Species richness (current) ####
#- - - - - - - - - - - - - - - - - - - - - -

load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack

# Calculate area with 19 species
species_stack %>% filter(Richness==19 & !is.na(Richness)) %>% count()

# extract most prominent species
View(as.data.frame(colSums(species_stack, na.rm=T)) %>% arrange(colSums(species_stack, na.rm=T)))

# species richness
world.inp <- map_data("world")

png(file=paste0(here::here(), "/figures/SpeciesRichness_cert0.1_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(-10, 30) +
  ylim(35, 70) +
  
  geom_tile(data=extent_df %>% inner_join(species_stack %>% filter(Richness>0), by=c("x","y")), 
            aes(x=x, y=y, fill=Richness))+
  ggtitle("Species richness (number of species)")+
  scale_fill_viridis_c()+
  geom_tile(data=extent_df %>% inner_join(species_stack %>% filter(Richness==0), by=c("x","y")), aes(x=x, y=y), fill="grey60")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(), legend.text = element_text(size=10),
        legend.position = c(0.1,0.9), legend.direction = "horizontal",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

while (!is.null(dev.list()))  dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Species distributions (current) ####
#- - - - - - - - - - - - - - - - - - - - - -

# map binary species distributions
plots <- lapply(3:(ncol(species_stack)-1), function(s) {try({
  print(s-2)
  temp_data <- extent_df %>% inner_join(species_stack[!is.na(species_stack[,s]),])
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-10, 30) +
    ylim(35, 70) +
    
    geom_tile(data=temp_data, 
              aes(x=x, y=y, fill=as.factor(temp_data[,s])))+
    ggtitle(colnames(species_stack)[s])+
    scale_fill_manual(values=c("1"="#440154","0"="grey60","NA"="lightgrey"))+
    theme_bw()+
    guides(fill = guide_legend(# title.hjust = 1, # adjust title if needed
      label.position = "bottom",
      label.hjust = 0.5))+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.9), legend.direction = "horizontal",
          legend.text = element_text(size=20),
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
})
})

require(gridExtra)
#pdf(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_", Taxon_name, ".pdf"))
png(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_cert0.1_", Taxon_name, ".png"),width=3000, height=3300)
do.call(grid.arrange, plots)
dev.off()

while (!is.null(dev.list()))  dev.off()


#- - - - - - - - - - - - - - - - - - - - - -
## Species richness and species distributions (future) ####
#- - - - - - - - - - - - - - - - - - - - - -

world.inp <- map_data("world")

for(no_future in scenarioNames){
  
  # only plot subclim scenario TP
  subclim <- "TP"
  
  load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_2041-2070_", no_future, "_", subclim, ".RData")) #species_stack
  
  # species richness
  png(file=paste0(here::here(), "/figures/SpeciesRichness_", "2041-2070_", no_future, "_", subclim, "_", Taxon_name, ".png"),width=1000, height=1000)
  print(ggplot()+
          geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
          xlim(-10, 30) +
          ylim(35, 70) +
          
          geom_tile(data=species_stack %>% filter(Richness>0), aes(x=x, y=y, fill=Richness))+
          ggtitle(paste0("Species richness (number of species) ",  no_future, "_", subclim))+
          scale_fill_viridis_c()+
          geom_tile(data=species_stack %>% filter(Richness==0), aes(x=x, y=y), fill="grey60")+ 
          theme_bw()+
          theme(axis.title = element_blank(), legend.title = element_blank(),
                legend.position = c(0.1,0.4)))
  dev.off()
  
  #while (!is.null(dev.list()))  dev.off()
  
  
  # map binary species distributions
  plots <- lapply(3:(ncol(species_stack)-1), function(s) {try({
    print(s-2)
    ggplot()+
      geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
      xlim(-10, 30) +
      ylim(35, 70) +
      
      geom_tile(data=species_stack[!is.na(species_stack[,s]),], 
                aes(x=x, y=y, fill=as.factor(species_stack[!is.na(species_stack[,s]),s])))+
      ggtitle(colnames(species_stack)[s])+
      scale_fill_manual(values=c("1"="#440154","0"="grey","NA"="lightgrey"))+
      theme_bw()+
      theme(axis.title = element_blank(), legend.title = element_blank(),
            legend.position = c(0.1,0.4))
  })
  })
  
  require(gridExtra)
  #pdf(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_2041-2070_", no_future, "_", subclim, "_", Taxon_name, ".pdf"))
  png(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_2041-2070_", no_future, "_", subclim, "_", Taxon_name, ".png"),width=3000, height=3000)
  do.call(grid.arrange, plots)
  dev.off()
  
  #while (!is.null(dev.list()))  dev.off()
  
}

#- - - - - - - - - - - - - - - - - - - - - -
## Range decline ####
#- - - - - - - - - - - - - - - - - - - - - -

range_sum <- read.csv(file=paste0(here::here(), "/results/Range_shift_", Taxon_name, ".csv"))

# plot species range decline
png(file=paste0(here::here(), "/figures/Range_shift_", "2041-2070_", Taxon_name, ".png"), height=600, width=800)
ggplot(range_sum %>% mutate("SpeciesID"=substr(range_sum$SpeciesID, 1, 10))) +
  geom_bar(aes(x=reorder(SpeciesID, area_km2), y=area_km2_change_p*1000,  fill="[-0.5, 0.75]",), 
           stat = "identity" ,alpha=0.4, col="grey60")+
  #geom_text(aes(x=reorder(SpeciesID, area_km2), y=-290, label=reorder(SpeciesID, area_km2)), size=7)+
  
  geom_segment( aes(x=reorder(SpeciesID, area_km2), xend=reorder(SpeciesID, area_km2), y=area_km2/1000, yend=area_km2_mean/1000), color="grey") +
  geom_point( aes(x=reorder(SpeciesID, area_km2), y=area_km2_mean/1000, color="area_km2_mean"), size=7)+
  geom_point( aes(x=reorder(SpeciesID, area_km2), y=area_km2/1000, color="area_km2"),  size=7)+
  
  scale_color_manual(values = c("black", "grey60"),
                     guide  = guide_legend(), 
                     name   = "Range size in 1,000 km²",
                     labels = c("Current", "Future (mean)")) +
  scale_fill_manual(values = c("grey90"), 
                    labels = c("[-0.5, 0.75]"),
                    name = "Proportional change in range size")+
  scale_y_continuous(
    # Features of the first axis
    name = "",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./1000)) + 
  coord_flip()+
  theme_bw() +
  theme(legend.position = "top", #axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size=20), legend.title = element_text(size=25),
        axis.text = element_text(size=15),
        legend.direction = "vertical", legend.box = "horizontal",
        panel.grid.major.y = element_blank(),  panel.grid.minor.y = element_blank()) +
  xlab("")
dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
# plot both current and future range per species in one plot
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_future_species_", Taxon_name, ".RData")) #future_stack
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack

world.inp <- map_data("world")

plots <- lapply(unique(speciesNames[speciesNames$NumCells_2km>=100, "SpeciesID"]), function(s) {try({
  print(s)
  col_future <- colnames(future_stack)[stringr::str_detect(colnames(future_stack), paste0(s, ".future_mean"))]
  col_current <- colnames(species_stack)[stringr::str_detect(colnames(species_stack), paste0(s, "_current"))]
  
  temp_data <- extent_df %>% inner_join(future_stack[,c(col_future, "x", "y")] %>% 
                                          full_join(species_stack[,c(col_current, "x", "y")]))
  temp_data[,paste(s, "_change")] <- temp_data[,col_future] - temp_data[,col_current]
  temp_data[,paste(s, "_change_f")] <- as.factor(round(temp_data[,paste(s, "_change")]))
  
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-23, 40) +
    ylim(31, 75) +
    annotate(geom="text", x=-10, y=72, label=s, color="black", size=15)+
    
    geom_tile(data=temp_data[!is.na(temp_data[,paste(s, "_change")]),], 
              aes(x=x, y=y, fill=temp_data[!is.na(temp_data[,paste(s, "_change")]),paste(s, "_change")]))+
    
    #ggtitle(s)+
    # scale_fill_manual(name="Change", breaks=c("-1", "0", "1"), values=c("brown2", "#440154", "gold2"))+
    scale_fill_gradient2(low="tan1", high="deepskyblue2", mid="#440154")+
    
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.4), legend.text = element_text(size = 15))
})})
require(gridExtra)
png(file=paste0(here::here(), "/figures/DistributionMap_2041-2070_change_", Taxon_name, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots)
dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Max and min future distribution ####
#- - - - - - - - - - - - - - - - - - - - - -
world.inp <- map_data("world")

plots <- lapply(unique(speciesNames[speciesNames$NumCells_2km>=100, "SpeciesID"]), function(s) {try({
  print(s)
  col_min <- colnames(future_stack)[stringr::str_detect(colnames(future_stack), paste0(s, ".future_min"))]
  col_max <- colnames(future_stack)[stringr::str_detect(colnames(future_stack), paste0(s, ".future_max"))]
  
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-23, 60) +
    ylim(31, 75) +
    
    geom_tile(data=future_stack[!is.na(future_stack[,col_min]),], 
              aes(x=x, y=y, fill=as.factor(future_stack[!is.na(future_stack[,col_min]),col_min])))+
    
    geom_tile(data=future_stack[!is.na(future_stack[,col_max]) & future_stack[,col_max]==1,], 
              aes(x=x, y=y), fill="blue")+
    scale_fill_manual(values=c("1"="#440154","0"="grey","NA"="lightgrey"))+
    ggtitle(s)+
    
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.4))
})})
require(gridExtra)
png(file=paste0(here::here(), "/figures/DistributionMap_2041-2070_", Taxon_name, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots)
dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Future distribution ####
#- - - - - - - - - - - - - - - - - - - - - -

load(file=paste0(here::here(), "/results/_Maps/SDM_stack_future_richness_change_", Taxon_name, ".RData")) #average_stack

# plot future mean distribution
png(file=paste0(here::here(), "/figures/SpeciesRichness_cert0.1_", "2041-2070_future_", Taxon_name, ".png"),width=1000, height=1000)
print(ggplot()+
        geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
        xlim(-10, 30) +
        ylim(35, 70) +
        
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0)), 
                  aes(x=x, y=y, fill=FutureRichness))+
        ggtitle(paste0("Future species richness (number of species)"))+
        scale_fill_viridis_c()+
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(Richness==0)), aes(x=x, y=y), fill="grey60")+
        theme_bw()+
        theme(axis.title = element_blank(), legend.title = element_blank(),
              legend.position = c(0.1,0.4)))
dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
# plot change in distribution
png(file=paste0(here::here(), "/figures/SpeciesRichness_cert0.1_", "2041-2070_change_", Taxon_name, ".png"),width=1000, height=1000)
print(ggplot()+
        geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
        xlim(-10, 30) +
        ylim(35, 70) +
        
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0)), 
                  aes(x=x, y=y, fill=Change_f))+
        ggtitle(paste0("Change in species richness (number of species)"))+
        scale_fill_manual(breaks=c("[10,15]", "[5,10]", "[0,5]", "[-5,0]", "[-10,-5]", "[-15,-10]"), 
                          values=c("steelblue4", "steelblue2", "lightblue","darksalmon", "brown2", "brown4"))+
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(Change==0)), aes(x=x, y=y), fill="linen")+
        theme_bw()+
        theme(axis.title = element_blank(), legend.title = element_blank(),
              legend.position = c(0.1,0.4)))
dev.off()

# ssp126 plot change in distribution
png(file=paste0(here::here(), "/figures/SpeciesRichness_cert0.1_", "2041-2070_change_ssp126_", Taxon_name, ".png"),width=1000, height=1000)
print(ggplot()+
        geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
        xlim(-10, 30) +
        ylim(35, 70) +
        
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(Richness==0)), aes(x=x, y=y), fill="grey60")+
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0)), 
                  aes(x=x, y=y, fill=Change_f_ssp126))+
        ggtitle(paste0("Change in species richness (number of species) SSP126"))+
        scale_fill_manual(breaks=c("[10,15]", "[5,10]", "[0,5]", "[-5,0]", "[-10,-5]", "[-15,-10]"), 
                          values=c("steelblue4", "steelblue2", "lightblue","darksalmon", "brown2", "brown4"))+
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(Change_ssp126==0)), aes(x=x, y=y), fill="linen")+
        theme_bw()+
        guides(fill = guide_legend(label.position = "left", label.hjust = 1))+
        theme(axis.title = element_blank(), legend.title = element_blank(),
              legend.position = "right",
              axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()))
dev.off()

# ssp370 plot change in distribution
png(file=paste0(here::here(), "/figures/SpeciesRichness_cert0.1_", "2041-2070_change_ssp370_", Taxon_name, ".png"),width=1000, height=1000)
print(ggplot()+
        geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
        xlim(-10, 30) +
        ylim(35, 70) +
        
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(Richness==0)), aes(x=x, y=y), fill="grey60")+
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0)), 
                  aes(x=x, y=y, fill=Change_f_ssp370))+
        ggtitle(paste0("Change in species richness (number of species) SSP370"))+
        scale_fill_manual(breaks=c("[10,15]", "[5,10]", "[0,5]", "[-5,0]", "[-10,-5]", "[-15,-10]"), 
                          values=c("steelblue4", "steelblue2", "lightblue","darksalmon", "brown2", "brown4"))+
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(Change_ssp370==0)), aes(x=x, y=y), fill="linen")+
        theme_bw()+
        theme(axis.title = element_blank(), legend.title = element_blank(),
              legend.position = "right",
              axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()))
dev.off()

# ssp585 plot change in distribution
png(file=paste0(here::here(), "/figures/SpeciesRichness_cert0.1_", "2041-2070_change_ssp585_", Taxon_name, ".png"),width=1000, height=1000)
print(ggplot()+
        geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
        xlim(-10, 30) +
        ylim(35, 70) +
        
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(Richness==0)), aes(x=x, y=y), fill="grey60")+
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0)), 
                  aes(x=x, y=y, fill=Change_f_ssp585))+
        ggtitle(paste0("Change in species richness (number of species) SSP585"))+
        scale_fill_manual(breaks=c("[10,15]", "[5,10]", "[0,5]", "[-5,0]", "[-10,-5]", "[-15,-10]"), 
                          values=c("steelblue4", "steelblue2", "lightblue","darksalmon", "brown2", "brown4"))+
        geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(Change_ssp585==0)), aes(x=x, y=y), fill= "linen")+
        theme_bw()+
        theme(axis.title = element_blank(), legend.title = element_blank(),
              legend.position = "right",
              axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()))
dev.off()

#- - - - - - - - - - - - - - - - - - - - -
## Agreement between scenarios ####
#- - - - - - - - - - - - - - - - - - - - -

load(file=paste0(here::here(), "/results/_Maps/SDM_stack_future_richness_agreement_", Taxon_name, ".RData")) #average_stack

# plot number scenarios that predict gain/loss/no change
png(file=paste0(here::here(), "/figures/SpeciesRichness_cert0.1_", "2041-2070_change_noScenarios_", Taxon_name, ".png"), width=1100, height=1000)
print({ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-10, 30) +
    ylim(35, 70) +
    
    #geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(Richness==0)), aes(x=x, y=y), fill="grey60")+
    geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(!is.na(Change))), 
              aes(x=x, y=y, fill=No_change))+
    ggtitle(paste0("Agreement between SSP scenarios"))+
    scale_fill_manual(values=c("steelblue4", "steelblue2", "lightblue", "linen", "sandybrown", "darksalmon", "brown2", "brown4"))+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = "right",
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())})
dev.off()

#- - - - - - - - - - - - - - - - - - - - -
## Which species are lost? ####
#- - - - - - - - - - - - - - - - - - - - -

temp_coords <- average_stack %>% filter(Change<0) %>% dplyr::select(x,y)
temp_coords <- future_stack %>% inner_join(temp_coords)

head(temp_coords)

# add current distributions
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
temp_coords <- temp_coords %>% left_join(species_stack)

head(temp_coords)

temp_species <- unique(substr(colnames(temp_coords %>% dplyr::select(-x, -y, -colnames(temp_coords)[str_detect(colnames(temp_coords),"Richness")])), 1, 10))
for(spID in temp_species){
  temp_coords[,paste0(temp_species, "_change")] <- temp_coords[,paste0(temp_species, ".future_mean")] - temp_coords[,paste0(temp_species, "_current")] 
}

head(temp_coords)
summary(temp_coords$Allol_chlo_change)

#colSums(temp_coords[, colnames(temp_coords)[str_detect(colnames(temp_coords), "_change")]], na.rm=T) %>% as.data.frame()
temp_coords %>% count()
colSums(temp_coords[, colnames(temp_coords)[str_detect(colnames(temp_coords), "_change")]] < 0 , na.rm=T)
temp_coords %>% filter(Allol_chlo_change==1) %>% count()
temp_coords %>% filter(Allol_chlo_change==-1) %>% count()


#- - - - - - - - - - - - - - - - - - - - -
## Model performance ####
#- - - - - - - - - - - - - - - - - - - - -

data_eval <- read.csv(paste0(here::here(), "/results/Model_evaluation_", Taxon_name, ".csv"))

data_eval %>% dplyr::select(-SpeciesID) %>% summarize_all(mean)
data_eval %>% dplyr::select(-SpeciesID) %>% summarize_all(sd)
