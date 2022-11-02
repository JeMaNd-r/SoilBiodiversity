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
covarsNames <- c("MAT", "MAP_Seas", "Dist_Coast", "Agriculture", "pH", 
                 "P", "CEC", "Elev", "Clay.Silt", "Pop_Dens")

# define future scenarios
scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
                               "mri-esm2-0", "ukesm1-0-ll"), "_",
                             rep(c("ssp126", "ssp370", "ssp585"),5)))

# geographic extent of Europe
extent_Europe <- c(-23, 60, 31, 75)

# load background map
world.inp <- map_data("world")

#- - - - - - - - - - - - - - - - - - - - -
## Occurrences ####
#- - - - - - - - - - - - - - - - - - - - -

# load raw data
occ_raw <- read.csv(file=paste0(here::here(), "/intermediates/Earthworm_occurrence_GBIF-sWorm-Edapho-SoilReCon-JM.csv"))

# load summary number of records during processing
occ_process <- read.csv(file=paste0(here::here(), "/results/NoRecords_summary_Crassiclitellata.csv"))
occ_process <- occ_process %>%
  filter(!str_detect(Subset,"species")) %>%
  filter(!str_detect(Subset,"[_]")) %>%
  filter(!str_detect(Subset,"merged")) %>%
  filter(!str_detect(Subset,"total")) 

## barplot number of occurrences per datasource
pdf(paste0(here::here(), "/figures/OccurrenceRaw_perDatasource_barplot.pdf"), width=10)
ggplot(data=occ_process, 
       aes(x=reorder(Subset, NumberRecords), y=NumberRecords, fill=ProcessingStep))+
  geom_bar(stat="identity", position="dodge")+
  geom_text(label=occ_process$NumberRecords, hjust=-0.5, position=position_dodge(width=1))+
  xlab("")+ ylab("Number of occurrence records")+
  scale_y_continuous(expand=c(0,0), limits=c(0,65000))+
  coord_flip()+
  theme_bw()
dev.off()

## plot raw occurrences colored by year
pdf(paste0(here::here(), "/figures/OccurrenceRaw_perYear.pdf"), width=10)
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "white")+
  geom_point(data=occ_raw, aes(x=longitude, y=latitude, col=year),cex=0.3)+ theme_bw()+
  xlim(min(extent_Europe[1], na.rm = T), max(extent_Europe[2], na.rm = T)) +
  ylim(min(extent_Europe[3], na.rm = T), max(extent_Europe[4], na.rm = T)) +
  scale_color_steps2(breaks=c(1970, 1980, 1990, 2000, 2010), midpoint=1995, 
                     high="#10a53dFF", mid="#ffcf20FF", low="#541352FF")+
  theme(panel.background = element_rect(fill = "grey80",
                                        colour = "grey80",
                                        size = 0.5, linetype = "solid"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey80"))
dev.off()

# load cleaned occurrence records 
occ_clean <- read.csv(file=paste0(here::here(), "/results/Occurrences_", Taxon_name, ".csv"))

# load matrix containing information on number of occurrence records in grid
occ_points <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))

# calculate number of records per datasource
full_join(occ_raw %>% group_by(datasource) %>% count(name="raw"), 
          occ_clean %>% group_by(datasource) %>% count(name="clean"))

# load number removed records during cleaning
read.csv(file=paste0(here::here(), "/results/NoRecords_cleaning_", Taxon_name, ".csv"))

# count occurrences per species & data source
count_data <- occ_raw %>% 
  full_join(speciesNames[,c("Species", "SpeciesID")], by=c("species"="Species")) %>%
  group_by(datasource, SpeciesID) %>% count() %>%
  pivot_wider(names_from = datasource, values_from = n) %>%
  full_join(occ_clean %>% group_by(datasource, SpeciesID) %>% count() %>%
              pivot_wider(names_from = datasource, values_from = n), suffix = c("_raw", "_clean"), by="SpeciesID")
count_data$RawOcc <- rowSums(count_data[,2:7], na.rm=T)
count_data$CleanOcc <- rowSums(count_data[,8:13], na.rm=T)

# add number of records after
count_data <- count_data %>% 
  full_join(speciesNames[,c("SpeciesID", #"Acc_name", "Species_final", "Species", #"NumCells_1km", 
                             "NumCells_2km")], by=c("SpeciesID")) %>%
  filter(!is.na(SpeciesID))

# replace NA with 0
#count_data[is.na(count_data)] <- 0

# add info if species was analysed or not
count_data$Included <- FALSE
count_data[count_data$NumCells_2km >=10, "Included"] <- TRUE 

# sort by included or not, and have a look
count_data <- count_data %>% arrange(desc(Included), SpeciesID) %>%
  dplyr::select(SpeciesID, RawOcc, CleanOcc, NumCells_2km, everything()) %>%
  unique()
count_data

# save
write.csv(count_data, file=paste0(here::here(), "/results/NoRecords_perSpecies_full_", Taxon_name, ".csv"), row.names = F)

count_data <- read.csv(file=paste0(here::here(), "/results/NoRecords_perSpecies_full_", Taxon_name, ".csv"))

# look at speciesID with most records
count_data %>% arrange(desc(NumCells_2km), SpeciesID)
# Eisen_tetr, Aporr_cali, Aporr_rose, Lumbr_terr, Lumbr_rube...

count_data2 <- count_data %>% full_join(speciesNames %>%
                                          dplyr::select(-NumCells_2km, -Group_name), by=c("SpeciesID"))
count_data2

## plot count occurrences
ggplot(count_data2 %>% filter(!is.na(SpeciesID), Included=TRUE), aes(x=RawOcc-CleanOcc, y=SpeciesID, fill=Ecogroup))+
  geom_bar(stat = "identity")

## plot total species' occurrences
plotOccRaw <- ggplot()+ 
  geom_map(data = world.inp, map = world.inp, 
           aes(map_id = region), fill = "white")+
  xlim(min(extent_Europe[1], na.rm = T), 40) +
  ylim(min(extent_Europe[3], na.rm = T), max(extent_Europe[4], na.rm = T)) +
  geom_point(data=occ_clean, 
             aes(x=longitude, y=latitude, color=datasource), 
             cex=0.3, shape=".")+
  theme_bw()+
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill = "grey80",
                                        colour = "grey80",
                                        size = 0.5, 
                                        linetype = "solid"))+
  guides(color = guide_legend(override.aes = list(size = 3))) #makes legend icons bigger

plotOccRaw

# pdf(paste0(here::here(), "/figures/CleanOccurrences_", Taxon_name, "_perDatasource.pdf")); plotOccRaw; dev.off()

rm(plotOccRaw)

occ_clean %>% group_by(datasource) %>% count() 
# edapho 13054, gbif 54883, jean 24860, jerome 732, soilrecon 175, sworm 5028

## plot in Germany
german.inp <- map_data("world", "Germany")

# plot total species' occurrences
plotOccRawGER <- ggplot()+ #, alpha=`Number of Records`
  #geom_polygon(data=bg.map)+
  geom_map(data=world.inp, map = world.inp, aes(map_id = region), fill="grey90") +
  geom_map(data=german.inp, map = german.inp, aes(map_id = region), fill="white")+
  xlim(5, 17) +
  ylim(46,57) +
  
  geom_point(data=occ_clean, aes(x=longitude, y=latitude, color=datasource, shape="."), cex=0.4)+
  #scale_x_continuous(limits=c(5, 17))+ 
  #scale_y_continuous(limits=c(46,57))+
  theme_bw()+
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill = "grey70",
                                        colour = "grey70",
                                        size = 0.5, 
                                        linetype = "solid"))+
  guides(color = guide_legend(override.aes = list(size = 3))) #makes legend icons bigger

plotOccRawGER

# pdf(paste0(here::here(), "/figures/RawOccurrences_", Taxon_name, "_perDatasource_GER.pdf")); plotOccRawGER; dev.off()

rm(plotOccRawGER, occ_clean)

# 
# # calculate raw species richness
# #### needs to be fixed ######
# occ_rich <- occ_points %>% 
#   group_by(Latitude = round(x,0), Longitude=round(y,0)) %>%
#   summarise_at(vars(colnames(occ_points %>% dplyr::select(-x, -y))), mean, na.rm=T)
# occ_rich$Richness <- apply(occ_rich > 0, 1, sum, na.rm=T)
# 
# # plot total species' occurrences
# plotOcc <- ggplot()+ 
#   geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
#   xlim(min(extent_Europe[1], na.rm = T), max(extent_Europe[2], na.rm = T)) +
#   ylim(min(extent_Europe[3], na.rm = T), max(extent_Europe[4], na.rm = T)) +
#   
#   geom_point(data=occ_rich %>%
#                dplyr::select(c(Latitude, Longitude, Richness)),
#              aes(x=Latitude, y=Longitude, color=Richness))+ #, alpha=`Number of Species`
#   
#   scale_color_gradient2(5,    # provide any number of colors
#                         low = "black", high="orange", mid= "blue",
#                         midpoint = 10,
#                         #values = scales::rescale(c(1, 2, 3, 5, 10, 30)), 
#                         breaks = c(1, 2, 5, 10, 20, 30), limits=c(0,30))+
#   theme_bw()+
#   theme(legend.position = "bottom", legend.text = element_text(size=8), legend.key.width = unit(2, "cm"))
# 
# plotOcc

# calculate individual species' occurrences
occ_points_species <- occ_points %>% 
  pivot_longer(cols=speciesNames$SpeciesID[speciesNames$SpeciesID %in% colnames(occ_points)], 
               names_to = "SpeciesID") %>% 
  mutate("Latitude"=round(x,0), "Longitude"=round(y,0)) %>%
  group_by(Latitude, Longitude, SpeciesID) %>%
  filter(!is.na(value)) %>%
  summarize("Number of Records"= n(), .groups="keep") %>%
  filter("Number of Records" > 0) 

# only keep species that will be analyzed (i.e., present in at least 5 grid cells)
occ_points_species <- occ_points_species[occ_points_species$SpeciesID %in%
                                           speciesNames[speciesNames$NumCells_2km >=5, "SpeciesID"],]
occ_points_species

## plot some of the individual species' occurrences
plotOccSpecies <- ggplot(occ_points_species, 
                         aes(x=Latitude, y=Longitude, color=`Number of Records`, group=SpeciesID))+
  #geom_polygon(data=bg.map)
  geom_point(cex=0.015, pch=15)+
  facet_wrap(vars(SpeciesID))+
  scale_color_gradientn(    # provide any number of colors
    colors = c("black", "blue", "orange"),
    values = scales::rescale(c(1, 5, 20, 30, 50, 100, 300)), 
    breaks = c(5, 20, 50, 100, 200))+
  
  # add number of grid cells in which the species is present
  geom_text(data=occ_points_species %>% group_by(SpeciesID) %>% summarize("n"=sum(`Number of Records`)), 
            aes(x=30, y=33, label=paste0("n=", n)), color="black", 
            inherit.aes=FALSE, parse=FALSE, cex=0.7)+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks.length = unit(0, "cm"),
        legend.position = "bottom", legend.text = element_text(size=5))
plotOccSpecies

# pdf(paste0(here::here(), "/figures/GriddedOccurrences_", Taxon_name, "_perSpecies.pdf")); plotOccSpecies; dev.off()

rm(plotOccSpecies, occ_points_species)

## Barplot: occurrence per year
summary(occ_points$year)

ggplot(data=occ_points, aes(x=year))+
  geom_bar()+theme_bw()


## Calculate number of occurrences per country
library(rworldmap)
m <- rworldmap::getMap()

occ_points_sp <- occ_points
occ_points_sp <- occ_clean

coordinates(occ_points_sp) <- ~x+y
proj4string(occ_points_sp) = proj4string(m)

# extract country names
occ_country <- droplevels(over(occ_points_sp,m)$NAME)
unique(occ_country)

table(occ_country)

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

# load uncertainty
load(file=paste0(here::here(), "/results/_Maps/SDM_Uncertainty_", Taxon_name, ".RData")) #uncertain_df


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


temp_thresh <- 0.1
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
    #ggtitle(colnames(uncertain_df)[s])+
    annotate(geom="text", x=-3, y=68, label=colnames(uncertain_df)[s], color="black", size=15)+
    scale_fill_viridis_c(option="E")+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.8), legend.direction = "horizontal",
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
  geom_errorbar(aes(x=reorder(SpeciesID, desc(SpeciesID)),
                    ymin=(area_km2_change_p*1000)-(area_km2_p_sd*1000),
                    ymax=(area_km2_change_p*1000)+(area_km2_p_sd*1000)),
                position=position_dodge(width=0.9), width=0.4,alpha=0.4)+
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
  temp_data[temp_data[,col_future]==0 & temp_data[,col_current]==0,paste(s, "_change")] <- NA
  temp_data[,paste(s, "_change_f")] <- as.factor(round(temp_data[,paste(s, "_change")]))
  
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-10, 30) +
    ylim(35, 70) +
    annotate(geom="text", x=-3, y=68, label=s, color="black", size=15)+
    
    geom_tile(data=temp_data[!is.na(temp_data[,paste(s, "_change")]),], 
              aes(x=x, y=y, fill=temp_data[!is.na(temp_data[,paste(s, "_change")]),paste(s, "_change")]))+
    
    #ggtitle(s)+
    # scale_fill_manual(name="Change", breaks=c("-1", "0", "1"), values=c("brown2", "#440154", "gold2"))+
    scale_fill_gradient2(low="tan1", high="deepskyblue2", mid="#440154")+
    
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.8), legend.text = element_text(size = 15), legend.direction = "horizontal",
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
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


# where do scenarios disagree?
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & Gain==3)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & Loss==-3)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & No_change=="mixed")))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & No_change=="mixed")))*5

nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & Unchanged==1)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & ssp126_unchanged==1)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & ssp370_unchanged==1)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & ssp585_unchanged==1)))/nrow(extent_df)

nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & Change_ssp126==0)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & Change_ssp370==0)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & Change_ssp585==0)))/nrow(extent_df)

nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & ssp126_loss!=0)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & ssp370_loss!=0)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & ssp585_loss!=0)))/nrow(extent_df)

nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & ssp126_gain!=0)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & ssp370_gain!=0)))/nrow(extent_df)
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & ssp585_gain!=0)))/nrow(extent_df)

# where are all species lost?
nrow(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & Richness>0 & FutureRichness==0)))*5
summary(extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & Richness>0 & FutureRichness==0)) %>% dplyr::select(Richness))

ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(-10, 30) +
  ylim(35, 70) +
  
  geom_tile(data=extent_df %>% inner_join(average_stack %>% filter(!is.na(Change) & Richness>0 & FutureRichness==0)), 
            aes(x=x, y=y, fill=Richness))+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

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

mod_eval <- read.csv(file=paste0(here::here(), "/results/ModelEvaluation_", Taxon_name, ".csv"))

# point plot with lables, tss over roc
pdf(paste0(here::here(), "/figures/Model_performance_", Taxon_name, "_tss-roc_perSpecies.pdf"))
ggplot(mod_eval %>% filter(!is.na(species)), aes(x=tss, y=roc, color=model))+
  geom_text(label=mod_eval[!is.na(mod_eval$species),"model"], nudge_x = 0, nudge_y = 0, check_overlap = F, cex=1)+
  facet_wrap(vars(species))+
  xlim(0,1)+
  theme_bw()
dev.off()

# boxplot, tss per algorithm
pdf(paste0(here::here(), "/figures/Model_performance", Taxon_name, "_boxplot.pdf"))
ggplot(mod_eval %>% filter(!is.na(species)), aes(x=tss, y=model))+
  geom_boxplot()+
  xlim(0,1)+
  theme_bw()
dev.off()

#- - - - - - - - - - - - - - - - - - - - -
## Protection of ranges ####
#- - - - - - - - - - - - - - - - - - - - -

cover_df <- read.csv(file=paste0(here::here(), "/results/ProtectionStatus_", Taxon_name, ".csv"))
cover_df$IUCNcat <- factor(cover_df$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))
load(file=paste0(here::here(), "/results/ProtectionStatus_SR_SSPs_", Taxon_name, ".csv")) #cover_sr
cover_sr$IUCNcat <- factor(cover_sr$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))

# boxplot of percent area covered by PA overall
a <- ggplot(data=cover_sr %>%  filter(IUCNcat!="Presence" & IUCNcat!="Protected"), aes(x=reorder(IUCNcat, IUCNcat), y=current_mean, fill=IUCNcat))+
  geom_violin(width=1.4, alpha=0.7)+
  ggtitle("Current")+
  # geom_boxplot(width=0.1, color="black", fill="white", alpha=1)+
  stat_summary(fun = "mean",geom = "point",color = "black", size=3.5, show.legend = F)+
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                position =  position_dodge(width = 0.9),
                width=0.1) +
  #geom_jitter(alpha=0.6, width=0.2)+
  theme_bw()+ 
  xlab("Type of protected area (IUCN categories)")+ ylab("Number of species")+
  scale_fill_manual(values=c("olivedrab1","olivedrab3", "olivedrab4", "darkolivegreen", "goldenrod3","goldenrod1", "gold1",
                             "lemonchiffon2","gainsboro", "grey","white" ))+
  scale_y_continuous(expand=c(0,0))+
  theme(legend.position="left", axis.text.x=element_text(angle=30, hjust=1),
        legend.text = element_text(size=5), legend.title = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

# bar chart of percent area covered by PA per species
b <- ggplot(data=data_barplot, 
            aes(y=coverage, x=reorder(SpeciesID, desc(SpeciesID)), fill=reorder(IUCNcat, desc(IUCNcat))))+
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  xlab("Species")+ ylab("Proportion of range covered by protected area network")+
  coord_flip()+
  #scale_fill_viridis_d()+
  scale_fill_manual(values=rev(c("olivedrab1","olivedrab3", "olivedrab4", "darkolivegreen", "goldenrod3","goldenrod1", "gold1",
                                 "lemonchiffon2","gainsboro", "grey" )))+
  scale_y_continuous(expand=c(0,0), limits=c(0,0.65))+
  geom_vline(xintercept=20.5, lty=2)+
  theme(legend.position="none",legend.text = element_text(size=5), legend.title = element_text(size=5),
        axis.text.y = element_text(size=15),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

png(paste0(here::here(), "/figures/ProtectionStatus_current_", Taxon_name, ".png"))
grid.arrange(a, b, heights=c(1,1))
dev.off()

## Load future protection
cover_df <- read.csv(file=paste0(here::here(), "/results/ProtectionStatus_SSPs_", Taxon_name, ".csv"))
cover_df_current <- read.csv(file=paste0(here::here(), "/results/ProtectionStatus_", Taxon_name, ".csv"))

load(file=paste0(here::here(), "/results/ProtectionStatus_SR_SSPs_", Taxon_name, ".csv")) #cover_sr

cover_df$SpeciesID <- substr(cover_df$SpeciesID, 1, 10)

cover_df <- cover_df %>% full_join(cover_df_current %>% mutate("SSP"="current"))

cover_df$IUCNcat <- factor(cover_df$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))
cover_sr$IUCNcat <- factor(cover_sr$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))

# boxplot of percent area covered by PA overall
boxplot_template <- function(data, col_name){
  data$layer <- c(as.numeric(data[,paste0(col_name, "_mean")] %>% unlist()))
  ggplot(data=data %>%  filter(IUCNcat!="Presence" & IUCNcat!="Protected"), aes(x=reorder(IUCNcat, IUCNcat), y=layer, fill=IUCNcat))+
    geom_violin(width=1.4, alpha=0.7)+
    ggtitle(col_name)+
    # geom_boxplot(width=0.1, color="black", fill="white", alpha=1)+
    stat_summary(fun = "mean",geom = "point",color = "black", size=3.5, show.legend = FALSE)+
    geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                  fun.args = list(mult = 1),
                  position =  position_dodge(width = 0.9),
                  width=0.1) +#geom_jitter(alpha=0.6, width=0.2)+
    theme_bw()+ 
    xlab("Type of protected area (IUCN categories)")+ ylab("Number of species")+
    scale_fill_manual(values=c("olivedrab1","olivedrab3", "olivedrab4", "darkolivegreen", "goldenrod3","goldenrod1", "gold1",
                               "lemonchiffon2","gainsboro", "grey","white" ))+
    scale_y_continuous(expand=c(0,0))+
    theme(legend.position="right", axis.text.x=element_text(angle=30, hjust=1),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
}

#a <- boxplot_template(cover_sr, "current")
#c <- boxplot_template(cover_sr, "ssp126")
#d <- boxplot_template(cover_sr, "ssp370")
#e <- boxplot_template(cover_sr, "ssp585")

cover_matrix <- cover_sr %>% pivot_longer(cols=ssp126_mean:ssp585_mean, names_to = "SSP", values_to = "SR") %>%
  mutate("name"="Species richness")
cover_matrix$SR_change <- cover_matrix$current_mean - cover_matrix$SR

cover_matrix %>% dplyr::select(-x, -y, -name) %>% filter(SSP!="current") %>% group_by(IUCNcat) %>% summarize_all(mean)
cover_matrix %>% dplyr::select(-x, -y, -name)%>% filter(SSP!="current") %>% group_by(IUCNcat) %>% summarize_all(sd)
cover_matrix %>% dplyr::select(-x, -y, -name) %>% filter(SSP!="current") %>% group_by(IUCNcat) %>% summarize_all(max)
cover_matrix %>% dplyr::select(-x, -y, -name) %>% filter(SSP!="current") %>% group_by(IUCNcat) %>% summarize_all(min)


# change in protection in species richness
c <- ggplot(cover_matrix %>% filter(IUCNcat!="Outside.PA" & SSP!="current" & IUCNcat!="Presence" & IUCNcat!="Protected"), aes(x=SSP, y=name))+
  geom_tile(aes(fill=SR_change))+
  scale_fill_gradient2(high="dodgerblue3", low="brown3", mid="white", name="No. species gained (+) or lost(-)        ")+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  facet_wrap(vars(IUCNcat), nrow=1, strip.position="bottom")+
  theme(axis.text.x = element_text(angle=30, hjust=1),
        legend.position="top", legend.text = element_text(size=10), legend.title = element_text(size=5),
        axis.text.y = element_text(size=10), axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(),strip.background = element_blank(),  panel.spacing = unit(0, "lines"), 
        panel.border = element_rect(color="grey"), strip.text = element_text(angle=30, hjust=1))

# calculate change in protection
cover_matrix <- cover_df %>% full_join(cover_df %>% filter(SSP=="current") %>% dplyr::select(-SSP), by=c("SpeciesID", "IUCNcat"),
                                       suffix = c("", "_current"))
cover_matrix$coverage_change <- cover_matrix$coverage_current - cover_matrix$coverage

#cover_matrix <- cover_matrix %>% pivot_wider(id_cols=c(SpeciesID, SSP), names_from=IUCNcat, values_from=coverage_change)

#temp_matrix <- cover_matrix %>% filter(SSP=="ssp126") %>% dplyr::select(-SSP) %>% as.data.frame()
#rownames(temp_matrix) <- temp_matrix$SpeciesID

# plot heatmap-like matrix to show change
cover_matrix$IUCNcat <- factor(cover_matrix$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))

f <- ggplot(cover_matrix %>% filter(IUCNcat!="Outside.PA" & SSP!="current" & IUCNcat!="Presence" & IUCNcat!="Protected"), aes(x=SSP, y=reorder(SpeciesID, desc(SpeciesID))))+
  geom_tile(aes(fill=coverage_change*100))+
  scale_fill_gradient2(high="dodgerblue3", low="brown3", mid="white", name="Change in coverage [%]    ")+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  facet_wrap(vars(IUCNcat), nrow=1, strip.position="top")+
  theme(axis.text.x = element_text(angle=30, hjust=1),
        legend.position="top", legend.text = element_text(size=10), legend.title = element_text(size=15),
        axis.text.y = element_text(size=10), axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(),strip.background = element_blank(),  panel.spacing = unit(0, "lines"), 
        panel.border = element_rect(color="grey"), strip.text = element_blank())

require(gridExtra)
png(paste0(here::here(), "/figures/ProtectionStatus_heatmap_", Taxon_name, ".png"), height=1000, width=1000)
grid.arrange(a, b, c, f, 
             layout_matrix = rbind(c(1,1,1,3,3),
                                   c(1,1,1,3,3),
                                   #c(1,1,1,4,4),
                                   c(2,2,2,4,4),
                                   c(2,2,2,4,4),
                                   c(2,2,2,4,4),
                                   c(2,2,2,4,4),
                                   c(2,2,2,4,4),
                                   c(2,2,2,4,4)))
dev.off()


#- - - - - - - - - - - - - - - - - - - - 
## Protection: some numbers ####
#- - - - - - - - - - - - - - - - - - - -

# protected area per species
cover_df %>% filter(SpeciesID!="_Mean") %>% dplyr::select(-SpeciesID) %>% filter(IUCNcat=="Protected") %>% arrange(coverage)
cover_df %>% filter(SpeciesID=="_Mean") %>% dplyr::select(-SpeciesID) %>% arrange(coverage_km2)

# categories Ia and Ib coverage
cover_df %>% filter(SpeciesID!="_Mean" & SSP=="current") %>% dplyr::select(-SpeciesID) %>% group_by(IUCNcat) %>% summarize_all(mean)
cover_df %>% filter(SpeciesID!="_Mean" & SSP=="current") %>% dplyr::select(-SpeciesID) %>% group_by(IUCNcat) %>% summarize_all(sd)
cover_df %>% filter(SpeciesID!="_Mean") %>% dplyr::select(-SpeciesID) %>% group_by(IUCNcat) %>% summarize_all(sd)
max(cover_df[cover_df$IUCNcat=="Ia",]$coverage)
max(cover_df[cover_df$IUCNcat=="Ib",]$coverage)

cover_sr %>% group_by(IUCNcat) %>% summarize_all(mean)
cover_sr %>% group_by(IUCNcat) %>% summarize_all(sd)

cover_matrix %>% dplyr::select(-SpeciesID) %>% group_by(IUCNcat) %>% summarize_all(mean)
cover_matrix %>% dplyr::select(-SpeciesID) %>% group_by(IUCNcat) %>% summarize_all(sd)

cover_matrix %>% filter(IUCNcat!="Presence") %>% dplyr::select(-SpeciesID, -IUCNcat) %>% summarize_all(mean)
cover_matrix %>% filter(IUCNcat!="Presence") %>% dplyr::select(-SpeciesID, -IUCNcat) %>% summarize_all(sd)
cover_matrix %>% filter(IUCNcat!="Presence") %>% dplyr::select(-SpeciesID, -IUCNcat) %>% summarize_all(median)

# save mean coverage in km2 per IUCN protected area type
write.csv(cover_df %>% 
            group_by(SpeciesID, IUCNcat) %>% 
            dplyr::select(SpeciesID, IUCNcat, coverage_km2) %>% 
            pivot_wider(names_from=IUCNcat, values_from=coverage_km2), 
          file=paste0(here::here(), "/results/ProtectionStatus_coveragePerCategory_", Taxon_name, ".csv"), row.names=F)




#- - - - - - - - - - - - - - - - - - - - 
## Map: BD high, PA high and BD low, PA high etc. ####
#- - - - - - - - - - - - - - - - - - - - 

protect_df$Protected <- round(rowSums(protect_df %>% dplyr::select(-x,-y,-Not.Reported, -Not.Assigned, -Not.Applicable), na.rm=T),2)
protect_df[protect_df$Protected>=1 & !is.na(protect_df$Protected), "Protected"] <- 1

biplot_df <- right_join(protect_df %>% dplyr::select(x,y,Protected), 
                        species_stack %>% dplyr::select(x,y,Richness))

# area (of grid) that is currently protected
sum(biplot_df$Protected)
sum(biplot_df$Protected)*5
sum(biplot_df$Protected) / nrow(biplot_df) * 100

biplot_df <- biplot_df %>% 
  full_join(data.frame("Richness"=c(NA,0:19),
                       "Richness_cont"=c(0,0,rep(1,4), rep(5,5), rep(10,5), rep(15,5)))) %>%
  mutate("Earthworm_richness"=as.factor(Richness_cont)) %>%
  filter(Richness>0)

biplot_df$Protection <- "0"
biplot_df[biplot_df$Protected>0 & !is.na(biplot_df$Protected), "Protection"] <- "1" 
biplot_df[biplot_df$Protected>=0.5 & !is.na(biplot_df$Protected), "Protection"] <- "2"
biplot_df[biplot_df$Protected>=0.75 & !is.na(biplot_df$Protected), "Protection"] <- "3"

biplot_df$Protection <- factor(biplot_df$Protection, levels=c(0,1,2,3))
biplot_df$Earthworm_richness <- factor(biplot_df$Earthworm_richness, levels=c(1,5,10,15))

head(biplot_df)

# define fill scale
biplot_input <- biscale::bi_class(biplot_df, x=Protection, y=Earthworm_richness, dim=4)

# load map
world.inp <- map_data("world")

biplot <- ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), col="grey90", fill="white") +
  xlim(-10, 30) +
  ylim(35, 70) +
  
  geom_tile(data=biplot_input, aes(x=x, y=y, fill=bi_class))+
  bi_scale_fill(pal = "BlueGold", dim=4, flip_axes = T)+
  theme_gray()+
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1),
        panel.background = element_rect(fill="#e8e8e8"))


legend <- bi_legend(pal = "BlueGold",
                    dim = 4,
                    xlab = "Protection",
                    ylab = "Earthworm diversity",
                    size = 12,
                    flip_axes = T)
#legend

finalPlot <- cowplot::ggdraw(biplot) +
  cowplot::draw_plot(legend, 0.03, 0.8, 0.28, 0.15)

png(paste0(here::here(), "/figures/Protection_vs_richness_", Taxon_name, ".png"), width=1500, height=1300)
finalPlot +
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "A", size = 16/.pt,
           fontface = "bold")
dev.off()


#- - - - - - - - - - - - - - - - - - - -
### DRAFT Plots ####
#- - - - - - - - - - - - - - - - - - - -
# 
# # bar chart of percent area covered by PA per species
# f <- ggplot(data=cover_df %>% filter(IUCNcat!="Presence" & IUCNcat!="Unprotected" & IUCNcat!="Outside.PA" & IUCNcat!="Protected"), 
#             aes(y=coverage, x=SSP, fill=IUCNcat))+
#   geom_bar(position="stack", stat="identity")+
#   theme_bw()+
#   xlab("Species")+ ylab("Proportion of range covered by protected area network")+
#   coord_flip()+
#   facet_wrap(vars(SpeciesID))+
#   #scale_fill_viridis_d()+
#   scale_fill_manual(values=c("darkgoldenrod4","darkgoldenrod3", "darkgoldenrod2", "darkgoldenrod1", "goldenrod2","goldenrod1", "gold1",
#                              "lightgoldenrod1","palegoldenrod", "lemonchiffon2","gainsboro" ))+
#   scale_y_continuous(expand=c(0,0), limits=c(0,0.65))+
#   #scale_x_discrete(labels=SpeciesID)+
#   geom_vline(xintercept=1.5, lty=2)+
#   theme(legend.position="none")



