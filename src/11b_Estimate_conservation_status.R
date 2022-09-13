#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Estimate conservation status         #
#          author: Romy Zeiss               #
#            date: 2022-05-04               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load stack of IUCN category coverage
load(file=paste0(here::here(), "/data/Shapefiles/WDPA_WDOECM_Dec2021_Public_EU_shp/WDPA_WDOECM_IUCNcat_df.RData")) #protect_df
head(protect_df)

# load species names
speciesNames <- read.csv(file=paste0(here::here(), "/results/Species_list_", Taxon_name, ".csv")) #number of records added

# load species distributions
load(paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
head(species_stack)

# merge protected and species stack
df <- dplyr::full_join(protect_df, species_stack, by=c("x", "y"))
head(df)

# create empty dataframe
cover_df <- data.frame("IUCNcat" = "I", "sumCell"=1, "SpeciesID"="species", "coverage"=1)[0,]

# calculate percent of coverage per species and IUCN category
for(sp in unique(speciesNames$SpeciesID)){ try({
	temp_df <- df[,c(names(protect_df %>% dplyr::select(-x, -y)), colnames(df)[stringr::str_detect(colnames(df), sp)])]
	temp_df$Presence <- temp_df[,colnames(temp_df)[stringr::str_detect(colnames(temp_df), sp)]]
	temp_df <- temp_df[,c(names(protect_df %>% dplyr::select(-x, -y)), "Presence")]

	# keep only presence rows
	temp_df <- temp_df[temp_df[,"Presence"]==1 & !is.na(temp_df[,"Presence"]),]

	# calculate sum of all columns (will give you coverage)
	temp_cover <- data.frame("IUCNcat" = names(temp_df), "sumCell"= as.numeric(colSums(temp_df)))
	temp_cover <- temp_cover %>%
	  add_row("IUCNcat"="Unprotected", 
	          "sumCell"=temp_cover[temp_cover$IUCNcat=="Presence","sumCell"]-sum(colSums(temp_df %>% dplyr::select(Ia, Ib, II, III, IV, V, VI), na.rm=T))) %>%
	  add_row("IUCNcat"="Outside.PA", 
	          "sumCell"=temp_cover[temp_cover$IUCNcat=="Presence","sumCell"]-sum(colSums(temp_df %>% dplyr::select(Ia, Ib, II, III, IV, V, VI, Not.Assigned, Not.Reported, Not.Applicable), na.rm=T))) %>%
	  add_row("IUCNcat"="Protected", 
	          "sumCell"=sum(colSums(temp_df %>% dplyr::select(Ia, Ib, II, III, IV, V, VI), na.rm=T)))
	
	temp_cover$SpeciesID <- sp
	temp_cover$coverage <- round(temp_cover$sumCell / sum(temp_df[,"Presence"], na.rm=T),4)

	cover_df <- rbind(cover_df, temp_cover)
	rm(temp_cover, temp_df)
	
	print(paste0("Species ", sp, " is ready."))
	
}, silent=TRUE)}

cover_df$coverage_km2 <- round(cover_df$sumCell * 5, 2)

cover_df <- cover_df %>% arrange(SpeciesID, IUCNcat) %>% filter(!is.na(coverage))

cover_df <- rbind(cover_df, 
                  cbind(cover_df %>% group_by(IUCNcat) %>% dplyr::select(-SpeciesID) %>% summarize_all(mean), "SpeciesID"="_Mean"))
cover_df

write.csv(cover_df, file=paste0(here::here(), "/results/ProtectionStatus_", Taxon_name, ".csv"), row.names=T)


## Calculate number of species per IUCN category ####
# load uncertainty extent for all maps
load(file=paste0(here::here(), "/results/_Maps/SDM_Uncertainty_extent_", Taxon_name, ".RData")) #extent_df

cover_sr <- protect_df %>% 
  mutate("Unprotected"=ifelse(rowSums(protect_df %>% dplyr::select(-x, -y))==0, 1, 0)) %>%
  full_join(species_stack %>% dplyr::select(x, y, Richness)) %>%
  pivot_longer(cols=II:Unprotected, names_to="IUCNcat", values_to = "PA_coverage") %>%
  inner_join(extent_df) %>%
  filter(PA_coverage>0)

head(cover_sr)

## Plotting ####
cover_df <- read.csv(file=paste0(here::here(), "/results/ProtectionStatus_", Taxon_name, ".csv"))
cover_df$IUCNcat <- factor(cover_df$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))
cover_sr$IUCNcat <- factor(cover_sr$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))

# boxplot of percent area covered by PA overall
a <- ggplot(data=cover_sr %>%  filter(IUCNcat!="Presence" & IUCNcat!="Protected"), aes(x=reorder(IUCNcat, desc(IUCNcat)), y=Richness, fill=IUCNcat))+
  geom_violin(width=1.4, alpha=0.7)+
  # geom_boxplot(width=0.1, color="black", fill="white", alpha=1)+
  stat_summary(fun = "mean",geom = "point",color = "black", size=4)+
  #geom_jitter(alpha=0.6, width=0.2)+
  theme_bw()+ 
  xlab("Type of protected area (IUCN categories)")+ ylab("Number of species")+
  scale_fill_manual(values=c("darkgoldenrod4","darkgoldenrod3", "darkgoldenrod2", "darkgoldenrod1", "goldenrod2","goldenrod1", "gold1",
                             "lightgoldenrod1","palegoldenrod", "lemonchiffon2","gainsboro" ))+
  scale_y_continuous(expand=c(0,0))+
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1))

# bar chart of percent area covered by PA per species
b <- ggplot(data=cover_df %>% filter(IUCNcat!="Presence" & IUCNcat!="Unprotected" & IUCNcat!="Protected"), aes(y=coverage, x=SpeciesID, fill=IUCNcat))+
	geom_bar(position="stack", stat="identity")+
	theme_bw()+
  xlab("Species")+ ylab("Proportion of range covered by protected area network")+
  coord_flip()+
  #scale_fill_viridis_d()+
  scale_fill_manual(values=c("darkgoldenrod4","darkgoldenrod3", "darkgoldenrod2", "darkgoldenrod1", "goldenrod2","goldenrod1", "gold1",
                              "lightgoldenrod1","palegoldenrod", "lemonchiffon2","gainsboro" ))+
  scale_y_continuous(expand=c(0,0))+
  geom_vline(xintercept=1.5, lty=2)+
	theme(legend.position="none")

png(paste0(here::here(), "/figures/ProtectionStatus_current_", Taxon_name, ".png"))
grid.arrange(a, b, heights=c(1,1.5))
dev.off()

## SSP126 ####
# load species distributions
load(paste0(here::here(), "/results/_Maps/SDM_stack_future_species_", Taxon_name, ".RData")) #future_stack
head(future_stack)

# merge protected and species stack
cover_sr <- protect_df %>% 
  mutate("Unprotected"=ifelse(rowSums(protect_df %>% dplyr::select(-x, -y))==0, 1, 0)) %>%
  full_join(future_stack[,c("x", "y", colnames(future_stack)[str_detect(colnames(future_stack), "Richness")])]) %>%
  pivot_longer(cols=II:Unprotected, names_to="IUCNcat", values_to = "PA_coverage") %>%
  inner_join(extent_df) %>%
  filter(PA_coverage>0)

head(cover_sr)

# calculate average per SSP
for(temp_ssp in c("ssp126", "ssp370", "ssp585")){
  temp_cols <- colnames(cover_sr)[stringr::str_detect(colnames(cover_sr), temp_ssp)]
  cover_sr[,as.character(paste0(temp_ssp, "_mean"))] <- rowMeans(cover_sr[,temp_cols])
  cover_sr[,as.character(paste0(temp_ssp, "_max"))] <- matrixStats::rowMaxs(as.matrix(cover_sr[,temp_cols]))
  cover_sr[,as.character(paste0(temp_ssp, "_min"))] <- matrixStats::rowMins(as.matrix(cover_sr[,temp_cols]))
  cover_sr[,as.character(paste0(temp_ssp, "_sd"))] <- matrixStats::rowSds(as.matrix(cover_sr[,temp_cols]))
}
colnames(cover_sr)

cover_sr <- cover_sr %>% dplyr::select(x,y,ssp126_mean, ssp370_mean, ssp585_mean, IUCNcat)

# merge protected and species stack
df <- dplyr::full_join(protect_df, future_stack[,c("x", "y", colnames(future_stack)[str_detect(colnames(future_stack), ".future_mean")])], by=c("x", "y"))
head(df)

# create empty dataframe
cover_df <- data.frame("IUCNcat" = "I", "sumCell"=1, "SpeciesID"="species", "coverage"=1)[0,]

# calculate percent of coverage per species and IUCN category
for(sp in unique(paste0(speciesNames$SpeciesID, ".future_mean"))){ try({
  temp_df <- df[,c(names(protect_df %>% dplyr::select(-x, -y)), colnames(df)[stringr::str_detect(colnames(df), sp)])]
  temp_df$Presence <- temp_df[,colnames(temp_df)[stringr::str_detect(colnames(temp_df), sp)]]
  temp_df <- temp_df[,c(names(protect_df %>% dplyr::select(-x, -y)), "Presence")]
  
  # keep only presence rows
  temp_df <- temp_df[temp_df[,"Presence"]==1 & !is.na(temp_df[,"Presence"]),]
  
  # calculate sum of all columns (will give you coverage)
  temp_cover <- data.frame("IUCNcat" = names(temp_df), "sumCell"= as.numeric(colSums(temp_df)))
  temp_cover <- temp_cover %>%
    add_row("IUCNcat"="Unprotected", 
            "sumCell"=temp_cover[temp_cover$IUCNcat=="Presence","sumCell"]-sum(colSums(temp_df %>% dplyr::select(Ia, Ib, II, III, IV, V, VI), na.rm=T))) %>%
    add_row("IUCNcat"="Outside.PA", 
            "sumCell"=temp_cover[temp_cover$IUCNcat=="Presence","sumCell"]-sum(colSums(temp_df %>% dplyr::select(Ia, Ib, II, III, IV, V, VI, Not.Assigned, Not.Reported, Not.Applicable), na.rm=T))) %>%
    add_row("IUCNcat"="Protected", 
            "sumCell"=sum(colSums(temp_df %>% dplyr::select(Ia, Ib, II, III, IV, V, VI), na.rm=T)))
  
  temp_cover$SpeciesID <- sp
  temp_cover$coverage <- round(temp_cover$sumCell / sum(temp_df[,"Presence"], na.rm=T),4)
  
  cover_df <- rbind(cover_df, temp_cover)
  rm(temp_cover, temp_df)
  
  print(paste0("Species ", sp, " is ready."))
  
}, silent=TRUE)}

cover_df$coverage_km2 <- round(cover_df$sumCell * 5, 2)
cover_df <- cover_df %>% arrange(SpeciesID, IUCNcat) %>% filter(!is.na(coverage))
cover_df <- rbind(cover_df, 
                  cbind(cover_df %>% group_by(IUCNcat) %>% dplyr::select(-SpeciesID) %>% summarize_all(mean), "SpeciesID"="_Mean"))
cover_df

cover_df$IUCNcat <- factor(cover_df$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))
cover_sr$IUCNcat <- factor(cover_sr$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))

# boxplot of percent area covered by PA overall
a <- ggplot(data=cover_sr %>%  filter(IUCNcat!="Presence" & IUCNcat!="Protected"), aes(x=reorder(IUCNcat, desc(IUCNcat)), y=ssp126_mean, fill=IUCNcat))+
  geom_violin(width=1.4, alpha=0.7)+
  # geom_boxplot(width=0.1, color="black", fill="white", alpha=1)+
  stat_summary(fun = "mean",geom = "point",color = "black", size=4)+
  #geom_jitter(alpha=0.6, width=0.2)+
  theme_bw()+ 
  xlab("Type of protected area (IUCN categories)")+ ylab("Number of species")+
  scale_fill_manual(values=c("darkgoldenrod4","darkgoldenrod3", "darkgoldenrod2", "darkgoldenrod1", "goldenrod2","goldenrod1", "gold1",
                             "lightgoldenrod1","palegoldenrod", "lemonchiffon2","gainsboro" ))+
  scale_y_continuous(expand=c(0,0))+
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1))

# bar chart of percent area covered by PA per species
b <- ggplot(data=cover_df %>% filter(IUCNcat!="Presence" & IUCNcat!="Unprotected" & IUCNcat!="Protected"), aes(y=coverage, x=SpeciesID, fill=IUCNcat))+
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  xlab("Species")+ ylab("Proportion of range covered by protected area network")+
  coord_flip()+
  #scale_fill_viridis_d()+
  scale_fill_manual(values=c("darkgoldenrod4","darkgoldenrod3", "darkgoldenrod2", "darkgoldenrod1", "goldenrod2","goldenrod1", "gold1",
                             "lightgoldenrod1","palegoldenrod", "lemonchiffon2","gainsboro" ))+
  scale_y_continuous(expand=c(0,0))+
  geom_vline(xintercept=1.5, lty=2)+
  theme(legend.position="none")

png(paste0(here::here(), "/figures/ProtectionStatus_ssp126_", Taxon_name, ".png"))
grid.arrange(a, b, heights=c(1,1.5))
dev.off()


#- - - - - - - - - - - - - - - - - - - - 
## Map: Calculate area BD high, PA high and BD low, PA high etc. ####

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


## some numbers
# protected area per species
cover_df %>% filter(SpeciesID!="_Mean") %>% dplyr::select(-SpeciesID) %>% filter(IUCNcat=="Protected") %>% arrange(coverage)
cover_df %>% filter(SpeciesID=="_Mean") %>% dplyr::select(-SpeciesID) %>% arrange(coverage_km2)

# categories Ia and Ib coverage
cover_df %>% filter(SpeciesID!="_Mean") %>% dplyr::select(-SpeciesID) %>% group_by(IUCNcat) %>% summarize_all(sd)
max(cover_df[cover_df$IUCNcat=="Ia",]$coverage)
max(cover_df[cover_df$IUCNcat=="Ib",]$coverage)

# save mean coverage in km2 per IUCN protected area type
write.csv(cover_df %>% 
            group_by(SpeciesID, IUCNcat) %>% 
            dplyr::select(SpeciesID, IUCNcat, coverage_km2) %>% 
            pivot_wider(names_from=IUCNcat, values_from=coverage_km2), 
          file=paste0(here::here(), "/results/ProtectionStatus_coveragePerCategory_", Taxon_name, ".csv"), row.names=F)

  