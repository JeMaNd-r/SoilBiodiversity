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


## Plotting ####
cover_df <- read.csv(file=paste0(here::here(), "/results/ProtectionStatus_", Taxon_name, ".csv"))
cover_df$IUCNcat <- factor(cover_df$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))

# bar chart of percent area covered by PA per species
png(paste0(here::here(), "/figures/ProtectionStatus_", Taxon_name, ".png"))
ggplot(data=cover_df %>% filter(IUCNcat!="Presence" & IUCNcat!="Unprotected" & IUCNcat!="Protected"), aes(y=coverage_km2, x=forcats::fct_reorder(.f=SpeciesID, .x=coverage_km2, .fun = mean), fill=IUCNcat))+
	geom_bar(position="stack", stat="identity")+
	theme_bw()+
  coord_flip()+
  #scale_fill_viridis_d()+
  scale_fill_manual(values=c("darkgoldenrod4","darkgoldenrod3", "darkgoldenrod2", "darkgoldenrod1", "goldenrod2","goldenrod1", "gold1",
                              "lightgoldenrod1","palegoldenrod", "lemonchiffon2","gainsboro" ))+
  geom_vline(xintercept=9.5, lty=2)+ geom_vline(xintercept=10.5, lty=2)+
	theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1))
dev.off()

# bar chart of percent area covered by PA overall
png(paste0(here::here(), "/figures/ProtectionStatus_total_", Taxon_name, ".png"))
ggplot(data=cover_df %>%  filter(IUCNcat!="Presence" & IUCNcat!="Unprotected" & IUCNcat!="Protected" & SpeciesID!="_Mean"), aes(x=reorder(IUCNcat, desc(IUCNcat)), y=coverage_km2, fill=IUCNcat))+
  geom_boxplot()+
  geom_jitter(alpha=0.6, height=0.2)+
  theme_bw()+ coord_flip()+
  scale_fill_manual(values=c("darkgoldenrod4","darkgoldenrod3", "darkgoldenrod2", "darkgoldenrod1", "goldenrod2","goldenrod1", "gold1",
                             "lightgoldenrod1","palegoldenrod", "lemonchiffon2","gainsboro" ))+
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
dev.off()


#- - - - - - - - - - - - - - - - - - - - 
## Calculate area BD high, PA high and BD low, PA high etc. ####

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

  