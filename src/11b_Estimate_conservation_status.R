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
	  add_row("IUCNcat"="Unprotected", "sumCell"=temp_cover[temp_cover$IUCNcat=="Presence","sumCell"]-sum(colSums(temp_df %>% dplyr::select(-Presence), na.rm=T)))
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

write.csv(cover_df, file=paste0(here::here(), "/results/ProtectionStatus_", Taxon_name, ".csv"))


## Plotting ####
cover_df <- read.csv(file=paste0(here::here(), "/results/ProtectionStatus_", Taxon_name, ".csv"))
cover_df$IUCNcat <- factor(cover_df$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected"))

# bar chart of percent area covered by PA per species
png(paste0(here::here(), "/figures/ProtectionStatus_", Taxon_name, ".png"))
ggplot(data=cover_df %>% filter(IUCNcat!="Presence"), aes(y=coverage_km2, x=forcats::fct_reorder(.f=SpeciesID, .x=coverage_km2, .fun = mean), fill=IUCNcat))+
	geom_bar(position="stack", stat="identity")+
	theme_bw()+
  coord_flip()+
  scale_fill_viridis_d()+
  geom_vline(xintercept=9.5, lty=2)+ geom_vline(xintercept=10.5, lty=2)+
	theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1))
dev.off()

# bar chart of percent area covered by PA overall
png(paste0(here::here(), "/figures/ProtectionStatus_total_", Taxon_name, ".png"))
ggplot(data=cover_df %>% filter(IUCNcat!="Presence" & SpeciesID!="_Mean"), aes(x=reorder(IUCNcat, desc(IUCNcat)), y=coverage_km2, fill=IUCNcat))+
  geom_boxplot()+
  geom_jitter(alpha=0.6, height=0.2)+
  theme_bw()+ coord_flip()+
  scale_fill_viridis_d()+
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
dev.off()


#- - - - - - - - - - - - - - - - - - - - 
## Calculate area BD high, PA high and BD low, PA high etc. ####






#- - - - - - - - - - - - - - - - - - - - 
## OLD CODE: other solution... 
#area_protected <- raster::crop(protect_stack[1], species_stack[3][species_stack[3]==1])

  