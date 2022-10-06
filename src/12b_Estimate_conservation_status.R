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

cover_sr_current <- protect_df %>% 
  mutate("Unprotected"=ifelse(rowSums(protect_df %>% dplyr::select(-x, -y))==0, 1, 0)) %>%
  full_join(species_stack %>% dplyr::select(x, y, Richness)) %>%
  pivot_longer(cols=II:Unprotected, names_to="IUCNcat", values_to = "PA_coverage") %>%
  inner_join(extent_df) %>%
  filter(PA_coverage>0)

head(cover_sr_current)

## SSP126, 370, 585 ####

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

cover_sr <- cover_sr %>% full_join(cover_sr_current %>% mutate("current_mean"=Richness) %>% dplyr::select(x,y,IUCNcat,current_mean))

save(cover_sr, file=paste0(here::here(), "/results/ProtectionStatus_SR_SSPs_", Taxon_name, ".csv"))


  