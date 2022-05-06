#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Estimate conservation status         #
#          author: Romy Zeiss               #
#            date: 2022-05-04               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load stack of IUCN category coverage
protect_stack <- raster::stack(paste0(here::here(), "/data/Shapefiles/WDPA_WDOECM_Dec2021_Public_EU_shp/WDPA_WDOECM_IUCNcat.grd"))

# load species names
speciesNames <- read.csv(file=paste0(here::here(), "/results/Species_list_", Taxon_name, ".csv")) #number of records added

# load species distributions
load(paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
head(species_stack)

# transform protected cover into dataframe
protect_df <- as.data.frame(raster::rasterToPoints(protect_stack))
head(protect_df)

# merge protected and species stack
df <- dplyr::full_join(protect_df, species_stack, by=c("x", "y"))
head(df)

# create empty dataframe
cover_df <- data.frame("IUCNcat" = "I", "sumCell"=1, "SpeciesID"="species", "coverage"=1)[0,]

# calculate percent of coverage per species and IUCN category
for(sp in unique(speciesNames$SpeciesID)){ try({
	temp_df <- df[,c(names(protect_stack), colnames(df)[stringr::str_detect(colnames(df), sp)])]
	temp_df$Presence <- temp_df[,colnames(temp_df)[stringr::str_detect(colnames(temp_df), sp)]]
	temp_df <- temp_df[,c(names(protect_stack), "Presence")]

	if(length(unique(stringr::str_detect(colnames(df), sp)))==1) next

	# keep only presence rows
	temp_df <- temp_df[temp_df[,"Presence"]==1 & !is.na(temp_df[,"Presence"]),]

	# calculate sum of all columns (will give you coverage)
	temp_cover <- data.frame("IUCNcat" = names(temp_df), "sumCell"= as.numeric(colSums(temp_df)))
	temp_cover$SpeciesID <- sp
	temp_cover$coverage <- round(temp_cover$sumCell / sum(temp_df[,"Presence"]),4)

	cover_df <- rbind(cover_df, temp_cover)
	rm(temp_cover, temp_df)
}, silent=TRUE)}

cover_df$coverage_km2 <- round(cover_df$sumCell * 4, 2)

cover_df

write.csv(cover_df, file=paste0(here::here(), "/results/ProtectionStatus_", Taxon_name, ".csv"))


## Plotting ####
cover_df$IUCNcat <- factor(cover_df$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported"))

# bar chart of percent area covered by PA per species
pdf(paste0(here::here(), "/figures/ProtectionStatus_", Taxon_name, ".pdf"))
ggplot(data=cover_df, aes(x=IUCNcat, y=coverage_km2, group=SpeciesID, fill=IUCNcat))+
	geom_bar(stat="identity")+
	facet_wrap(vars(SpeciesID))
dev.off()


#- - - - - - - - - - - - - - - - - - - - 
## OLD CODE: other solution... 
#area_protected <- raster::crop(protect_stack[1], species_stack[3][species_stack[3]==1])

  