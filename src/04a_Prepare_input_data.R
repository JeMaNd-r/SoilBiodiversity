#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare input data for SDMs          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)
library(raster)

library(biomod2) # also to create pseudo-absences

library(parallel)
library(doParallel)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0("./results/Species_list_", Taxon_name, ".csv"))
speciesSub <- speciesNames %>% filter(NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique() %>% c()
#speciesSub <- speciesNames %>% filter(family == "Lumbricidae" & NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique()
speciesSub <- c(speciesSub$SpeciesID)

#- - - - - - - - - - - - - - - - - - - - -

# load environmental variables
Env_clip <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km_clipped.grd"))

#- - - - - - - - - - - - - - - - - - - - -
## Prepare data ####
mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))

registerDoParallel(3)
foreach(spID = speciesSub, 
        .combine = rbind,
        .export = c("Env_clip", "mySpeciesOcc"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({
          
          myResp <- mySpeciesOcc[!is.na(mySpeciesOcc[,spID]), c("x","y",spID)]
          
          myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = as.numeric(myResp[,spID]),
                                                        expl.var = Env_clip,
                                                        resp.xy = myResp[,c("x", "y")],
                                                        resp.name = spID,
                                                        PA.nb.rep = 1,
                                                        PA.nb.absences = 10000,
                                                        PA.strategy = "random")
          
          # save data
          save(myBiomodData, file=paste0(here::here(), "/intermediates/BIOMOD_data/BiomodData_", Taxon_name,"_", spID, ".RData"))
          
          rm(myBiomodData, myResp, spID)
        })}

stopImplicitCluster()

#- - - - - - - - - - - - - - - - - - - - -
## Calculate number of records per species ####
records <- data.frame("x"=12, "y"=12,"occ"=1, "SpeciesID"="species")[0,]

for(spID in speciesSub){
  
  print(paste0(spID, " will be added."))
  
  # load biomod data
  load(file=paste0(here::here(), "/intermediates/BIOMOD_data/BiomodData_", Taxon_name,"_", spID, ".RData")) #myBiomodData

  # extract occurrences & pseudo-absences
  myData <- cbind(myBiomodData@data.species, myBiomodData@coord, myBiomodData@data.env.var)
  myData$SpeciesID <- spID
  myData <- myData %>% rename("occ" = "myBiomodData@data.species")
  myData[is.na(myData$occ),"occ"] <- 0
  
  records <- rbind(records, myData[,c("x", "y","occ", "SpeciesID")])
  
}

head(records)
nrow(records) #522,637
nrow(records %>% filter(occ==1)) # 22,637
nrow(records %>% filter(occ==0)) #500,000

records_species <- records %>% group_by(SpeciesID) %>% summarize(across("occ", sum)) %>%
  full_join(records %>% filter(occ==0) %>% group_by(SpeciesID) %>% count(name="Pseudoabsences"))
records_species

records_species %>% filter(occ>=10) %>% count() #41 species
records_species %>% filter(occ>=100) %>% count() #19 species

write.csv(records, file=paste0(here::here(), "/results/Occurrence_rasterized_2km_BIOMOD_", Taxon_name, ".csv"))

#- - - - - - - - - - - - - - - - - - - - -
## Add number of records to Species List table
speciesNames$NumCells_2km_biomod <- 0

for(spID in unique(speciesNames$SpeciesID)){ try({
  speciesNames[speciesNames$SpeciesID==spID,"NumCells_2km_biomod"] <- records_species[records_species$SpeciesID==spID, "occ"]
}, silent=T)}

write.csv(speciesNames, file=paste0(here::here(), "/results/Species_list_", Taxon_name, ".csv"), row.names = F)

#- - - - - - - - - - - - - - - - - - - - -
## Check how many species can be included
count_cuts <- c(0,1,5,10,15,20,50,100,200,300,500,1000)
data <- data.frame("NumOcc" = count_cuts, "NumSpecies"=NA, "NumSpeciesID"=NA, "OccType"="NumCells_2km")
for(i in count_cuts){
  data[data$NumOcc==i,"NumSpecies"] <- speciesNames %>% filter(NumCells_2km>=i) %>% dplyr::select("Species_final") %>% unique() %>% count()
  data[data$NumOcc==i,"NumSpeciesID"] <- speciesNames %>% filter(NumCells_2km>=i) %>% dplyr::select("SpeciesID") %>% unique() %>% count()
}
data

data2 <- data.frame("NumOcc" = count_cuts, "NumSpecies"=NA, "NumSpeciesID"=NA, "OccType"="Records")
for(i in count_cuts){
  data2[data2$NumOcc==i,"NumSpecies"] <- speciesNames %>% filter(Records>=i) %>% dplyr::select("Species_final") %>% unique() %>% count()
  data2[data2$NumOcc==i,"NumSpeciesID"] <- speciesNames %>% filter(Records>=i) %>% dplyr::select("SpeciesID") %>% unique() %>% count()
}
data2

data3 <- data.frame("NumOcc" = count_cuts, "NumSpecies"=NA, "NumSpeciesID"=NA, "OccType"="NumCells_2km_biomod")
for(i in count_cuts){
  data3[data3$NumOcc==i,"NumSpecies"] <- speciesNames %>% filter(NumCells_2km_biomod>=i) %>% dplyr::select("Species_final") %>% unique() %>% count()
  data3[data3$NumOcc==i,"NumSpeciesID"] <- speciesNames %>% filter(NumCells_2km_biomod>=i) %>% dplyr::select("SpeciesID") %>% unique() %>% count()
}
data3

# merge all datasets
data <- rbind(data, data2)
data <- rbind(data, data3)
data

write.csv(data, file=paste0(here::here(), "/results/NumSpecies_perOcc.csv"), row.names=F)


