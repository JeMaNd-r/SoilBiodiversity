#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("I:/eie/==PERSONAL/RZ_SoilBON/DELETE_SoilBiodiversity_RStudio")
#setwd("D:/_students/Romy/SoilBiodiversity")

#packageurl <- "https://cran.r-project.org/src/contrib/Archive/biomod2/biomod2_3.1-64.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
library(biomod2) 

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8g")) # expand Java memory
gc()
library(tidyverse)
library(here)

# read functions for ensemble of small models (ESM)
#devtools::source_url("https://raw.githubusercontent.com/cran/ecospat/master/R/ecospat.ESM.R")
library(ecospat)
#source("ecospat_function_fixed.R")

library(raster)

library(usdm) # for variable inflation factor calculation

# for model performance:
library(precrec)
library(ggplot2) # for plotting the curves
# devtools::install_github("meeliskull/prg/R_package/prg")
library(prg)
library(ggpubr)
library (ROCR)
library(sdm) # to calculate kappa
library(parallel)
library(doParallel)

# plotting
library(GGally) #for correlations with ggpairs
library(gridExtra)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0("./results/Species_list_", Taxon_name, ".csv"))
#speciesSub <- speciesNames %>% filter(NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique() %>% c()
speciesSub <- speciesNames %>% filter(family == "Lumbricidae" & NumCells_2km >=10) %>% filter(NumCells_2km < 100) %>%
  dplyr::select(SpeciesID) %>% unique()
speciesSub <- c(speciesSub$SpeciesID)

# covariates in order of importance (top 10 important)
covarsNames <- c("Dist_Coast", "Forest_Deci", "Pop_Dens", "Shrubland", 
                 "Forest_Coni", "Agriculture", "Pasture", "MAP_Seas", "Clay.Silt", "CEC")

# define future scenarios
scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
                               "mri-esm2-0", "ukesm1-0-ll"), "_",
                             rep(c("ssp126", "ssp370", "ssp585"),5)))

# load data with sampling year information
occ_points <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))
str(occ_points)

#- - - - - - - - - - - - - - - - - - - - -
# note: we will load the datasets before each individual model

# load environmental variables (for projections)
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km_normalized.grd"))
#Env_norm <- stack(Env_norm)

# as dataframe
load(paste0(here::here(),"/results/EnvPredictor_2km_df_normalized.RData")) #Env_norm_df

# define formula for GLM (and biomod)
form <- paste0("occ ~ ", paste0(paste0("s(", covarsNames, ")"), collapse=" + "))

# Calculate the number of cores
no.cores <-  parallel::detectCores()/2 

## function to get Pseudo-absence dataset
get_PAtab <- function(bfd){
  dplyr::bind_cols(
    x = bfd@coord[, 1],
    y = bfd@coord[, 2],
    status = bfd@data.species,
    bfd@PA
  )
}

#- - - - - - - - - - - - - - - - - - - - -
## Prepare data ####
mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))

registerDoParallel(no.cores)
foreach(spID = speciesSub, 
        .export = c("Env_norm", "Env_norm_df", "form", "mySpeciesOcc"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({
          
          # subset coordinates with presence of species spID
          coords <- mySpeciesOcc %>% filter(!is.na(mySpeciesOcc[,spID])) %>% dplyr::select(x,y)
          
          # transform to SpatialPoints object
          myResp <- sp::SpatialPoints(coords, proj4string=CRS(projection(Env_norm)))
          
          # transform into BiomodData
          myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                                        expl.var = Env_norm,
                                                        resp.name = spID,
                                                        PA.nb.rep = 1,
                                                        PA.nb.absences = 10000,
                                                        PA.strategy = "random")
          
          # save data
          save(myBiomodData, file=paste0(here::here(), "/results/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData"))
          
          rm(myBiomodData, myResp, myRespCoord, spID, na.id)
        })}
stopImplicitCluster()


registerDoParallel(no.cores)
foreach(spID = speciesSub, 
        .export = c("Env_norm", "Env_norm_df", "form", "mySpeciesOcc"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({
          
          load(paste0(here::here(), "/results/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData")) #myBiomodData
          
          myData <- data.frame(occ = myBiomodData@data.species,
                               x = myBiomodData@coord$x,
                               y = myBiomodData@coord$y)
          myData[is.na(myData$occ),"occ"] <- 0
          
          myESMData <- biomod2::BIOMOD_FormatingData(resp.var = as.numeric(myData$occ),
                                                        expl.var = Env_norm,
                                                        resp.xy = myData %>% dplyr::select(x,y),
                                                        resp.name = spID)
          
          # save data
          save(myESMData, file=paste0(here::here(), "/results/", Taxon_name, "/ESMData_", Taxon_name,"_", spID, ".RData"))
          
          rm(myBiomodData, myData, spID)
        })}
stopImplicitCluster()

#- - - - - - - - - - - - - - - - - - - - -
# if more than 10 but less than 100 occurrences   
myBiomodOption <- BIOMOD_ModelingOptions(
  MAXENT.Phillips = list(path_to_maxent.jar = paste0(here::here(), "/results"), # change it to maxent directory
                         memory_allocated = NULL, # use default from Java defined above
                         visible = FALSE, 	# don't make maxEnt user interface visible
                         linear = TRUE, 	# linear features allowed
                         quadratic = TRUE, # quadratic allowed
                         product = TRUE,	# product allowed
                         threshold = TRUE,	# threshold allowed
                         hinge = TRUE,	# hinge allowed
                         lq2lqptthreshold = 80, # default
                         l2lqthreshold = 10, # default
                         hingethreshold = 15, # default
                         beta_threshold = -1, # default
                         beta_categorical = -1, # default
                         beta_lqp = -1, # default
                         beta_hinge = -1, # default
                         betamultiplier = 1, # default
                         defaultprevalence = 0.5 ), #default   
)

registerDoParallel(no.cores)
foreach(spID = speciesSub, 
        .export = c("Env_norm", "Env_norm_df", "form"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({
          
          load(paste0(here::here(), "/results/", Taxon_name, "/ESMData_", Taxon_name,"_", spID, ".RData")) #myESMData
          
          # subset covarsNames
          myESMData@data.env.var <- myESMData@data.env.var[,colnames(myESMData@data.env.var) %in% covarsNames]
          
          # model fitting
          tmp <- proc.time()[3]
          setwd(paste0("./results"))
          
          set.seed(32639)
          
          myBiomodModelOut <- ecospat::ecospat.ESM.Modeling(data = myESMData,
                                        models = "MAXENT.Phillips",
                                        models.options = myBiomodOption,
                                        Prevalence = NULL,
                                        tune = TRUE, # TRUE: estimate optimal parameters for the models
                                        NbRunEval = 10,   # 3-fold crossvalidation evaluation run
                                        DataSplit = 80, # use subset of the data for training
                                        weighting.score = "TSS",
                                        #ESM_Projection = FALSE, #no projections now (will be done later)
                                        cleanup = 2, #when to delete temporary unused files, in hours
                                        modeling.id = paste(spID,"_Modeling", sep = ""))
          
          
          # ensemble modeling
          myBiomodEM <- ecospat.ESM.EnsembleModeling(ESM.modeling.output = myBiomodModelOut,
                                                     weighting.score = "TSS")
          
          temp_model_time <- proc.time()[3] - tmp
          
        })}
stopImplicitCluster()

#- - - - - - - - - - - - - - - - -
## Predict in current climate at 5km ####

# load environmental variables (for projections)
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_5km_normalized.grd"))
#Env_norm <- stack(Env_norm)

# as dataframe
load(paste0(here::here(),"/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df

registerDoParallel(no.cores)
foreach(spID = speciesSub, 
        .export = c("Env_norm", "Env_norm_df", "form"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({   
          
          # list files in species-specific BIOMOD folder
          temp_files <- list.files(paste0(here::here(), "/results/", Taxon_name, "/", stringr::str_replace(spID, "_", ".")), full.names = TRUE)
          
          setwd(paste0(here::here(), "/results/", Taxon_name))
          
          # load model output
          myBiomodModelOut <- temp_files[stringr::str_detect(temp_files,"Modeling.models.out")]
          print(myBiomodModelOut)
          myBiomodModelOut <-get(load(myBiomodModelOut))
          
          # load ensemble model output
          myBiomodEM <- temp_files[stringr::str_detect(temp_files,"Modelingensemble.models.out")]
          myBiomodEM <- get(load(myBiomodEM))
          
          tmp <- proc.time()[3]
          ## NOTE: because biomod output can hardly be stored in list file, we will do calculations based on model output now
          # project single models (also needed for ensemble model)
          myBiomodProj <- biomod2::BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                     new.env = Env_norm_df[,colnames(Env_norm_df) %in% covarsNames],        #column/variable names have to perfectly match with training
                                                     proj.name = "modeling",  #name of the new folder being created
                                                     selected.models = "all", #use all models
                                                     binary.meth = NULL,     #binary transformation according to criteria, or no transformation if NULL
                                                     compress = TRUE,         #compression format of objects stored on hard drive
                                                     build.clamping.mask = TRUE, #TRUE: clamping mask will be saved on hard drive different
                                                     do.stack = TRUE,         #save output projections as rasterstack (if not too heavy)
                                                     output.format = ".RData", #what format should projections have: RData, grd or img
                                                     keep.in.memory = TRUE)  #FALSE: only story link to copy to projection file
          
          # project ensemble of all models
          myBiomodEnProj <- biomod2::BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                                                EM.output = myBiomodEM,
                                                                #... same arguments as above could be added but are not necessary when loading myBiomodProj
                                                                selected.models = "all")
          
          temp_predict_time <- proc.time()[3] - tmp
          
          # extracting the values for ensemble validation
          myEnProjDF <- as.data.frame(get_predictions(myBiomodEM)[,2]) #for weighted probability mean
          
          # see the first few validations
          # note: the prediction scale of biomod is between 0 and 1000
          #head(myEnProjDF)
          
          # Get model evaluation values for later
          myBiomodModelEval <- as.data.frame(biomod2::get_evaluations(myBiomodEM))
          
          # Calculate variable importance across all PA sets, eval runs and algorithms
          # and extract only the one for weighed mean predictions (for later)
          temp_varImp <- biomod2::get_variables_importance(myBiomodEM)[, , 2]
          # average across 3 runs
          temp_varImp <- temp_varImp %>% as.data.frame() %>%
            mutate(mean_vi = as.numeric(rowMeans(temp_varImp, na.rm=T)),
                   Predictor = rownames(temp_varImp))
          colnames(temp_varImp)[colnames(temp_varImp) == "mean_vi"] <- "biomod"
          temp_varImp <- temp_varImp[,c("biomod", "Predictor")]
          
          # save predictions as raster file
          temp_prediction <- myBiomodEnProj@proj@val[,2]
          temp_prediction <- as.numeric(temp_prediction)
          # add names of grid cell (only for those that have no NA in any layer)
          names(temp_prediction) <- rownames(Env_norm_df)
          temp_prediction <- as.data.frame(temp_prediction)
          temp_prediction$x <- Env_norm_df$x
          temp_prediction$y <- Env_norm_df$y
          temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
            rename("layer" = temp_prediction)
          temp_prediction$layer <- temp_prediction$layer / 1000
          
          temp_runs <- 1
          
          biomod_list <- list(time_predict=temp_predict_time, validation=myBiomodModelEval, prediction=temp_prediction, varImp=temp_varImp, evaluation=myBiomodModelEval)
          save(biomod_list, file=paste0("./SDMs/SDM_biomod_", spID, ".RData"))
          
          rm(biomod_list, temp_predict_time, temp_runs, temp_prediction, temp_varImp, myBiomodEnProj, myBiomodProj, myBiomodModelEval, myEnProjDF, myBiomodModelOut, myBiomodEM)
          
          setwd(here::here())
          
        })}
stopImplicitCluster()


#- - - - - - - - - - - - - - - - -
## Predict in future climate ####

setwd(here::here())

registerDoParallel(no.cores)
foreach(spID = speciesSub, 
        .export = c("Env_norm", "Env_norm_df", "form"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({ 
          
          # list files in species-specific BIOMOD folder
          temp_files <- list.files(paste0(here::here(), "/results/", Taxon_name, "/", stringr::str_replace(spID, "_", ".")), full.names = TRUE)
          
          # load model output
          myBiomodModelOut <- temp_files[stringr::str_detect(temp_files,"Modeling.models.out")]
          print(myBiomodModelOut)
          myBiomodModelOut <-get(load(myBiomodModelOut))
          
          # load ensemble model output
          myBiomodEM <- temp_files[stringr::str_detect(temp_files,"Modelingensemble.models.out")]
          myBiomodEM <- get(load(myBiomodEM))
          
          for(no_future in scenarioNames){
            
            # load env. data with future climate (MAP, MAT, MAP_Seas)
            load(paste0(here::here(), "/results/_FutureEnvironment/EnvPredictor_2041-2070_", no_future, "_5km_df_normalized.RData")) #temp_Env_df
            
            # one loop per future climate subset, one with both future, each one with only 1 future and 1 current climate
            for(subclim in c("TP", "T", "P")){
              
              if(subclim=="TP"){
                temp_Env_sub <- temp_Env_df[,c("x", "y", colnames(temp_Env_df)[colnames(temp_Env_df) %in% covarsNames])]
              }
              
              if(subclim=="T"){
                temp_Env_sub <- temp_Env_df[,c("x", "y", colnames(temp_Env_df)[colnames(temp_Env_df) %in% covarsNames])]
                temp_Env_sub$MAP_Seas <- Env_norm_df$MAP_Seas
              }
              
              if(subclim=="P"){
                temp_Env_sub <- temp_Env_df[,c("x", "y", colnames(temp_Env_df)[colnames(temp_Env_df) %in% covarsNames])]
                temp_Env_sub$MAT <- Env_norm_df$MAT
              }
              
              setwd(paste0(here::here(), "/results/", Taxon_name))
              
              ## NOTE: because biomod output can hardly be stored in list file, we will do calculations based on model output now
              # project single models (also needed for ensemble model)
              myBiomodProj <- biomod2::BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                         new.env = temp_Env_sub[,colnames(temp_Env_sub) %in% covarsNames],  #column/variable names have to perfectly match with training
                                                         proj.name = "modeling",  #name of the new folder being created
                                                         selected.models = "all", #use all models
                                                         binary.meth = NULL,     #binary transformation according to criteria, or no transformation if NULL
                                                         compress = TRUE,         #compression format of objects stored on hard drive
                                                         build.clamping.mask = TRUE, #TRUE: clamping mask will be saved on hard drive different
                                                         do.stack = TRUE,         #save output projections as rasterstack (if not too heavy)
                                                         output.format = ".RData", #what format should projections have: RData, grd or img
                                                         keep.in.memory = TRUE) #FALSE: : only story link to copy to projection file
              
              # project ensemble of all models
              myBiomodEnProj <- biomod2::BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                                                    EM.output = myBiomodEM,
                                                                    #... same arguments as above could be added but are not necessary when loading myBiomodProj
                                                                    selected.models = "all")
              
              # save predictions as raster file
              temp_prediction <- myBiomodEnProj@proj@val[,2]
              temp_prediction <- as.numeric(temp_prediction)
              # add names of grid cell (only for those that have no NA in any layer)
              names(temp_prediction) <- rownames(temp_Env_sub)
              temp_prediction <- as.data.frame(temp_prediction)
              temp_prediction$x <- temp_Env_sub$x
              temp_prediction$y <- temp_Env_sub$y
              temp_prediction <- temp_prediction %>% full_join(temp_Env_sub %>% dplyr::select(x,y)) %>%
                rename("layer" = temp_prediction)
              temp_prediction$layer <- temp_prediction$layer / 1000
              
              setwd(here::here())
              save(temp_prediction, file=paste0(here::here(), "/results/_Maps/SDM_2041-2070_", no_future, "_", subclim, "_biomod_", spID,  ".RData")) 
              rm(temp_prediction, temp_Env_sub, myBiomodEnProj, myBiomodProj)
            }
          }
          
        })}
stopImplicitCluster()



#- - - - - - - - - - - - - - - - - - - - - -
## Create maps and calculate richness ####
#- - - - - - - - - - - - - - - - - - - - - -
# create empty data frame
species_stack <- Env_norm_df %>% dplyr::select(x, y)

# for loop through all species
for(spID in speciesSub){ try({
  
  ## Load probability maps 
  load(file=paste0(here::here(), "/SDMs/SDM_biomod_", spID, ".RData")) #biomod_list
  
  
  print(paste0(spID, " successfully loaded."))
  
  ## Transform to binary maps ####
  
  # extract threshold to define presence/absence
  temp_thresh <- mod_eval[mod_eval$model==temp_model,"thres.maxTSS"]
  if(is.na(temp_thresh)) temp_tresh <- 0.9
  
  # change to binary
  best_pred[best_pred$layer>=temp_thresh & !is.na(best_pred$layer), "layer"] <- 1
  best_pred[best_pred$layer<temp_thresh & !is.na(best_pred$layer), "layer"] <- 0
  
  best_pred[,paste0(spID,"_", temp_model)] <- best_pred$layer
  best_pred <- best_pred[,c("x","y",paste0(spID,"_", temp_model))]
  
  # save binary
  save(best_pred, file=paste0(here::here(), "/results/_Maps/SDM_bestPrediction_binary_", Taxon_name, "_", spID, ".RData"))
  
  print(paste0("Saved binary prediction of ", spID))
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Stack species binary maps ####
  
  # add species dataframe to stacked dataframe
  species_stack <- species_stack %>% full_join(best_pred, by=c("x","y"))
  
  print(paste0("Added binary prediction of ", spID, " to the species stack"))
  
  rm(temp_thresh, best_pred, temp_model)
}, silent=T)}  

head(species_stack)


#- - - - - - - - - - - - - - - - - - - - - -
## Calculate richness ####
species_stack$Richness <- rowSums(species_stack %>% dplyr::select(-x, -y), na.rm=T)

#- - - - - - - - - - - - - - - - - - - - - -
## Save species stack ####
save(species_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData"))


#- - - - - - - - - - - - - - - - - - - - - -
## View individual binary maps and species stack ####
# species richness
png(file=paste0(here::here(), "/figures/SpeciesRichness_", Taxon_name, ".png"),width=1000, height=1000)
ggplot(data=species_stack, aes(x=x, y=y, fill=Richness))+
  geom_tile()+
  ggtitle("Species richness (number of species)")+
  scale_fill_viridis_c()+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()

while (!is.null(dev.list()))  dev.off()


# map binary species distributions
plots <- lapply(3:(ncol(species_stack)-1), function(s) {try({
  print(s-2)
  ggplot(data=species_stack, aes(x=x, y=y, fill=as.factor(species_stack[,s])))+
    geom_tile()+
    ggtitle(colnames(species_stack)[s])+
    scale_fill_manual(values=c("1"="#fde725","0"="#440154","NA"="lightgrey"))+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.4))
})
})

require(gridExtra)
#pdf(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_", Taxon_name, ".pdf"))
png(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_", Taxon_name, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots)
dev.off()

while (!is.null(dev.list()))  dev.off()



#- - - - - - - - - - - - - - - - - - - - -
## TRASH?
#- - - - - - - - - - - - - - - - - - - - -
# if less than 100 occurrences (but more than 10)
if(speciesNames[speciesNames$SpeciesID==spID, "NumCells_2km"]<100){ 
  
} 

#- - - - - - - - - - - - - - - - - - - - -
## Check if everything went well ####
#- - - - - - - - - - - - - - - - - - - - -

for(spID in speciesSub){ try({
  
  check_files <- list.files(paste0("./results/", Taxon_name, "/temp_files"))
  check_files <- check_files[stringr::str_detect(check_files, spID)]
  
  print(check_files)
  if(length(check_files)==17){ print("Everything looks well. Please continue :)")
  }else{ print("Some algorithms weren't saved; please check which ones are missing.")}
  
})}

stopImplicitCluster()
