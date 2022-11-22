#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Sensitivity analysis              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# Sensitivity analysis: We will check, if SDMs for species with many occurrence 
# records (i.e., present in >100 grid cells) are similar to the SDMs for the 
# same species made with only a subset of available records (i.e., with only 50, 
# 10, and 5 records).

#setwd("D:/_students/Romy/SoilBiodiversity")

library(here)
library(tidyverse)
library(doParallel)
library(biomod2)

# change temporary directory for files
raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - - -
## Create background data ####
#- - - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"

# load species list containing information on number of occurrence records
speciesNames <- read.csv(file=paste0(here::here(), "/results/Species_list_", Taxon_name, ".csv"))

# load environmental space
Env_clip <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km_clipped.grd"))

# as dataframe
load(paste0(here::here(),"/results/EnvPredictor_2km_df_clipped.RData")) #Env_clip_df

# response variable (i.e., species occurrences) in wide format
occ_points <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))
occ_points <- occ_points %>% rename("x"="ï..x")

# get species with more than equal to 200 records:
if(is.null(speciesNames$NumCells_2km_biomod)) print("Please use the species list in the results folder!")
speciesSub <- unique(speciesNames[speciesNames$NumCells_2km_biomod >= 100,]$SpeciesID)

# subset species' records
occ_points <- occ_points[,c("x", "y", "year", speciesSub)]

# covariates in order of importance (top 10 important)
covarsNames <- c("MAT", "Dist_Coast", "MAP_Seas", "Elev", "Agriculture",
                 "pH", "MAP", "Clay.Silt", "CEC","P" )

## parallelize
# Calculate the number of cores
#no.cores <- length(speciesSub)/2; no.cores
no.cores <- 10

# Initiate cluster used in foreach function
doParallel::registerDoParallel(no.cores)

#- - - - - - - - - - - - - - - - - - - - - - -
## Create data subsets ####
data_sens <- occ_points %>% dplyr::select(x,y)

for(no_replicate in c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10")){ 
  for(no_subset in c(50, 75, 90)){ # 50,
    for(spID in speciesSub) {
      temp_records <- occ_points[!is.na(occ_points[,spID]),] %>%
        dplyr::select(x,y,spID) %>%
        dplyr::slice_sample(prop = no_subset/100) %>%
        filter(!is.na(spID)) 
      temp_records[,paste0(spID, "_", no_subset, "_", no_replicate)] <- temp_records[,spID]
      data_sens <- data_sens %>% 
        full_join(temp_records %>% dplyr::select(-spID), 
                  by=c("x", "y"))
    }
  }
}

write.csv(data_sens, paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_data/SensAna_records.csv"),
          row.names = FALSE)

data_sens <- read.csv(paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_data/SensAna_records.csv"))

speciesSub <- colnames(data_sens %>% dplyr::select(-x, -y))
#speciesSub <- speciesSub[str_detect(speciesSub, "50")]

#- - - - - - - - - - - - - - - - - - - - - - - 
## For loop through all selected species ####
foreach(spID = speciesSub, .export = c("data_sens"), 
        .packages = c("biomod2", "tidyverse")) %dopar% {
   
  setwd(paste0(here::here(), "/results/_Sensitivity_percent"))         
          
  # subset occurrence records
  myResp <- data_sens[!is.na(data_sens[,spID]),]
 
  myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = as.numeric(myResp[,spID]),
                                                expl.var = Env_clip,
                                                resp.xy = myResp[,c("x","y")],
                                                resp.name = spID,
                                                PA.nb.rep = 1,
                                                PA.nb.absences = 10000,
                                                PA.strategy = "random")
  
  # save data
  save(myBiomodData, file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_data/BiomodData_", spID, ".RData"))

  rm(myBiomodData, myResp, myRespCoord, na.id)

 rm(spID)

}

# at the end
doParallel::stopImplicitCluster()

#- - - - - - - - - - - - - - - - - - - - -
## Build models ####
# define parameters of the algorithms
myBiomodOption <- BIOMOD_ModelingOptions(
  GLM = list (type = "quadratic",
              interaction.level = 0,
              myFormula = NULL,
              test = "AIC",
              family = binomial(link = "logit") ),
  
  GAM = list (algo = "GAM_mgcv",
              myFormula = NULL,
              type = "s_smoother",
              interaction.level = 0,
              family =  binomial(link = "logit"),
              method = "GCV.Cp",
              optimizer = c("outer","newton"),
              select = FALSE,
              knots = NULL,
              paraPen = NULL,
              k = -1 ), 		#avoid error messages
  
  MARS = list(myFormula = NULL,
              nk = NULL, 		# maximum number of model terms, NULL: max(21, 2*nb_expl_var+1)
              penalty = 2, 	# default
              thresh = 0.001, 	# default
              nprune = 1+length(covarsNames), # max. number of terms including intercept
              pmethod = "backward" ), #pruning method
  
  MAXENT.Phillips = list(path_to_maxent.jar = paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_sdm/"), # change it to maxent directory
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
  
  GBM = list( distribution = "bernoulli",
              n.trees = 2500,	# default
              interaction.depth = 7, # default
              n.minobsinnode = 5, # default
              shrinkage = 0.001, # default, learning rate
              bag.fraction = 0.5, # default, proportion of observations used in selecting variables
              train.fraction = 0.8, # default 1, train.fraction * nrows(data) observations are used to fit the gbm 
              cv.folds = 10,	# default 3
              keep.data = FALSE, # default
              verbose = FALSE,	# default
              perf.method = "cv", # default
              n.cores = 1 ),	# default
  
  CTA = list(	method = "class", # default, response is factor
              parms = "default", # default
              cost = NULL ),	# default
  
  ANN = list(	NbCV = 10, 		# default, number CV
              size = NULL, 	# default, size parameter will be optimised by cross validation based on model AUC
              decay = NULL, 	# default, decay parameter will be optimised by cross validation
              rang = 0.1, 	# default, initial random weights on [-rang, rang] 
              maxit = 200 ), 	# default, maximum number of iterations
  
  SRE = list(quant = 0.025),	# default
  
  FDA = list(	method = "mars",	# default, regression method used in optimal scaling
              add_args = NULL ),# default
  
  RF = list(	do.classif = TRUE, # default classification random.forest computed, else regression random.forest 
             ntree = 500,	# default
             mtry = 10,		# number of variables randomly sampled as candidates at each split
             nodesize = 1,	# default 5, but 1 for classification, minimum size of terminal nodes
             maxnodes = NULL ) # default, maximum number of terminal nodes trees in the forest
)

# models to predict with
mymodels <- c("GLM","GBM","GAM","CTA","ANN", "SRE", "FDA","MARS","RF","MAXENT.Phillips")
 
## function to get Pseudo-absence dataset
get_PAtab <- function(bfd){
  dplyr::bind_cols(
    x = bfd@coord[, 1],
    y = bfd@coord[, 2],
    status = bfd@data.species,
    bfd@PA
  )
}

#- - - - - - - - - - - - - - - - - - - - - -
# Calculate the number of cores
no.cores <- 10

# Initiate cluster used in foreach function
doParallel::registerDoParallel(no.cores)

# for loop
foreach(spID = speciesSub, 
        .export = c("Taxon_name", "covarsNames"), 
        .packages = c("tidyverse")) %dopar% {
         
          try({
            
            # load background data (pseudo-absences) for each modeling approach
            load(paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_data/BiomodData_", spID, ".RData")) #myBiomodData


		myBiomodData@data.env.var <- myBiomodData@data.env.var[,colnames(myBiomodData@data.env.var) %in% covarsNames]
 
	    # define weights of presence records based on sampling year
 	    temp_weights <- occ_points %>% dplyr::select(x, y, year, substr(spID, 1, 10)) %>% unique()
 	    temp_weights <- temp_weights[!is.na(temp_weights[,4]),]
	    temp_weights <- get_PAtab(myBiomodData) %>% left_join(temp_weights, by=c("x","y"))
	    temp_weights$weight <- 0.1
	    try(temp_weights[!is.na(temp_weights$status),]$weight <- 0.2) #includes NA in year
	    try(temp_weights[!is.na(temp_weights$status) & temp_weights$year >= 1980 & !is.na(temp_weights$year),]$weight <- 0.3)
	    try(temp_weights[!is.na(temp_weights$status) & temp_weights$year >= 1990 & !is.na(temp_weights$year),]$weight <- 0.4)
	    try(temp_weights[!is.na(temp_weights$status) & temp_weights$year >= 2000 & !is.na(temp_weights$year),]$weight <- 0.5)
	    try(temp_weights[!is.na(temp_weights$status) & temp_weights$year >= 2010 & !is.na(temp_weights$year),]$weight <- 0.6)
	    try(temp_weights[!is.na(temp_weights$status) & temp_weights$year >= 2020 & !is.na(temp_weights$year),]$weight <- 0.7) 
         
          # model fitting
          #tmp <- proc.time()[3]
          setwd(paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_sdm/"))
          
          set.seed(32639)
          myBiomodModelOut <- biomod2::BIOMOD_Modeling(myBiomodData,
                                                       models = mymodels,
                                                       models.options = myBiomodOption,
                                                       NbRunEval = 10,   # 3-fold crossvalidation evaluation run
                                                       DataSplit = 80, # use subset of the data for training
                                                       Yweights = temp_weights$weight, # weight to observations, here based on year
                                                       models.eval.meth = "TSS",
                                                       SaveObj = TRUE, #save output on hard drive?
                                                       rescal.all.models = FALSE, #scale all predictions with binomial GLM?
                                                       do.full.models = FALSE, # do evaluation & calibration with whole dataset
                                                       modeling.id = paste(spID))
          
          # ensemble modeling using mean probability
          myBiomodEM <- biomod2::BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                                         chosen.models = "all",  # all algorithms
                                                         em.by = "all",    #evaluated over evaluation data if given (it is not, see Prepare_input.R)
                                                         # note: evaluation not that important as we will calculate measures on independent data
                                                         eval.metric = "TSS", # 'all' would takes same as above in BIOMOD_Modelling
                                                         eval.metric.quality.threshold = NULL, # since some species's auc are naturally low
                                                         prob.mean = FALSE, #estimate mean probabilities across predictions
                                                         prob.cv = TRUE,   #estimate coefficient of variation across predictions
                                                         prob.ci = FALSE,  #estimate confidence interval around the prob.mean
                                                         prob.median = FALSE, #estimate the median of probabilities
                                                         committee.averaging = TRUE, #estimate committee averaging across predictions
                                                         prob.mean.weight = TRUE, #estimate weighted sum of predictions
                                                         prob.mean.weight.decay = "proportional", #the better a model (performance), the higher weight
                                                         VarImport = 0)    #number of permutations to estimate variable importance
          #temp_model_time <- proc.time()[3] - tmp
          
          setwd(here::here())

          }) # end of try loop
          
          print(paste0("The model was now built for ", Taxon_name, ": ", spID, "."))
          
        }

# at the end
stopImplicitCluster()


#- - - - - - - - - - - - - - - - -
## Predict in current climate at 5km ####

# load environmental variables (for projections)
Env_clip <- raster::stack(paste0(here::here(), "/results/EnvPredictor_5km_clipped.grd"))
#Env_clip <- stack(Env_clip)

# as dataframe
load(paste0(here::here(),"/results/EnvPredictor_5km_df_clipped.RData")) #Env_clip_df

registerDoParallel(3)
foreach(spID = speciesSub,
        .export = c("Env_clip", "Env_clip_df"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({   
          
          # list files in species-specific BIOMOD folder
          temp_files <- list.files(paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_sdm/", stringr::str_replace_all(spID, "_", ".")), full.names = TRUE)
          
  	    setwd(paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_sdm/"))

          # load model output
          myBiomodModelOut <- temp_files[stringr::str_detect(temp_files,"[:digit:].models.out")]
          print(myBiomodModelOut)
          myBiomodModelOut <-get(load(myBiomodModelOut))
          
          # load ensemble model output
          myBiomodEM <- temp_files[stringr::str_detect(temp_files,"ensemble.models.out")]
          myBiomodEM <- get(load(myBiomodEM))
          
          tmp <- proc.time()[3]
          ## NOTE: because biomod output can hardly be stored in list file, we will do calculations based on model output now
          # project single models (also needed for ensemble model)
          myBiomodProj <- biomod2::BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                     new.env = Env_clip_df[,colnames(Env_clip_df) %in% covarsNames],        #column/variable names have to perfectly match with training
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
          
          # Get model evaluation values for later
          myBiomodModelEval <- as.data.frame(biomod2::get_evaluations(myBiomodEM))
            
          # save predictions as raster file
          temp_prediction <- myBiomodEnProj@proj@val[,2]
          temp_prediction <- as.numeric(temp_prediction)
          # add names of grid cell (only for those that have no NA in any layer)
          names(temp_prediction) <- rownames(Env_clip_df)
          temp_prediction <- as.data.frame(temp_prediction)
          temp_prediction$x <- Env_clip_df$x
          temp_prediction$y <- Env_clip_df$y
          temp_prediction <- temp_prediction %>% full_join(Env_clip_df %>% dplyr::select(x,y)) %>%
            rename("layer" = temp_prediction)
          temp_prediction$layer <- temp_prediction$layer / 1000
          
          temp_runs <- 1
          
          biomod_list <- list(time_predict=temp_predict_time, validation=myBiomodModelEval, prediction=temp_prediction)
          save(biomod_list, file=paste0("../_SensAna_output/SDM_biomod_", spID, ".RData"))
          
          rm(biomod_list, temp_predict_time, temp_runs, temp_prediction, myBiomodEnProj, myBiomodProj, myBiomodModelEval, myBiomodModelOut, myBiomodEM)
          
          setwd(here::here())

})}
stopImplicitCluster()


#- - - - - - - - - - - - - - - - - - - - - -
## Create maps and calculate richness ####
#- - - - - - - - - - - - - - - - - - - - - -
# load stack created with >=100 records
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
colnames(species_stack) <- substr(colnames(species_stack), 1, 10)

species_stack$no_subset <- 100

# for loop through all species
for(spID in speciesSub){ try({

  ## Load probability maps 
  load(file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_biomod_", spID, ".RData")) #biomod_list
  best_pred <- biomod_list$prediction
  
  print(paste0(spID, " successfully loaded."))
  
  ## Transform to binary maps ####
  
  # extract threshold to define presence/absence: TSS [row 2]
  temp_thresh <- biomod_list$validation[2,str_detect(colnames(biomod_list$validation), "EMcaByTSS_mergedAlgo_mergedRun_mergedData.Cutoff")]/1000
  if(is.na(temp_thresh)) temp_tresh <- 0.9
  
  # change to binary
  best_pred[best_pred$layer>=temp_thresh & !is.na(best_pred$layer), "layer"] <- 1
  best_pred[best_pred$layer<temp_thresh & !is.na(best_pred$layer), "layer"] <- 0
  
  best_pred[,substr(spID, 1,10)] <- best_pred$layer
  best_pred <- best_pred[,c("x","y",substr(spID, 1,10))]%>% mutate("no_subset"=as.numeric(str_split(spID, "_")[[1]][3]))
  
  # save binary
  save(best_pred, file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_bestPrediction_binary_", Taxon_name, "_", spID, ".RData"))
  
  print(paste0("Saved binary prediction of ", spID))
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Stack species binary maps ####

  # add species dataframe to stacked dataframe
  species_stack <- species_stack %>% full_join(best_pred)
  
  print(paste0("Added binary prediction of ", spID, " to the species stack"))
  
  rm(temp_thresh, best_pred, biomod_list)
}, silent=T)} 

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate richness ####
species_stack$Richness <- rowSums(species_stack %>% dplyr::select(Allol_chlo:Satch_mamm), na.rm=T)

#species_stack <- species_stack %>% group_by(x,y, no_subset) %>% summarize_all(sum, na.rm=T)

#- - - - - - - - - - - - - - - - - - - - - -
## Save species stack ####
species_stack_full <- species_stack

species_stack <- species_stack_full %>% filter(no_subset==50)
save(species_stack, file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_stack_binary_", Taxon_name, "_50.RData"))

species_stack <- species_stack_full %>% filter(no_subset==75)
save(species_stack, file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_stack_binary_", Taxon_name, "_75.RData"))

species_stack <- species_stack_full %>% filter(no_subset==90)
save(species_stack, file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_stack_binary_", Taxon_name, "_90.RData"))

species_stack <- species_stack_full %>% filter(no_subset==100)
save(species_stack, file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_stack_binary_", Taxon_name, "_100.RData"))

rm(species_stack_full, species_stack)

#- - - - - - - - - - - - - - - - - - - - - -
## Uncertainty ####

# load environmental space as data frame
load(paste0(here::here(),"/results/EnvPredictor_5km_df_clipped.RData")) #Env_clip_df

for(no_subset in c(50, 75, 90)){
  uncertain_df <- Env_clip_df %>% dplyr::select(x, y)
  
  for(spID in speciesSub[str_detect(speciesSub, as.character(no_subset))]){try({
    
    print(paste0("Species: ", spID))
    
    # list files in species-specific BIOMOD folder
    temp_files <- list.files(paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_sdm/", stringr::str_replace_all(spID, "_", "."), "/proj_modeling"), full.names = TRUE)
    
    #setwd(paste0(here::here(), "/results/", Taxon_name))
    
    myBiomodEnProj <- get(load(temp_files[stringr::str_detect(temp_files,"ensemble.RData")]))
    
    # save predictions as raster file
    temp_prediction <- myBiomodEnProj[,1]
    temp_prediction <- as.numeric(temp_prediction)
    # add names of grid cell (only for those that have no NA in any layer)
    names(temp_prediction) <- rownames(Env_clip_df)
    temp_prediction <- as.data.frame(temp_prediction)
    temp_prediction$x <- Env_clip_df$x
    temp_prediction$y <- Env_clip_df$y
    temp_prediction <- temp_prediction %>% full_join(Env_clip_df %>% dplyr::select(x,y)) %>%
      rename("layer" = temp_prediction)
    temp_prediction$layer <- temp_prediction$layer / 1000
    temp_prediction[,substr(spID, 1, 10)] <- temp_prediction$layer
    
    # add layer to stack
    uncertain_df <- full_join(uncertain_df, temp_prediction %>% dplyr::select(x,y, substr(spID, 1, 10)))
  
  })}
  
  uncertain_df$Mean <- rowMeans(uncertain_df %>% dplyr::select(-x, -y), na.rm=T)
  
  # calculate sd of predictions
  uncertain_df$SD <- apply(uncertain_df %>% dplyr::select(-x, -y, -Mean), 1, sd, na.rm = TRUE)
  
  head(uncertain_df)
  
  # save species' uncertainty map
  save(uncertain_df, file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_Uncertainty_", Taxon_name, "_", no_subset,  ".RData"))
}

## Plotting

for(no_subset in c(50, 75, 90)){
  load(file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_Uncertainty_", Taxon_name,"_", no_subset,  ".RData")) #uncertain_df
  
  # view uncertainty in map 
  world.inp <- map_data("world")
  
  png(file=paste0(here::here(), "/figures/Uncertainty_", Taxon_name,"_", no_subset, ".png"), width=1000, height=1000)
  print({ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-10, 30) +
    ylim(35, 70) +
    
    geom_tile(data=uncertain_df %>% filter(Mean!=0), aes(x=x, y=y, fill=Mean))+
    ggtitle("Coefficient of variation averaged across SDMs", subtitle=paste(no_subset, "% of records"))+
    scale_fill_viridis_c(option="E")+
    theme_bw()+  
    
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position ="bottom",legend.direction = "horizontal")})
  dev.off()
  
  # extract area with uncertainty lower than threshold
  summary(uncertain_df$Mean)
  
  extent_df <- uncertain_df %>% filter(Mean<0.1 & !is.na(Mean)) %>% dplyr::select(x,y)
  save(extent_df, file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_Uncertainty_extent_", Taxon_name, "_", no_subset, ".RData"))
}

#- - - - - - - - - - - - - - - - - - - - - -
## View individual binary maps and species stack ####

# species richness
world.inp <- map_data("world")

for(no_subset in c(50, 75, 90)){ try({
  # load uncertainty extent
  load(file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_Uncertainty_extent_", Taxon_name, "_", no_subset, ".RData")) #extent_df
  
	load(file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_stack_binary_", Taxon_name, "_", no_subset, ".RData")) #species_stack
	
	png(file=paste0(here::here(), "/figures/SensAna_SpeciesRichness_", Taxon_name, "_", no_subset, ".png"), width=1000, height=1000)
	print({ggplot()+
	  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
	    xlim(-10, 30) +
	    ylim(35, 70) +

 	 geom_tile(data=extent_df %>% inner_join(species_stack %>% filter (Richness>0), by=c("x","y")), 
			aes(x=x, y=y, fill=Richness))+
	 ggtitle(paste0(no_subset, "% of records"))+
 	 scale_fill_viridis_c()+
	    geom_tile(data=extent_df %>% inner_join(species_stack %>% filter(Richness==0), by=c("x","y")), aes(x=x, y=y), fill="grey60")+
 	 theme_bw()+
	   theme(axis.title = element_blank(), legend.title = element_blank(), legend.text = element_text(size=10),
	          legend.position = c(0.1,0.9), legend.direction = "horizontal")})
	dev.off()
})}

while (!is.null(dev.list()))  dev.off()


# map binary species distributions
for(no_subset in c(50,75, 90)){
load(file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_stack_binary_", Taxon_name, "_", no_subset, ".RData")) #species_stack
species_stack <- species_stack %>% dplyr::select(-no_subset)
species_stack <- extent_df %>% inner_join(species_stack, by=c("x","y"))
  plots <- lapply(3:(ncol(species_stack)-1), function(s) {try({
  print(s-2)
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-10, 30) +
    ylim(35, 70) +
    
    geom_tile(data=species_stack[!is.na(species_stack[,s]),], 
              aes(x=x, y=y, fill=factor(as.numeric(unlist(species_stack[,s])), levels=c("0", "1", "NA"))))+
    ggtitle(colnames(species_stack)[s])+
    scale_fill_manual(values=c("1"="#440154","0"="grey60","NA"="lightgrey"))+
    theme_bw()+
    guides(fill = guide_legend(# title.hjust = 1, # adjust title if needed
      label.position = "bottom",
      label.hjust = 0.5))+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.9), legend.direction = "horizontal",
          legend.text = element_text(size=20))
})
})

require(gridExtra)
#pdf(file=paste0(here::here(), "/figures/SensAna_DistributionMap_bestBinary_", Taxon_name, "_", no_subset, ".pdf"))
print(png(file=paste0(here::here(), "/figures/SensAna_DistributionMap_bestBinary_", Taxon_name, "_", no_subset, ".png"),width=3000, height=3000))
do.call(grid.arrange, plots)
dev.off()
}

while (!is.null(dev.list()))  dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Difference between no_subset results ####

# calculate difference between subsets and full model output

# check how strong no_subset influences lm output
full_stack <- get(load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData"))) #species_stack
full_stack <- full_stack %>% dplyr::select(x,y,Richness)
full_stack$subset <- 100

for(no_subset in c(50, 75, 90)){		
	load(file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_stack_binary_", Taxon_name, "_", no_subset, ".RData")) #species_stack
	
	species_stack$subset <- no_subset

	full_stack <- full_stack %>% full_join(species_stack)
}
colnames(full_stack)
nrow(full_stack)

load(paste0(here::here(),"/results/EnvPredictor_5km_df_clipped.RData")) #Env_clip_df

data_stack <- full_stack %>% full_join(Env_clip_df)

lm1 <- lm(data=data_stack, Richness~subset+MAT+Dist_Coast+MAP_Seas+CEC+Elev+P+Pop_Dens+Agriculture+pH+Clay.Silt)
summary(lm1)

lm_varImp <- data.frame("t_value"=summary(lm1)[["coefficients"]][,"t value"])
lm_varImp$Predictor <- rownames(lm_varImp)
lm_varImp <- lm_varImp %>% filter(Predictor != "(Intercept)")
lm_varImp$t_abs <- abs(lm_varImp$t_value)
lm_varImp$Direction <- factor(sign(lm_varImp$t_value), 1:(-1), c("positive", "neutral", "negative"))

# transform to long format and add variable categories
pred_tab <- readr::read_csv(file=paste0(here::here(), "/doc/Env_Predictors_table.csv"))

lm_varImp <- lm_varImp%>%
  left_join(pred_tab %>% dplyr::select(Predictor, Category), by="Predictor")

# add category for clay.silt and subset
lm_varImp[lm_varImp$Predictor=="Clay.Silt","Category"] <- "Soil"

plotTopVI <- lm_varImp %>% dplyr::select(t_abs, Predictor, Category, Direction) %>% arrange(desc(t_abs)) %>%
  ggplot(aes(x=reorder(Predictor, t_abs), y=t_abs, fill=Category)) + 
  geom_segment(aes(x=reorder(Predictor, t_abs), xend=reorder(Predictor, t_abs), y=0, yend=t_abs, lty=Direction), color="black") +
  geom_point(aes(color=Category), size=4, alpha=1) +
  coord_flip() +
  xlab("Predictors")+ylab("Variable importance (SR)")+
  theme_bw()+theme(aspect.ratio=1/1)
plotTopVI

png(paste0(here::here(), "/figures/SensAna2_VariableImportance_biomod_top10_lm_", Taxon_name, ".png")); plotTopVI; dev.off()

# save model summary
sink(paste0(here::here(), "/results/SensAna2_Summary_lm1_Crassiclitellata_varImp.txt"))
print(summary(lm1))
sink()

#- - - - - - - - - - - - - - - - - - - - - -
## Sensitive species ####

# calculate species ranges (area)
no_subset <- 75
load(file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_stack_binary_", Taxon_name, "_", no_subset, ".RData")) #species_stack
species_stack <- species_stack %>% dplyr::select(-no_subset)
range_df <- data.frame("SpeciesID"=colnames(species_stack), 
                       "cells"=colSums(species_stack, na.rm=T), 
                       "area_km2"=colSums(species_stack, na.rm=T)*5) %>% 
  filter(SpeciesID!="x", SpeciesID!="y", SpeciesID!="Richness")
rownames(range_df) <- NULL
range_df$subset <- no_subset
head(range_df)

no_subset <- 50
load(file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_stack_binary_", Taxon_name, "_", no_subset, ".RData")) #species_stack
species_stack <- species_stack %>% dplyr::select(-no_subset)
temp_df <- data.frame("SpeciesID"=colnames(species_stack), 
                      "cells"=colSums(species_stack, na.rm=T), 
                      "area_km2"=colSums(species_stack, na.rm=T)*5) %>% 
  filter(SpeciesID!="x", SpeciesID!="y", SpeciesID!="Richness")
rownames(temp_df) <- NULL
temp_df$subset <- no_subset
head(temp_df)

range_df <- range_df %>% full_join(temp_df)

no_subset <- 90
load(file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_stack_binary_", Taxon_name, "_", no_subset, ".RData")) #species_stack
species_stack <- species_stack %>% dplyr::select(-no_subset)
temp_df <- data.frame("SpeciesID"=colnames(species_stack), 
                       "cells"=colSums(species_stack, na.rm=T), 
                       "area_km2"=colSums(species_stack, na.rm=T)*5) %>% 
  filter(SpeciesID!="x", SpeciesID!="y", SpeciesID!="Richness")
rownames(temp_df) <- NULL
temp_df$subset <- no_subset
head(temp_df)

range_df <- range_df %>% full_join(temp_df)
head(range_df)

# add full ranges
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
range_df_current <- data.frame("SpeciesID"=substr(colnames(species_stack), 1, 10), 
                               "cells_full"=colSums(species_stack, na.rm=T), 
                               "area_km2_full"=colSums(species_stack, na.rm=T)*5) %>%
  filter(SpeciesID!="x" & SpeciesID!="y")
rownames(range_df_current) <- NULL
head(range_df_current)  

range_sum <- range_df %>% full_join(range_df_current %>% filter(SpeciesID!="Richness"))
range_sum

range_sum$area_km2_change <- range_sum$area_km2 - range_sum$area_km2_full
range_sum$area_km2_change_p <- range_sum$area_km2_change / range_sum$area_km2_full

range_sum %>% left_join(speciesNames %>% dplyr::select(SpeciesID, NumCells_2km, Species_final)) %>% arrange(subset, area_km2_change_p)

range_sum %>% dplyr::select(-subset) %>% group_by(SpeciesID) %>% summarize_all(mean) %>% arrange(area_km2_change_p)
range_sum %>% dplyr::select(-subset) %>% group_by(SpeciesID) %>% summarize_all(sd)

range_sum %>% arrange(subset, area_km2_change_p)
mean(range_sum$area_km2_change_p); sd(range_sum$area_km2_change_p)

write.csv(range_sum, file=paste0(here::here(), "/results/_Sensitivity_percent/Range_shift_", Taxon_name, ".csv"), row.names=F)
range_sum <- read.csv(file=paste0(here::here(), "/results/_Sensitivity_percent/Range_shift_", Taxon_name, ".csv"))

#- - - - - - - - - - - - - - - - - - - - - -
## Sensitive area ####

full_stack$Richness <- rowSums(full_stack %>% dplyr::select(Allol_chlo:Satch_mamm), na.rm=T)

# full records
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack

full_stack <- full_stack %>% full_join(species_stack %>% mutate("Richness_full"=Richness) %>% dplyr::select(x,y,Richness_full))

# Calculate percent change in distribution
full_stack$Change <- (full_stack$Richness - full_stack$Richness_full)
full_stack$Change_f <- cut(full_stack$Change, 
                              breaks=c(-20, -10, -5, 0, 5, 10, 20),
                              labels=c("[-20,-10]", "[-10,-5]", "[-5,0]", "[0,5]", "[5,10]", "[10,20]"))


save(full_stack, file=paste0(here::here(), "/results/_Sensitivity_percent/SDM_stack_future_richness_change_", Taxon_name, ".RData"))
load(file=paste0(here::here(), "/results/_Sensitivity_percent/SDM_stack_future_richness_change_", Taxon_name, ".RData")) #full_stack

full_stack <- full_stack %>% filter(subset<100)


for(no_subset in c(50, 75, 90)){
  
  # load uncertainty extent
  load(file=paste0(here::here(), "/results/_Sensitivity_percent/_SensAna_output/SDM_Uncertainty_extent_", Taxon_name, "_", no_subset, ".RData")) #extent_df
  
  # plot change in distribution
  png(file=paste0(here::here(), "/figures/SensAna_SpeciesRichness_cert0.1_change_", Taxon_name, "_", no_subset, ".png"),width=1000, height=1000)
  print(ggplot()+
          geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
          xlim(-10, 30) +
          ylim(35, 70) +
          
          geom_tile(data=extent_df %>% inner_join(full_stack %>% filter(!is.na(Change)) %>% filter(Richness!=0 & Richness_full!=0)), 
                    aes(x=x, y=y, fill=Change_f))+
          ggtitle(paste0("Change in species richness (number of species)"), subtitle=no_subset)+
          scale_fill_manual(breaks=c("[10,20]", "[5,10]", "[0,5]", "[-5,0]", "[-10,-5]", "[-20,-10]"), 
                            values=c("steelblue4", "steelblue2", "lightblue","darksalmon", "brown2", "brown4"))+
          geom_tile(data=extent_df %>% inner_join(full_stack %>% filter(Change==0)), aes(x=x, y=y), fill="linen")+
          theme_bw()+
          theme(axis.title = element_blank(), legend.title = element_blank(),
                legend.position = c(0.1,0.95)))
  dev.off()

}





