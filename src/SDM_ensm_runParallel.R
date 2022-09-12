#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

setwd("I:/eie/==PERSONAL/RZ_SoilBON/SoilBiodiversity_RStudio")
#setwd("D:/_students/Romy/SoilBiodiversity")

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8g")) # expand Java memory
gc()
library(tidyverse)
library(here)

# read functions for ensemble of small moels (ESM)
#devtools::source_url("https://raw.githubusercontent.com/cran/ecospat/master/R/ecospat.ESM.R")
library(ecospat)
source("ecospat_function_fixed.R")

library(raster)

library(biomod2) # also to create pseudo-absences

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
speciesSub <- speciesNames %>% filter(family == "Lumbricidae" & NumCells_2km >=10) %>% dplyr::select(SpeciesID) %>% unique()
speciesSub <- c(speciesSub$SpeciesID)

# covariates in order of importance (top 10 important)
covarsNames <- c("MAT", "Dist_Coast", "MAP_Seas", "CEC", "Elev",
                 "P", "Pop_Dens", "Agriculture", "pH", "Clay.Silt")

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
          
          myResp <- as.numeric(mySpeciesOcc[,spID])
          
          # get NAs id
          na.id <- which(is.na(myResp))
          # remove NAs to enforce PA sampling to be done on explanatory rasters
          myResp <- myResp[-na.id]
          
          myRespCoord <- mySpeciesOcc[-na.id,c('x','y')]
          
          myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                                        expl.var = Env_norm,
                                                        resp.xy = myRespCoord,
                                                        resp.name = spID,
                                                        PA.nb.rep = 1,
                                                        PA.nb.absences = 10000,
                                                        PA.strategy = "random")
          
          # save data
          save(myBiomodData, file=paste0(here::here(), "/results/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData"))

          rm(myBiomodData, myResp, myRespCoord, spID, na.id)
})}
stopImplicitCluster()

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
         
#- - - - - - - - - - - - - - - - - - - - -
# if more than 100 occurrences   
registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 100,]$SpeciesID), 
        .export = c("Env_norm", "Env_norm_df", "form"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({
          
          load(paste0(here::here(), "/results/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData")) #myBiomodData
          
          # subset covarsNames
          myBiomodData@data.env.var <- myBiomodData@data.env.var[,colnames(myBiomodData@data.env.var) %in% covarsNames]
 
	    # define weights of presence records based on sampling year
 	    temp_weights <- mySpeciesOcc %>% dplyr::select(x, y, year, spID) %>% unique()
 	    temp_weights <- temp_weights[!is.na(temp_weights[,4]),]
	    temp_weights <- get_PAtab(myBiomodData) %>% left_join(temp_weights, by=c("x","y"))
	    temp_weights$weight <- 0.1
	    temp_weights[!is.na(temp_weights$status),]$weight <- 0.2 #includes NA in year
	    temp_weights[!is.na(temp_weights$status) & temp_weights$year >= 1980 & !is.na(temp_weights$year),]$weight <- 0.3
	    temp_weights[!is.na(temp_weights$status) & temp_weights$year >= 1990 & !is.na(temp_weights$year),]$weight <- 0.4
	    temp_weights[!is.na(temp_weights$status) & temp_weights$year >= 2000 & !is.na(temp_weights$year),]$weight <- 0.5
	    temp_weights[!is.na(temp_weights$status) & temp_weights$year >= 2010 & !is.na(temp_weights$year),]$weight <- 0.6
	    temp_weights[!is.na(temp_weights$status) & temp_weights$year >= 2020 & !is.na(temp_weights$year),]$weight <- 0.7 
         
          # model fitting
          #tmp <- proc.time()[3]
          setwd(paste0(here::here(), "/results/", Taxon_name))
          
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
                                                       modeling.id = paste(spID,"_Modeling", sep = ""))
          
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
                                                         VarImport = 5)    #number of permutations to estimate variable importance
          #temp_model_time <- proc.time()[3] - tmp
          
          setwd(here::here())
})}
stopImplicitCluster()

#- - - - - - - - - - - - - - - - -
## Predict in current climate at 5km ####

# load environmental variables (for projections)
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_5km_normalized.grd"))
#Env_norm <- stack(Env_norm)

# as dataframe
load(paste0(here::here(),"/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df

registerDoParallel(3)
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
          # average across 10 runs
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
          
          biomod_list <- list(time_predict=temp_predict_time, validation=myBiomodModelEval, prediction=temp_prediction, varImp=temp_varImp)
          save(biomod_list, file=paste0("./_SDMs/SDM_biomod_", spID, ".RData"))
          
          rm(biomod_list, temp_predict_time, temp_runs, temp_prediction, temp_varImp, myBiomodEnProj, myBiomodProj, myBiomodModelEval, myEnProjDF, myBiomodModelOut, myBiomodEM)
          
          setwd(here::here())

        })}
stopImplicitCluster()


#- - - - - - - - - - - - - - - - -
## Predict in future climate ####

setwd(here::here())
setwd(paste0(here::here(), "/results/", Taxon_name))

# load current environmental variables (for projections)
load(paste0(here::here(),"/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df

no.cores <- 3
registerDoParallel(no.cores)
foreach(spID = speciesSub,
        .export = c("Env_norm", "Env_norm_df", "form", "Taxon_name"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({ 

for(spID in speciesSub){ try({ 

          # list files in species-specific BIOMOD folder
          temp_files <- list.files(paste0(here::here(), "/results/", Taxon_name, "/", stringr::str_replace(spID, "_", ".")), full.names = TRUE)
          
	    temp_files

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


          ## NOTE: because biomod output can hardly be stored in list file, we will do calculations based on model output now
          # project single models (also needed for ensemble model)
          myBiomodProj <- biomod2::BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                     new.env = temp_Env_sub[,colnames(temp_Env_sub) %in% covarsNames],        #column/variable names have to perfectly match with training
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
          
          #temp_predict_time <- proc.time()[3] - tmp
          
          # save predictions as raster file
          temp_prediction <- myBiomodEnProj@proj@val[,2]
          temp_prediction <- as.numeric(temp_prediction)
          # add names of grid cell (only for those that have no NA in any layer)
          names(temp_prediction) <- rownames(temp_Env_sub)
          temp_prediction <- as.data.frame(temp_prediction)
          temp_prediction$x <- temp_Env_sub$x
          temp_prediction$y <-temp_Env_sub$y
          temp_prediction <- temp_prediction %>% full_join(temp_Env_sub %>% dplyr::select(x,y))
          colnames(temp_prediction)[1] <- "layer"
          temp_prediction$layer <- temp_prediction$layer / 1000
          
          save(temp_prediction, file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_2041-2070_", no_future, "_", subclim, "_biomod_", spID,  ".RData")) 
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
  load(file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_biomod_", spID, ".RData")) #biomod_list
  best_pred <- biomod_list$prediction
  
  print(paste0(spID, " successfully loaded."))
  
  ## Transform to binary maps ####
  
  # extract threshold to define presence/absence: TSS [row 2]
  temp_thresh <- biomod_list$validation[2,str_detect(colnames(biomod_list$validation), "EMcaByTSS_mergedAlgo_mergedRun_mergedData.Cutoff")]/1000
  if(is.na(temp_thresh)) temp_tresh <- 0.9
  
  # change to binary
  best_pred[best_pred$layer>=temp_thresh & !is.na(best_pred$layer), "layer"] <- 1
  best_pred[best_pred$layer<temp_thresh & !is.na(best_pred$layer), "layer"] <- 0
  
  best_pred[,paste0(spID,"_current")] <- best_pred$layer
  best_pred <- best_pred[,c("x","y",paste0(spID,"_current"))]
  
  # save binary
  save(best_pred, file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_bestPrediction_binary_", Taxon_name, "_", spID, ".RData"))
  
  print(paste0("Saved binary prediction of ", spID))
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Stack species binary maps ####
  
  # add species dataframe to stacked dataframe
  species_stack <- species_stack %>% full_join(best_pred, by=c("x","y"))
  
  print(paste0("Added binary prediction of ", spID, " to the species stack"))
  
  rm(temp_thresh, best_pred)
}, silent=T)}  

head(species_stack)


#- - - - - - - - - - - - - - - - - - - - - -
## Calculate richness ####
species_stack$Richness <- rowSums(species_stack %>% dplyr::select(-x, -y), na.rm=F)


#- - - - - - - - - - - - - - - - - - - - - -
## Save species stack ####
save(species_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData"))

#load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack

#- - - - - - - - - - - - - - - - - - - - - -
# Calculate area with 19 species
species_stack %>% filter(Richness==19 & !is.na(Richness)) %>% count()
3914*5

#- - - - - - - - - - - - - - - - - - - - - -
## Beta diversity ####
# vegan::vegdist() cannot handle such big vectors
# alternatively, we calculate clusters
temp_matrix <- species_stack[species_stack$Richness>0 & !is.na(species_stack$Richness),colnames(species_stack)[stringr::str_detect("_current", string = colnames(species_stack))]]
temp_matrix <- temp_matrix[complete.cases(temp_matrix),]

temp_clust <- kmeans(temp_matrix, 10) 
temp_clust <- cbind(species_stack[species_stack$Richness>0 & !is.na(species_stack$Richness),],
                   "Cluster"=temp_clust$cluster)

species_stack <- species_stack %>% full_join(temp_clust)


#- - - - - - - - - - - - - - - - - - - - - -
## View individual binary maps and species stack ####

# extract most prominent species
View(as.data.frame(colSums(species_stack, na.rm=T)) %>% arrange(colSums(species_stack, na.rm=T)))

# species richness
world.inp <- map_data("world")

png(file=paste0(here::here(), "/figures/SpeciesRichness_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(-23, 60) +
  ylim(31, 75) +

  geom_tile(data=species_stack %>% filter(Richness>0), aes(x=x, y=y, fill=Richness))+
  ggtitle("Species richness (number of species)")+
  scale_fill_viridis_c()+
  geom_tile(data=species_stack %>% filter(Richness==0), aes(x=x, y=y), fill="grey60")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()

while (!is.null(dev.list()))  dev.off()


# map binary species distributions
plots <- lapply(3:(ncol(species_stack)-1), function(s) {try({
  print(s-2)
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-23, 60) +
    ylim(31, 75) +
    
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
#pdf(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_", Taxon_name, ".pdf"))
png(file=paste0(here::here(), "/figures/DistributionMap_bestBinary_", Taxon_name, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots)
dev.off()

while (!is.null(dev.list()))  dev.off()

## View clusters
world.inp <- map_data("world")

png(file=paste0(here::here(), "/figures/Clusters_K10_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(-23, 60) +
  ylim(31, 75) +
  
  geom_tile(data=species_stack %>% filter(Richness>0), aes(x=x, y=y, fill=as.factor(Cluster)))+
  ggtitle("Clusters (k=10)")+
  scale_fill_viridis_d()+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()


#- - - - - - - - - - - - - - - - - - - - - -
## Create future maps and calculate richness ####
for(no_future in scenarioNames){
  
  # one loop per future climate subset, one with both future, each one with only 1 future and 1 current climate
  for(subclim in c("TP", "T", "P")){
        
    #- - - - - - - - - - - - - - - - - - - - - -
    # create empty data frame
    species_stack <- Env_norm_df %>% dplyr::select(x, y)
    
    # for loop through all species
    for(spID in speciesSub){ try({
      
      ## Load probability maps 
      load(file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_2041-2070_", no_future, "_", subclim, "_biomod_", spID,  ".RData")) #biomod_list
      best_pred <- temp_prediction
      
	# load model information 
  	load(file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_biomod_", spID, ".RData")) #biomod_list
  
      print(paste0(spID, " successfully loaded."))
      
      ## Transform to binary maps ####
      
      # extract threshold to define presence/absence
      temp_thresh <- biomod_list$validation[2,str_detect(colnames(biomod_list$validation), "EMcaByTSS_mergedAlgo_mergedRun_mergedData.Cutoff")]/1000
      if(is.na(temp_thresh)) temp_tresh <- 0.9
      
      # change to binary
      best_pred[best_pred$layer>=temp_thresh & !is.na(best_pred$layer), "layer"] <- 1
      best_pred[best_pred$layer<temp_thresh & !is.na(best_pred$layer), "layer"] <- 0
      
      best_pred[,paste0(spID,"_future")] <- best_pred$layer
      best_pred <- best_pred[,c("x","y",paste0(spID,"_future"))]
      
      # save binary
      save(best_pred, file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_bestPrediction_binary_2041-2070_", no_future, "_", subclim, "_biomod_", spID,  ".RData"))
      
      print(paste0("Saved binary prediction of ", spID))
      
      #- - - - - - - - - - - - - - - - - - - - - -
      ## Stack species binary maps ####
      
      # add species dataframe to stacked dataframe
      species_stack <- species_stack %>% full_join(best_pred, by=c("x","y"))
      
      print(paste0("Added binary prediction of ", spID, " to the species stack"))
      
      rm(temp_thresh, best_pred, temp_prediction)
    }, silent=T)}  
    
    head(species_stack)
    
    
    #- - - - - - - - - - - - - - - - - - - - - -
    ## Calculate richness ####
    species_stack$Richness <- rowSums(species_stack %>% dplyr::select(-x, -y), na.rm=F)
    
    #- - - - - - - - - - - - - - - - - - - - - -
    ## Save species stack ####
    save(species_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", "2041-2070_", no_future, "_", subclim, ".RData"))
    
  }
}


#- - - - - - - - - - - - - - - - - - - - - -
## View individual binary maps and species stack ####
world.inp <- map_data("world")

for(no_future in scenarioNames){
  
  # only plot subclim scenario TP
  subclim <- "TP"
    
  load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_2041-2070_", no_future, "_", subclim, ".RData")) #species_stack
  
  # species richness
  png(file=paste0(here::here(), "/figures/SpeciesRichness_", "2041-2070_", no_future, "_", subclim, "_", Taxon_name, ".png"),width=1000, height=1000)
  print(ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-23, 60) +
    ylim(31, 75) +
    
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
      xlim(-23, 60) +
      ylim(31, 75) +
      
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
  
#- - - - - - - - - - - - - - - - - - - - -
## Species-specific future predictions ####

future_stack <- get(load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", "2041-2070_", scenarioNames[1], "_TP.RData"))) #species_stack

for(no_future in scenarioNames[2:length(scenarioNames)]){
  
  load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", "2041-2070_", no_future, "_TP.RData")) #species_stack

  # add layer to stack
  future_stack <- full_join(future_stack, species_stack, suffix=c("", paste0(".", no_future)), by=c("x", "y"))
  
}
colnames(future_stack)[3:22] <- paste0(colnames(future_stack)[3:22], ".", scenarioNames[1])
colnames(future_stack)

# calculate average future prediction per species
for(spID in unique(speciesNames[speciesNames$NumCells_2km>=100,]$SpeciesID)){ try({
  future_stack[,as.character(paste0(spID, ".future_mean"))] <- rowSums(future_stack[,stringr::str_detect(colnames(future_stack), spID)])
  future_stack[,as.character(paste0(spID, ".future_max"))] <- matrixStats::rowMaxs(as.matrix(future_stack[,stringr::str_detect(colnames(future_stack), spID)]))
  future_stack[,as.character(paste0(spID, ".future_min"))] <- matrixStats::rowMins(as.matrix(future_stack[,stringr::str_detect(colnames(future_stack), spID)]))
})}
colnames(future_stack)

save(future_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_future_species_", Taxon_name, ".RData"))
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_future_species_", Taxon_name, ".RData")) #future_stack

# calculate species ranges (area)
range_df <- data.frame("scenario"=colnames(future_stack), 
                       "cells"=colSums(future_stack, na.rm=T), 
                       "area_km2"=colSums(future_stack, na.rm=T)*5) %>%
  filter(scenario!="x" & scenario!="y") %>%
  tidyr::separate(scenario, c("SpeciesID", "scenario"), "[.]")
rownames(range_df) <- NULL
head(range_df)

range_sum <- range_df[range_df$SpeciesID %in% range_df$SpeciesID[stringr::str_detect(range_df$SpeciesID, "_future")],] %>% group_by(SpeciesID) %>% dplyr::select(-scenario) %>% 
  summarise(across(everything(), list(min = min, max = max, mean = mean, sd = sd)))

# add current ranges
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
range_df_current <- data.frame("scenario"=colnames(species_stack), 
                       "cells"=colSums(species_stack, na.rm=T), 
                       "area_km2"=colSums(species_stack, na.rm=T)*5) %>%
  filter(scenario!="x" & scenario!="y")
rownames(range_df_current) <- NULL
head(range_df_current)  

range_sum <- range_sum %>% cbind(range_df_current %>% filter(scenario!="Richness"))
range_sum

range_sum$area_km2_change <- range_sum$area_km2_mean - range_sum$area_km2
range_sum$area_km2_change_p <- range_sum$area_km2_change / range_sum$area_km2
range_sum %>% arrange(area_km2_change_p)
mean(range_sum$area_km2_change_p); sd(range_sum$area_km2_change_p)

write.csv(range_sum, file=paste0(here::here(), "/results/Range_shift_", Taxon_name, ".csv"), row.names=F)

# plot species range decline
png(file=paste0(here::here(), "/figures/Range_shift_", "2041-2070_", Taxon_name, ".png"))
ggplot(range_sum) +
  geom_segment( aes(x=reorder(SpeciesID, area_km2), xend=reorder(SpeciesID, area_km2), y=area_km2/1000, yend=area_km2_mean/1000), color="grey") +
  geom_point( aes(x=reorder(SpeciesID, area_km2), y=area_km2/1000, color="area_km2"),  size=5)+
  geom_point( aes(x=reorder(SpeciesID, area_km2), y=area_km2_mean/1000, color="area_km2_mean"), size=5)+
  scale_color_manual(values = c("deepskyblue4", "orange"),
                   guide  = guide_legend(), 
                   name   = "Group",
                   labels = c("Current", "Future (mean)")) +
  coord_flip()+
  theme_bw() +
  theme(legend.position = c(0.8, 0.2)) +
  xlab("") +
  ylab("Range size in 1,000 km²")
dev.off()

# plot both current and future range per species in one plot
for(spID in unique(speciesNames[speciesNames$NumCells_2km>=100, "SpeciesID"])){ 
  temp_cols <- colnames(future_stack)[stringr::str_detect(colnames(future_stack), paste0(spID, "_future."))]
  plots <- lapply(temp_cols, function(s) {try({
    print(s)
    ggplot()+
      geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
      xlim(-23, 60) +
      ylim(31, 75) +
      
      geom_tile(data=future_stack[!is.na(future_stack[,s]),], 
                aes(x=x, y=y, fill=as.factor(future_stack[!is.na(future_stack[,s]),s])))+
      ggtitle(s)+
      scale_fill_manual(values=c("1"="#440154","0"="grey","NA"="lightgrey"))+
      theme_bw()+
      theme(axis.title = element_blank(), legend.title = element_blank(),
            legend.position = c(0.1,0.4))
  })})
  require(gridExtra)
  png(file=paste0(here::here(), "/figures/DistributionMap_2041-2070_", spID, ".png"),width=3000, height=3000)
  do.call(grid.arrange, plots)
  dev.off()
}

# plot max and min future distribution
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


#- - - - - - - - - - - - - - - - - - - - -
## Average future predictions ####
world.inp <- map_data("world")
average_stack <- Env_norm_df %>% dplyr::select(x, y)

for(no_future in scenarioNames){
  
  # only plot subclim scenario TP
  subclim <- "TP"
    
  load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", "2041-2070_", no_future, "_", subclim, ".RData")) #species_stack
  species_stack[,as.character(no_future)] <- species_stack$Richness
  species_stack <- species_stack[,c("x","y",no_future)]
      
  # add layer to stack
  average_stack <- full_join(average_stack, species_stack)
 
}

average_stack$Mean <- rowMeans(average_stack %>% dplyr::select(-x, -y), na.rm=T)

# calculate average per SSP
for(temp_ssp in c("ssp126", "ssp370", "ssp585")){
    temp_cols <- colnames(average_stack)[stringr::str_detect(colnames(average_stack), temp_ssp)]
    average_stack[,as.character(paste0(temp_ssp, "_mean"))] <- rowMeans(average_stack[,temp_cols])
    average_stack[,as.character(paste0(temp_ssp, "_max"))] <- matrixStats::rowMaxs(as.matrix(average_stack[,temp_cols]))
    average_stack[,as.character(paste0(temp_ssp, "_min"))] <- matrixStats::rowMins(as.matrix(average_stack[,temp_cols]))
    average_stack[,as.character(paste0(temp_ssp, "_sd"))] <- matrixStats::rowSds(as.matrix(average_stack[,temp_cols]))
}
colnames(average_stack)

# Calculate percent change in distribution
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
average_stack <- average_stack %>% dplyr::rename(FutureRichness=Mean) %>%
				full_join(species_stack %>% dplyr::select(x,y,Richness))

average_stack$Change <- (average_stack$FutureRichness - average_stack$Richness)
average_stack$Change_f <- cut(average_stack$Change, 
					breaks=c(-15, -10, -5, 0, 5, 10),
					labels=c("[-15,-10]", "[-10,-5]", "[-5,0]", "[0,5]", "[5,10]"))

average_stack$Change_ssp126 <- (average_stack$ssp126_mean - average_stack$Richness)
average_stack$Change_f_ssp126 <- cut(average_stack$Change_ssp126, 
                              breaks=c(-15, -10, -5, 0, 5, 10, 15),
                              labels=c("[-15,-10]", "[-10,-5]", "[-5,0]", "[0,5]", "[5,10]", "[10,15]"))
average_stack$Change_ssp370 <- (average_stack$ssp370_mean - average_stack$Richness)
average_stack$Change_f_ssp370 <- cut(average_stack$Change_ssp370, 
                                     breaks=c(-15, -10, -5, 0, 5, 10, 15),
                                     labels=c("[-15,-10]", "[-10,-5]", "[-5,0]", "[0,5]", "[5,10]", "[10,15]"))
average_stack$Change_ssp585 <- (average_stack$ssp585_mean - average_stack$Richness)
average_stack$Change_f_ssp585 <- cut(average_stack$Change_ssp585, 
                                     breaks=c(-15, -10, -5, 0, 5, 10, 15),
                                     labels=c("[-15,-10]", "[-10,-5]", "[-5,0]", "[0,5]", "[5,10]", "[10,15]"))

save(average_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_future_richness_change_", Taxon_name, ".RData"))
load(file=paste0(here::here(), "/results/_Maps/SDM_stack_future_richness_change_", Taxon_name, ".RData")) #average_stack

# plot future mean distribution
png(file=paste0(here::here(), "/figures/SpeciesRichness_", "2041-2070_future_", Taxon_name, ".png"),width=1000, height=1000)
print(ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-23, 60) +
    ylim(31, 75) +
    
    geom_tile(data=average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0), 
		aes(x=x, y=y, fill=FutureRichness))+
    ggtitle(paste0("Future species richness (number of species)"))+
    scale_fill_viridis_c()+
		geom_tile(data=average_stack %>% filter(Richness==0), aes(x=x, y=y), fill="grey60")+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.4)))
dev.off()

# plot change in distribution
png(file=paste0(here::here(), "/figures/SpeciesRichness_", "2041-2070_change_", Taxon_name, ".png"),width=1000, height=1000)
print(ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    xlim(-23, 60) +
    ylim(31, 75) +
    
    geom_tile(data=average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0), 
		aes(x=x, y=y, fill=Change_f))+
    ggtitle(paste0("Change in species richness (number of species)"))+
    scale_fill_viridis_d(breaks=c("[5,10]", "[0,5]", "[-5,0]", "[-10,-5]", "[-15,-10]"), option="B")+
		geom_tile(data=average_stack %>% filter(Change==0), aes(x=x, y=y), fill="lightblue")+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.4)))
dev.off()

# ssp126 plot change in distribution
png(file=paste0(here::here(), "/figures/SpeciesRichness_", "2041-2070_change_ssp126_", Taxon_name, ".png"),width=1000, height=1000)
print(ggplot()+
        geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
        xlim(-23, 60) +
        ylim(31, 75) +
        
        geom_tile(data=average_stack %>% filter(Richness==0), aes(x=x, y=y), fill="grey60")+
        geom_tile(data=average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0), 
                  aes(x=x, y=y, fill=Change_f_ssp126))+
        ggtitle(paste0("Change in species richness (number of species) SSP126"))+
        scale_fill_viridis_d(breaks=c("[10,15]", "[5,10]", "[0,5]", "[-5,0]", "[-10,-5]", "[-15,-10]"), option="B")+
        geom_tile(data=average_stack %>% filter(Change_ssp126==0), aes(x=x, y=y), fill="lightblue")+
        theme_bw()+
        theme(axis.title = element_blank(), legend.title = element_blank(),
              legend.position = c(0.1,0.4)))
dev.off()

# ssp370 plot change in distribution
png(file=paste0(here::here(), "/figures/SpeciesRichness_", "2041-2070_change_ssp370_", Taxon_name, ".png"),width=1000, height=1000)
print(ggplot()+
        geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
        xlim(-23, 60) +
        ylim(31, 75) +
        
        geom_tile(data=average_stack %>% filter(Richness==0), aes(x=x, y=y), fill="grey60")+
        geom_tile(data=average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0), 
                  aes(x=x, y=y, fill=Change_f_ssp370))+
        ggtitle(paste0("Change in species richness (number of species) SSP370"))+
        scale_fill_viridis_d(breaks=c("[10,15]", "[5,10]", "[0,5]", "[-5,0]", "[-10,-5]", "[-15,-10]"), option="B")+
        geom_tile(data=average_stack %>% filter(Change_ssp370==0), aes(x=x, y=y), fill="lightblue")+
        theme_bw()+
        theme(axis.title = element_blank(), legend.title = element_blank(),
              legend.position = c(0.1,0.4)))
dev.off()

# ssp585 plot change in distribution
png(file=paste0(here::here(), "/figures/SpeciesRichness_", "2041-2070_change_ssp370_", Taxon_name, ".png"),width=1000, height=1000)
print(ggplot()+
        geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
        xlim(-23, 60) +
        ylim(31, 75) +
        
        geom_tile(data=average_stack %>% filter(Richness==0), aes(x=x, y=y), fill="grey60")+
        geom_tile(data=average_stack %>% filter(!is.na(Change)) %>% filter(FutureRichness!=0 & Richness!=0), 
                  aes(x=x, y=y, fill=Change_f_ssp585))+
        ggtitle(paste0("Change in species richness (number of species) SSP585"))+
        scale_fill_viridis_d(breaks=c("[10,15]", "[5,10]", "[0,5]", "[-5,0]", "[-10,-5]", "[-15,-10]"), option="B")+
        geom_tile(data=average_stack %>% filter(Change_ssp585==0), aes(x=x, y=y), fill="lightblue")+
        theme_bw()+
        theme(axis.title = element_blank(), legend.title = element_blank(),
              legend.position = c(0.1,0.4)))
dev.off()


#- - - - - - - - - - - - - - - - - - - - -
## Variable importance ####
#- - - - - - - - - - - - - - - - - - - - -

## Calculate variable importance (VI) ####

# create result data frame
var_imp <- data.frame("Predictor"= c("Aridity", "MAP", "MAP_Seas", "MAT", 
                                              "MAT_Seas", "Snow", "Agriculture", "Dist_Urban",
                                              "Forest_Coni", "Forest_Deci", "NDVI", 
                                              "Pastures", "Pop_Dens", "Shrubland", "Aspect",
                                              "Dist_Coast", "Dist_River", "Elev", 
                                              "Slope", "CEC", "Clay.Silt", "Cu", "Hg",
                                              "Moisture", "N", "P", "pH", "SOC", "SoilT"),
                      "biomod"=NA, "Species"=NA)

for(spID in unique(speciesNames[speciesNames$NumCells_5km >= 10,]$SpeciesID)){ try({
  
  print("=====================================")
  print(spID)
  
  load(file=paste0(here::here(), "/results/", Taxon_name, "/_SDMs/SDM_biomod_", spID, ".RData")) #biomod_list
  temp_vi <- biomod_list$varImp
  
  # round variable importance to 3 decimals
  temp_vi[,"biomod"] <- round(temp_vi[,"biomod"], 3)
      
  # add species name
  temp_vi$Species <- spID
  
  var_imp <- rbind(var_imp, temp_vi)

  rm(temp_vi, biomod_list)
    
  #var_imp
          
})}

var_imp <- var_imp[!is.na(var_imp$Species),] %>% unique()
rownames(var_imp) <- 1:nrow(var_imp)

var_imp
str(var_imp)

## Save ####
write_csv(var_imp, file=paste0(here::here(), "/results/Variable_importance_biomod_", Taxon_name, ".csv"))
 

var_imp <- read.csv(file=paste0(here::here(), "/results/Variable_importance_biomod_", Taxon_name, ".csv"))
var_imp

# load predictor table to get classification of variables
# load the predictor table containing the individual file names
pred_tab <- readr::read_csv(file=paste0(here::here(), "/doc/Env_Predictors_table.csv"))

# transform to long format and add variable categories
var_imp <- var_imp %>%
  left_join(pred_tab %>% dplyr::select(Predictor, Category), by="Predictor")

# add category for clay.silt
var_imp[var_imp$Predictor=="Clay.Silt","Category"] <- "Soil"

# plot VIF
plotVarImp <- ggplot(data=var_imp, aes(x=biomod, y=reorder(Predictor, biomod), fill=Category))+
  geom_boxplot(cex=0.2, outlier.size=1.5)+
  geom_jitter(height=0.2, alpha=0.3)+
  xlab("Variable importance")+
  ylab("Predictor")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 10))
plotVarImp

png(paste0(here::here(), "/figures/VariableImportance_biomod_", Taxon_name, ".png")); plotVarImp; dev.off()

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

#- - - - - - - - - - - - - - - - - - - - -
## Calculate variable importance for richness plot ####
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
