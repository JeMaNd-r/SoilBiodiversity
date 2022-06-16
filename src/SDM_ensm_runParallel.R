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

# covariates in order of importance (top 10 important)
covarsNames <- c("MAT", "Dist_Coast", "MAP_Seas", "Pop_Dens", "Hg", 
                 "Agriculture", "Elev", "pH", "Pastures", "MAP")

# define future scenarios
scenarioNames <- sort(paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
                               "mri-esm2-0", "ukesm1-0-ll"), "_",
                             rep(c("ssp126", "ssp370", "ssp585"),5)))

#- - - - - - - - - - - - - - - - - - - - -
# note: we will load the datasets before each individual model

# load environmental variables (for projections)
Env_norm <- raster::stack(paste0("./results/EnvPredictor_2km_normalized.grd"))
#Env_norm <- stack(Env_norm)

# as dataframe
load(paste0("./results/EnvPredictor_2km_df_normalized.RData")) #Env_norm_df

# define formula for GLM (and biomod)
form <- paste0("occ ~ ", paste0(paste0("s(", covarsNames, ")"), collapse=" + "))

# Calculate the number of cores
no.cores <-  parallel::detectCores()/2 

#- - - - - - - - - - - - - - - - - - - - -
## Prepare data ####
mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))

registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 10,]$SpeciesID), 
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
 
          
          # model fitting
          tmp <- proc.time()[3]
          setwd(paste0("./results/", Taxon_name))
          
          set.seed(32639)
          myBiomodModelOut <- biomod2::BIOMOD_Modeling(myBiomodData,
                                                       models = mymodels,
                                                       models.options = myBiomodOption,
                                                       NbRunEval = 10,   # 3-fold crossvalidation evaluation run
                                                       DataSplit = 80, # use subset of the data for training
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
          temp_model_time <- proc.time()[3] - tmp
          
          setwd(here::here())
})}
stopImplicitCluster()

#- - - - - - - - - - - - - - - - - - - - -
# if more than 10 but less than 100 occurrences   
registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 10 & speciesNames$NumCells_2km < 100,]$SpeciesID), 
        .export = c("Env_norm", "Env_norm_df", "form"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({
          
          load(paste0(here::here(), "/results/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData")) #myBiomodData
          
          # subset covarsNames
          myBiomodData@data.env.var <- myBiomodData@data.env.var[,colnames(myBiomodData@data.env.var) %in% covarsNames]
          
	        models.options <- BIOMOD_ModelingOptions()
          models.options@GBM$n.trees <- 1000
          models.options@GBM$interaction.depth <- 4
          models.options@GBM$shrinkage <- 0.005
          models.options@GAM$select <- TRUE
          models.options@CTA$control$cp <- 0
          models.options@ANN$size <- 8
          models.options@ANN$decay <- 0.001
          models.options@MARS$interaction.level <- 0
          models.options@MARS$nprune <- 2
          models.options@MAXENT.Phillips$product <- FALSE
          models.options@MAXENT.Phillips$threshold <- FALSE
          models.options@MAXENT.Phillips$betamultiplier <- 0.5
	        models.options@MAXENT.Phillips$path_to_maxent.jar <- paste0(here::here())
          models.options@GLM$test <- "none"

          trControl <- caret::trainControl(method = "repeatedcv", summaryFunction = caret::twoClassSummary, 
                                           classProbs = T, returnData = F, repeats = 10)
          
          models.options <- BIOMOD_tuning(data = myBiomodData, 
                                          models = mymodels[mymodels != "RF"], 
                                          models.options = models.options)$models.options
          
          # model fitting
          tmp <- proc.time()[3]
          setwd(paste0("./results"))
          
          set.seed(32639)
          
          myBiomodModelOut <- BIOMOD_Modeling(data = myBiomodData, models = mymodels, models.options = models.options, 
                models.eval.meth = "TSS", DataSplit = 80, NbRunEval = 10, rescal.all.models = FALSE, 
                do.full.models = TRUE, VarImport = 5, modeling.id = paste(spID,"_Modeling", sep = ""))

          myBiomodModelOut <- ecospat.ESM.Modeling.fixed(data = myBiomodData,
                                                   models = mymodels,
                                                   models.options = models.options,
                                                   Prevalence = NULL,
                                                   tune = TRUE, # TRUE: estimate optimal parameters for the models
                                                   NbRunEval = 10,   # 3-fold crossvalidation evaluation run
                                                   DataSplit = 80, # use subset of the data for training
                                                   weighting.score = "TSS",
                                                   #ESM_Projection = FALSE, #no projections now (will be done later)
                                                   cleanup = 2, #when to delete temporary unused files, in hours
                                                   modeling.id = paste(myRespName,"_Modeling", sep = ""))
          
          
          # ensemble modeling
          myBiomodEM <- ecospat.ESM.EnsembleModeling(ESM.modeling.output = myBiomodModelOut,
                                                     weighting.score = "TSS")
          
          temp_model_time <- proc.time()[3] - tmp
          
})}
stopImplicitCluster()

#- - - - - - - - - - - - - - - - -
## Predict in current climate ####

registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 10,]$SpeciesID), 
        .export = c("Env_norm", "Env_norm_df", "form"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({   
          
          # list files in species-specific BIOMOD folder
          temp_files <- list.files(paste0("./results/", Taxon_name, "/", stringr::str_replace(spID, "_", ".")), full.names = TRUE)
          
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
          
          temp_validation <- myEnProjDF[,1]
          temp_validation <- as.data.frame(temp_validation)
          temp_validation$x <- training$bg.biomod@coord$x
          temp_validation$y <- training$bg.biomod@coord$y
          temp_validation <- temp_validation %>% rename("layer"=temp_validation)
          temp_validation$layer <- temp_validation$layer / 1000
          
          # Get model evaluation values for later
          myBiomodModelEval <- as.data.frame(biomod2::get_evaluations(myBiomodEM)[2])
          
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
          
          biomod_list <- list(time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp, evaluation=myBiomodModelEval)
          save(biomod_list, file=paste0("./SDMs/SDM_biomod_", spID, ".RData"))
          
          rm(biomod_list, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp, myBiomodEnProj, myBiomodProj, myBiomodModelEval, myEnProjDF)
          
        })}
stopImplicitCluster()


#- - - - - - - - - - - - - - - - -
## Predict in future climate ####

setwd(here::here())

registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 10,]$SpeciesID), 
        .export = c("Env_norm", "Env_norm_df", "form"),
        .packages = c("tidyverse","biomod2")) %dopar% { try({ 
          
          # list files in species-specific BIOMOD folder
          temp_files <- list.files(paste0("./results/", Taxon_name, "/", stringr::str_replace(spID, "_", ".")), full.names = TRUE)
          
          # load model output
          myBiomodModelOut <- temp_files[stringr::str_detect(temp_files,"Modeling.models.out")]
          print(myBiomodModelOut)
          myBiomodModelOut <-get(load(myBiomodModelOut))
          
          # load ensemble model output
          myBiomodEM <- temp_files[stringr::str_detect(temp_files,"Modelingensemble.models.out")]
          myBiomodEM <- get(load(myBiomodEM))
          
          for(no_future in scenarioNames){
            
            # load env. data with future climate (MAP, MAT, MAP_Seas)
            load(paste0("./results/EnvPredictor_2041-2070_", no_future, "_2km_df_normalized.RData")) #temp_Env_df
            
            setwd(paste0(here::here(), "/results/", Taxon_name))
            
            ## NOTE: because biomod output can hardly be stored in list file, we will do calculations based on model output now
            # project single models (also needed for ensemble model)
            myBiomodProj <- biomod2::BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                                       new.env = temp_Env_df[,colnames(temp_Env_df) %in% covarsNames],        #column/variable names have to perfectly match with training
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
            
            # save predictions as raster file
            temp_prediction <- myBiomodEnProj@proj@val[,2]
            temp_prediction <- as.numeric(temp_prediction)
            # add names of grid cell (only for those that have no NA in any layer)
            names(temp_prediction) <- rownames(temp_Env_df)
            temp_prediction <- as.data.frame(temp_prediction)
            temp_prediction$x <- temp_Env_df$x
            temp_prediction$y <- temp_Env_df$y
            temp_prediction <- temp_prediction %>% full_join(temp_Env_df %>% dplyr::select(x,y)) %>%
              rename("layer" = temp_prediction)
            temp_prediction$layer <- temp_prediction$layer / 1000
            
            setwd(here::here())
            save(temp_prediction, paste0("./results/_Maps/SDM_2041-2070_", no_future, "_biomod_", spID,  ".RData")) 
            rm(temp_prediction, temp_Env_df, myBiomodEnProj, myBiomodProj)
          }
          
        
        
})}
stopImplicitCluster()


#- - - - - - - - - - - - - - - - - - - - -
# if less than 100 occurrences (but more than 10)
if(speciesNames[speciesNames$SpeciesID==spID, "NumCells_2km"]<100){ 
  
  
  
  
  #- - - - - - - - - - - - - - - - - - - - -
  ## Check if everything went well ####
  #- - - - - - - - - - - - - - - - - - - - -
  
  for(spID in unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID)){ try({
    
    check_files <- list.files(paste0("./results/", Taxon_name, "/temp_files"))
    check_files <- check_files[stringr::str_detect(check_files, spID)]
    
    print(check_files)
    if(length(check_files)==17){ print("Everything looks well. Please continue :)")
    }else{ print("Some algorithms weren't saved; please check which ones are missing.")}
    
})}
  
stopImplicitCluster()
  
  