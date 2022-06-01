#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

setwd("~/share/groups/eie/==PERSONAL/RZ_SoilBON/SoilBiodiversity_RStudio")

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8g")) # expand Java memory
gc()
library(tidyverse)
library(here)

# read functions for ensemble of small moels (ESM)
devtools::source_url("https://raw.githubusercontent.com/cran/ecospat/master/R/ecospat.ESM.R")
#library(ecospat)

library(raster)

library(biomod2) # also to create pseudo-absences

library(usdm) # for variable inflation factor calculation

library(mgcv) # for GAM
library(gam)  # for GLM (!)
#library(remotes) #to download package from github


# download maxent.jar 3.3.3k, and place the file in the
# desired folder
# utils::download.file(url = "https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar", 
#     destfile = paste0(system.file("java", package = "dismo"), 
#         "/maxent.jar"), mode = "wb")  ## wb for binary file, otherwise maxent.jar can not execute

library(scales) # for Ensemble model

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

corMatPearson <- as.matrix(read.csv(file=paste0("./results/corMatPearson_predictors.csv")))
dimnames(corMatPearson)[[1]] <- dimnames(corMatPearson)[[2]]
# based on Valavi et al. 2021: Pearson 0.8
env_exclude <- caret::findCorrelation(corMatPearson, cutoff = 0.8, names=TRUE)
covarsNames <- dimnames(corMatPearson)[[1]][!(dimnames(corMatPearson)[[1]] %in% env_exclude)]
covarsNames <- covarsNames[covarsNames != "x" & covarsNames != "y"]
# exclude based on VIF
env_vif <- read.csv(file=paste0("./results/VIF_predictors.csv"))
env_exclude <- env_vif %>% filter(is.na(VIF)) %>% dplyr::select(Variables) %>% as.character()
covarsNames <- covarsNames[!(covarsNames %in% env_exclude)]
# excluded:
print("=== We excluded the following variables based on VIF and Pearson correlation: ===")
setdiff(env_vif$Variables, covarsNames)

# final predictor variables
print("=== And we kept the following, final predictor variables: ===")
covarsNames

# define future scenarios
scenarioNames <- paste0(c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", 
                          "mri-esm2-0", "ukesm1-0-ll"), "_",
                          rep(c("ssp126", "ssp370", "ssp585"),3))

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


registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 10,]$SpeciesID), 
         .export = c("Env_norm", "Env_norm_df", "form"),
         .packages = c("tidyverse","biomod2")) %dopar% { try({
             
  #- - - - - - - - - - - - - - - - - - - - -
  ## biomod ####
  #- - - - - - - - - - - - - - - - - - - - -
  
  modelName <- "bg.biomod"
  
  # identify and load all relevant data files
  temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                           pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
  
  lapply(temp.files, load, .GlobalEnv)
  
  # create biomod data format
  myBiomodData <- training$bg.biomod
  
  # subset covarsNames
  myBiomodData@data.env.var <- myBiomodData@data.env.var[,colnames(myBiomodData@data.env.var) %in% covarsNames]
  
  myRespName <- "occ"
  
  # # create own function for GAM and MAXENT.Phillips
  # myBiomodOption <- BIOMOD_ModelingOptions(
  #   GAM = list( algo = 'GAM_mgcv',
  #               type = 's_smoother',
  #               k = -1,
  #               interaction.level = 0,
  #               myFormula = as.formula(form),
  #               family = binomial(link = 'logit'),
  #               method = 'GCV.Cp',
  #               optimizer = c('outer','newton'),
  #               select = FALSE,
  #               knots = NULL,
  #               paraPen = NULL,
  #               control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07, maxit = 200, trace = FALSE,
  #                              mgcv.tol = 1e-07, mgcv.half = 15, rank.tol = 1.49011611938477e-08,
  #                              nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0),
  #                              optim = list(factr=1e+07),
  #                              newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0), outerPIsteps = 0,
  #                              idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE, scale.est = "fletcher",
  #                              edge.correct = FALSE)),
  #   MAXENT.Phillips = list( path_to_maxent.jar = paste0("./results/", Taxon_name, "/maxent_files"), # change it to maxent directory
  #                           memory_allocated = 512,
  #                           background_data_dir = 'default',
  #                           maximumbackground = 50000,
  #                           maximumiterations = 200,
  #                           visible = FALSE,
  #                           linear = TRUE,
  #                           quadratic = TRUE,
  #                           product = TRUE,
  #                           threshold = TRUE,
  #                           hinge = TRUE,
  #                           lq2lqptthreshold = 80,
  #                           l2lqthreshold = 10,
  #                           hingethreshold = 15,
  #                           beta_threshold = -1,
  #                           beta_categorical = -1,
  #                           beta_lqp = -1,
  #                           beta_hinge = -1,
  #                           betamultiplier = 1,
  #                           defaultprevalence = 0.5)
  # )
  
  # alternative: default algorithms
  # Note: try this out to see if GAM works
  myBiomodOption <- BIOMOD_ModelingOptions(
    GAM = list (k = -1), #avoid error messages
    #MAXENT.Phillips = list( path_to_maxent.jar =paste0("./results" )), # change it to maxent directory
  )
  
  # models to predict with
  mymodels <- c("GLM","GBM","GAM","CTA","ANN","FDA","MARS","RF","MAXENT.Phillips")
  
  # if less than 100 occurrences (but more than 10)
  if(speciesNames[speciesNames$SpeciesID==spID, "NumCells_2km"]<100){ 
    
    # model fitting
    tmp <- proc.time()[3]
    setwd(paste0("./results/", Taxon_name))
    
    set.seed(32639)
  
    myBiomodModelOut <- ecospat.ESM.Modeling(data = myBiomodData,
                                             models = mymodels,
                                             models.options = myBiomodOption,
                                             Prevalence = NULL,
                                             tune = TRUE, # TRUE: estimate optimal parameters for the models
                                             NbRunEval = 10,   # 3-fold crossvalidation evaluation run
                                             DataSplit = 80, # use subset of the data for training
                                             weighting.score = "TSS",
                                             ESM_Projection = FALSE, #no projections now (will be done later)
                                             cleanup = 2, #when to delete temporary unused files, in hours
                                             modeling.id = paste(myRespName,"_Modeling", sep = ""))
    
    # ensemble modeling
    myBiomodEM <- ecospat.ESM.EnsembleModeling(ESM.modeling.output = myBiomodModelOut,
                                               weighting.score = "TSS")
    
    temp_model_time <- proc.time()[3] - tmp
    
    
    tmp <- proc.time()[3]
    
    
  #- - - - - - - - - - - - - - - - - - - - -
  # if more than 100 occurrences    
  }else{
    
    # model fitting
    tmp <- proc.time()[3]
    setwd(paste0("./results/", Taxon_name))
    
    set.seed(32639)
    myBiomodModelOut <- biomod2::BIOMOD_Modeling(myBiomodData,
                                                 models = mymodels,
                                                 models.options = myBiomodOption,
                                                 NbRunEval = 10,   # 3-fold crossvalidation evaluation run
                                                 DataSplit = 80, # use subset of the data for training
                                                 models.eval.meth = c("TSS"),
                                                 SaveObj = FALSE, #save output on hard drive?
                                                 rescal.all.models = FALSE, #scale all predictions with binomial GLM?
                                                 do.full.models = FALSE, # do evaluation & calibration with whole dataset
                                                 modeling.id = paste(myRespName,"_Modeling", sep = ""))
    
    # ensemble modeling using mean probability
    myBiomodEM <- biomod2::BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                                   chosen.models = "all",  # all algorithms
                                                   em.by = "all",    #evaluated over evaluation data if given (it is not, see Prepare_input.R)
                                                   # note: evaluation not that important as we will calculate measures on independent data
                                                   eval.metric = c("TSS"), # 'all' would takes same as above in BIOMOD_Modelling
                                                   eval.metric.quality.threshold = NULL, # since some species's auc are naturally low
                                                   prob.mean = FALSE, #estimate mean probabilities across predictions
                                                   prob.cv = TRUE,   #estimate coefficient of variation across predictions
                                                   prob.ci = FALSE,  #estimate confidence interval around the prob.mean
                                                   prob.median = FALSE, #estimate the median of probabilities
                                                   committee.averaging = FALSE, #estimate committee averaging across predictions
                                                   prob.mean.weight = TRUE, #estimate weighted sum of predictions
                                                   prob.mean.weight.decay = "proportional", #the better a model (performance), the higher weight
                                                   VarImport = 1)    #number of permutations to estimate variable importance
    temp_model_time <- proc.time()[3] - tmp
     
    
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
     
    
    # extracting the values for ensemble prediction
    myEnProjDF <- as.data.frame(get_predictions(myBiomodEM)[,2]) #for weighted probability mean
    
    # see the first few predictions
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
    
    biomod_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp, evaluation=myBiomodModelEval)
    save(biomod_list, file=paste0("./results/", Taxon_name, "/temp_files/SDM_biomod_", spID, ".RData"))
    
    rm(biomod_list, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp, myBiomodEnProj, myBiomodProj, myBiomodModelEval, myEnProjDF)
    
    for(no_future in scenarioNames){
      
      # load env. data with future climate (MAP, MAT, MAP_Seas)
      load(paste0("./results/EnvPredictor_", no_future, "_2km_df_normalized.RData")) #temp_Env_df
      
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
      
      save(temp_prediction, paste0("./results/_Maps/SDM_", no_future, "_biomod_", spID,  ".RData"))
    }
    
   rm(temp_prediction, temp_Env_df, myBiomodEnProj, myBiomodProj)
    
  }
   
  setwd(here::here())
})}
stopImplicitCluster()


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

