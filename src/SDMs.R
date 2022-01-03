#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# 1. Selecting data corresponding to a species

# 2. Putting this data in the biomod2 format (i.e. BIOMOD_FormatingData)

# 3. Building individual models (i.e. BIOMOD_Modeling)

# 4. Building ensemble-models (i.e. BIOMOD_EnsembleModeling)

# 5. Making model projections (i.e. BIOMOD_Projection and BIOMOD_EnsembleForecasting)

#- - - - - - - - - - - - - - - - - - - - - -#
## Parallelize using snowfall package

## Biodiversity data (occurrences, response variable)
mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrences_wide_", Taxon_name, ".csv"))

#- - - - - - - - - - - - - - - - - - - - - -#
## Define species names
# check first if there are enough occurrence data for each species
sp.occ <- colSums(mySpeciesOcc[,speciesNames$SpeciesID[speciesNames$SpeciesID %in% colnames(mySpeciesOcc)]], na.rm=T) >=5

# define species 
sp.names <- speciesNames[sp.occ, "SpeciesID"]

#- - - - - - - - - - - - - - - - - - - - - -#
## Environmental (explanatory) variables as raster file
myExpl <- stack(paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))

# crop to Europe
for(i in 1:nlayers(myExpl)) {myExpl[[i]] <- raster::mask(myExpl[[i]], as(extent(extent_Europe), 'SpatialPolygons'))}

# set working directory
setwd(paste0(here::here(), "/results/SDMs"))

#- - - - - - - - - - - - - - - - - - - - - -#
# Define function to create SDMs
MyBiomodSF <- function(sp.n){
  
  cat('\n',sp.n,'modelling...') #tell which species is being modelled
  
  ## Define data for this run
  # i.e keep only the column of our species (response)
  myResp <- as.numeric(mySpeciesOcc[,sp.n])
  
  # get NAs id
  na.id <- which(is.na(myResp))
  
  # remove NAs to enforce PA sampling to be done on explanatory rasters
  myResp <- myResp[-na.id]
  myRespCoord = mySpeciesOcc[-na.id,c('x','y')]
  myRespName = sp.n
  
  ## Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName,
                                       PA.nb.rep = 2,
                                       PA.nb.absences = 10*sum(myResp==1,
                                                               na.rm=TRUE),
                                       PA.strategy = 'random')
  
  # options definition
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  ## Modelling
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c("GLM", "GBM", "GAM", "CTA", "ANN", "SRE", "FDA", "MARS", "RF",
               "MAXENT.Phillips", "MAXENT.Phillips.2"),
    models.options = myBiomodOption,
    NbRunEval=1,
    DataSplit=80,
    Yweights=NULL,
    VarImport=3,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = TRUE)
  
  # add an error catch if all models are removed due to threshold filtering
  t <- try({
  
    ## Building ensemble-models
    myBiomodEM <- BIOMOD_EnsembleModeling(
      modeling.output = myBiomodModelOut,
      eem.by = 'all'
      chosen.models = 'all',
      eval.metric = 'TSS',
      eval.metric.quality.threshold = c(0.85),
      prob.mean = T,
      prob.cv = T,
      prob.ci = T,
      prob.ci.alpha = 0.05,
      prob.median = T,
      committee.averaging = T,
      prob.mean.weight = T,
      prob.mean.weight.decay = 'proportional',
      VarImport=3)
  })
  
  # if ensemble modelling gets an error, stop here
  if("try-error" %in% class(t)){
    
    # clean environment
    rm(myBiomodData, myBiomodModelOut)
    
    # print error message
    print("There are no models kept for EnsembleModeling due to threshold filtering. 
          Please check your environmental predictors or adapt the threshold if necessary.")
  
  }else{ #if there is no error, continue
    
  # do projections on current varaiable
  myBiomomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = stack(myExpl),
    proj.name = 'current',
    selected.models = 'all',
    binary.meth= 'ROC',
    compress = 'xz',
    clamping.mask = F)
  
  # try if all models work
  t2 <- try({
      
    # Do ensemble-models projections on current varaiable
    myBiomodEF <- BIOMOD_EnsembleForecasting(
      projection.output = myBiomomodProj,
      EM.output = myBiomodEM,
      binary.meth = 'TSS',
      total.consensus = TRUE)
  })
  
  while("try-error" %in% class(t2)){
    t2 <- try({
      sample.models <- sample(c(get_built_models(myBiomodEM)), length(c(get_built_models(myBiomodEM)))-1)
      
      # Do ensemble-models projections on current varaiable
      myBiomodEF <- BIOMOD_EnsembleForecasting(
        projection.output = myBiomomodProj,
        selected.models = sample.models,
        EM.output = myBiomodEM,
        binary.meth = 'TSS',
        total.consensus = TRUE)
    })
    
  writeLines(sample.models, paste0("Models_used_in_EnsembleModel_", sp.n, ".txt"))  
  writeLines(setdiff(c(get_built_models(myBiomodEM)), sample.models), paste0("Models_NOT_used_in_EnsembleModel_", sp.n, ".txt")) 
    
  }
  
  # remove all temporal model outputs
  rm(myBiomodData, myBiomodModelOut, myBiomodEM, myBiomomodProj, myBiomodEF, t, t2, sample.models)
  
  } #end of error alternative
}

#- - - - - - - - - - - - - - - - - - - - - -#
## Initalize the snowfall package
sfInit(parallel=TRUE, cpus=as.numeric(setCPU))

# export packages
sfLibrary('biomod2', character.only=TRUE)

# export variables
sfExport('mySpeciesOcc')
sfExport('myExpl')
sfExport('sp.names')
# you may also use sfExportAll() to exprt all your workspace variables

## Run the actual models
mySFModelsOut <- sfLapply( sp.names, MyBiomodSF)

# stop snowfall
sfStop( nostop=FALSE )

# Clean your environment
rm(list=ls())

#- - - - - - - - - - - - - - - - - - - - - -#
## other approach based on Valavi et al. 2021
#- - - - - - - - - - - - - - - - - - - - - -#

# # specifying the species id
# spID <- "SpeciesID"
# 
# # re-loading the species data 
# pr <- disPo("NSW")
# bg <- disBg("NSW")
# pr <- pr[pr$spid == spID, ] # subset the target species 
# training <- rbind(pr, bg)
# training$vegsys <- as.factor(training$vegsys)
# 
# myRespName <- "occ"
# myResp <- as.numeric(training[, myRespName])
# myResp[which(myResp == 0)] <- NA
# myExpl <- data.frame(training[, covars])
# myRespXY <- training[, c("x", "y")]
# 
# # create biomod data format
# biomod.data <- BIOMOD_FormatingData(resp.var = myResp,
#                                      expl.var = myExpl,
#                                      resp.name = myRespName,
#                                      resp.xy = myRespXY,
#                                      PA.nb.absences = 10000,
#                                      PA.strategy = 'random',
#                                      na.rm = TRUE)
# 
# # using the default options
# # you can change the mentioned parameters by changes this
# myBiomodOption <- BIOMOD_ModelingOptions()
# 
# # models to predict with
# mymodels <- c("GLM","GBM","GAM","CTA","ANN","FDA","MARS","RF","MAXENT.Phillips")
# 
# # model fitting
# tmp <- Sys.time()
# set.seed(32639)
# myBiomodModelOut <- BIOMOD_Modeling(biomod.data,
#                                     models = mymodels,
#                                     models.options = myBiomodOption,
#                                     NbRunEval = 1, 
#                                     DataSplit = 100, # use all the data for training
#                                     models.eval.meth = c("ROC"), 
#                                     SaveObj = TRUE,
#                                     rescal.all.models = FALSE,
#                                     do.full.models = TRUE,
#                                     modeling.id = paste(myRespName,"NCEAS_Modeling", sep = ""))
# # ensemble modeling using mean probability
# myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
#                                       chosen.models = 'all',
#                                       em.by = 'all',
#                                       eval.metric = c("ROC"),
#                                       eval.metric.quality.threshold = NULL, # since some species's auc are naturally low
#                                       prob.mean = TRUE,
#                                       prob.cv = FALSE,
#                                       prob.ci = FALSE,
#                                       prob.median = FALSE,
#                                       committee.averaging = FALSE,
#                                       prob.mean.weight = FALSE)
# 
# Sys.time() - tmp
# 
# Time difference of 4.412951 mins
# 
# # project single models
# myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
#                                   new.env = as.data.frame(testing_env[, covars]),
#                                   proj.name = "nceas_modeling",
#                                   selected.models = "all",
#                                   binary.meth = "ROC",
#                                   compress = TRUE,
#                                   clamping.mask = TRUE)
# 
# # project ensemble of all models
# myBiomodEnProj <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
#                                              EM.output = myBiomodEM,
#                                              selected.models = "all")
# 
# # extracting the values for ensemble prediction
# myEnProjDF <- as.data.frame(get_predictions(myBiomodEnProj))
# 
# # see the first few pridictions
# # the prediction scale of biomod is between 0 and 1000
# head(myEnProjDF)

