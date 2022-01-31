#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#             Predict SDMs                  #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# note: we will load the datasets and models before each individual model

# load models
load(file=paste0(here::here(), "/results/", Taxon_name, "/Predicted_SDMs_", spID, ".RData"))

# load environmental space
Env <- stack(paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))
Env <- crop(Env, extent_Europe) # crop to Europe

# function to rescale raster between 0 and 1
fct.rescale <- function(x, x.min, x.max, new.min = 0, new.max = 1) {
  if(is.null(x.min)) {x.min = min(x)}
  if(is.null(x.max)) {x.max = max(x)}
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

#- - - - - - - - - - - - - - - - - - - - - -
## Start predicting ####

par(mfrow=c(3,3))

## GAM
gm <- SDMs[["gm_pred"]][[1]]
gm_prediction <- raster::predict(Env, gm)

# rescale between 0 and 1
gm_prediction <- fct.rescale(gm_prediction, x.min = gm_prediction@data@min, x.max = gm_prediction@data@max)

plot(gm_prediction, main = "GAM")


## GLM
lm1 <- SDMs[["lm1_pred"]][[1]]
lm1_prediction <- raster::predict(Env, lm1)
# rescale between 0 and 1
lm1_prediction <- fct.rescale(lm1_prediction, x.min = lm1_prediction@data@min, x.max = lm1_prediction@data@max)

plot(lm1_prediction, main = "GLM")

# second GLM
lm_subset <- SDMs[["lm_subset_pred"]][[1]]
lm_subset_prediction <- raster::predict(Env, lm_subset)
lm_subset_prediction <- fct.rescale(lm_subset_prediction, x.min = lm_subset_prediction@data@min, x.max = lm_subset_prediction@data@max)

plot(lm_subset_prediction, main = "GLM subset")


## Lasso
lasso_cv <- SDMs[["lasso_pred"]][[1]]

# we have to load again the training data set
# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0("bg.glm", "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# generating the quadratic terms for all continuous variables
# function to creat quadratic terms for lasso and ridge
quad_obj <- make_quadratic(training, cols = covarsNames)

# now we can predict this quadratic object on the training and testing data
# this make two columns for each covariates used in the transformation
testing_quad <- predict(quad_obj, newdata = Env)

# convert the data.frames to sparse matrices
# select all quadratic (and non-quadratic) columns, except the y (occ)
testing_sparse <- sparse.model.matrix( ~. -1, testing_quad[, covarsNames])

lasso_prediction <- predict(lasso_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
#head(lasso_prediction)

lasso_prediction <- as.numeric(lasso_prediction)
names(lasso_prediction) <- rownames(Env) #add site names


## Ridge
ridge_cv <- SDMs[["ridge_pred"]][[1]]
ridge_prediction <- predict(ridge_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
#head(ridge_prediction)

ridge_prediction <- as.numeric(ridge_prediction)
names(ridge_prediction) <- rownames(Env) #add site names


## MARS
mars_fit <- SDMs[["mars_pred"]][[1]][[sample(1:10, 1)]]
mars_pred <- raster::predict(Env, mars_fit)
# transform occurrence back into numeric
mars_pred <- as.character(mars_pred)
mars_pred[mars_pred=="C0"] <- 0
mars_pred[mars_pred=="C1"] <- 1
mars_pred <- as.numeric(mars_pred)
names(mars_pred) <- rownames(Env) #add site names
#head(mars_pred)

# define background dataset (for testing data)
modelName <- SDMs[["mars_pred"]][[3]]

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# load model
model <- SDMs[[temp.model]][[1]]
  

## MAXENT
maxmod <- SDMs[["maxmod_pred"]][[1]]
maxmod_prediction <- dismo::predict(maxmod, Env)
maxmod_prediction <- fct.rescale(maxmod_prediction, x.min = maxmod_prediction@data@min, x.max = maxmod_prediction@data@max)

plot(maxmod_prediction, main = "Maxent")


# predicting with MAXNET
mxnet <- SDMs[["maxnet_pred"]][[1]]
maxnet_prediction <- raster::predict(Env, mxnet)
maxnet_prediction <- fct.rescale(maxnet_prediction, x.min = maxnet_prediction@data@min, x.max = maxnet_prediction@data@max)

plot(maxnet_prediction, main = "MaxNet")


## BRT
brt <- SDMs[["brt_pred"]][[1]][,sample(1, 1:ncol(SDMs[["brt_pred"]][[1]]))]
# predicting with the best trees
brt_pred <- raster::predict(Env, brt)
#head(brt_pred)

brt_pred <- as.numeric(brt_pred)
names(brt_pred) <- rownames(Env) #add site names


## XGBoost
xbg_fit <- SDMs[["xgb_pred"]][[1]]
xgb_prediction <- raster::predict(Env, xgb_fit)

# transform occurrence back into numeric
xgb_prediction <- as.character(xgb_prediction)
xgb_prediction[xgb_prediction=="C0"] <- 0
xgb_prediction[xgb_prediction=="C1"] <- 1
xgb_prediction <- as.numeric(xgb_prediction)
names(xgb_prediction) <- rownames(Env) #add site names


xgb_prediction <- fct.rescale(xgb_prediction, x.min = maxmod_prediction@data@min, x.max = maxmod_prediction@data@max)

plot(xgb_prediction, main = "XGBoost")



## RF
rf <- SDMs[["rf_pred"]][[1]][,sample(1, 1:ncol(SDMs[["brt_pred"]][[1]]))]
rf_prediction <- raster::predict(Env, rf, type = "prob") # prob = continuous prediction
#head(rf_prediction)

rf_prediction <- as.numeric(rf_prediction)
names(rf_prediction) <- rownames(Env) #add site names

# RF downsampled
rf_downsample <- SDMs[["rf_downsample_pred"]][[1]]
rf_downsample_prediction <- predict(rf_downsample, Env, type = "prob")[, "1"] # prob = continuous prediction
#head(rf_downsample_prediction)

rf_downsample_prediction <- as.numeric(rf_downsample_prediction)
names(rf_downsample_prediction) <- rownames(Env) #add site names


## SVM
svm_e <- SDMs[["svm_pred"]][[1]]
svm_prediction <- raster::predict(Env, svm_e, probability = TRUE)
svm_prob <- attr(svm_prediction, "probabilities")[,"1"]

# see the first few predictions
#head(svm_prob)

svm_prob <- as.numeric(svm_prob)
names(svm_prob) <- rownames(Env) #add site names


## BIOMOD
myBiomodEM <- SDMs[["biomod_pred"]][[1]][[1]]
myBiomodModelOut <- SDMs[["biomod_pred"]][[1]][[2]]

# project single models
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                  new.env = as.data.frame(Env[, covarsNames]),
                                  proj.name = "modeling",
                                  selected.models = "all",
                                  binary.meth = "ROC",
                                  compress = TRUE,
                                  clamping.mask = TRUE)

# project ensemble of all models
myBiomodEnProj <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                             EM.output = myBiomodEM,
                                             selected.models = "all")


## simple Ensemble model
# create Raster stack of the models to be merged
stack_prediction <- raster::stack(gm_prediction, maxmod_prediction) #lasso_prediction, brt2_prediction, rf2_prediction

# average the predictions
ensm_prediction <- calc(stack_prediction, fun = mean, na.rm = T)


#- - - - - - - - - - - - - - - - - - - - - -
## save ####

# save predictions
#ave(model_pred, file=paste0(here::here(), "/results/SDM_Predictions_", Taxon_name, "_", spID, ".RData"))



## -------------------------------------------------------------------------------------- ##
## Obtain spatiotemporal projections ----
## -------------------------------------------------------------------------------------- ##
# Goncalas et al. ...
# https://github.com/joaofgoncalves/GoncalvesAna_et_al_2021/tree/master/RCODE/PostModelAnalyses


# Models to consider in the ensemble and projection
modelsToUse <- get_kept_models(myBiomodEM, 1)



for(projName in projNames){
  
  # Obtain spatiotemporal projections
  myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                    new.env = get(projName),
                                    proj.name = projName, ## Name of the projection from above variable proj.name
                                    selected.models = modelsToUse,
                                    filtered.meth = NULL,
                                    binary.meth = NULL,
                                    compress = 'gzip',
                                    clamping.mask = TRUE,
                                    output.format = '.grd',
                                    do.stack = TRUE)
  
  
  # Perform the ensembling of projections
  myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                           binary.meth = c('TSS','ROC','KAPPA'),
                                           EM.output = myBiomodEM,
                                           output.format = '.grd')
  
  # Convert all output raster files to GeoTIFF
  inFolder <- paste(getwd(),"/",sp,"/proj_",projName,sep="")
  outFolder <- paste(inFolder,"/","GeoTIFF", sep="")
  dir.create(outFolder)
  
  convertToGeoTIFF(inFolder, outFolder)
  
} 

save.image(file=paste(sp,"ModObjects.RData",sep="_"))
