#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#             Predict SDMs                  #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# Note: we will load the datasets and models before predicting each individual model

# load models
load(file=paste0(here::here(), "/sdm/SDM_Models_", spID, ".RData"))

# load environmental space
Env <- stack(paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))
Env <- crop(Env, extent_Europe) # crop to Europe

# define function to rescale raster (for predicted occurrence) between 0 and 1
fct.rescale <- function(x, x.min, x.max, new.min = 0, new.max = 1) {
  if(is.null(x.min)) {x.min = min(x)}
  if(is.null(x.max)) {x.max = max(x)}
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

#- - - - - - - - - - - - - - - - - - - - - -
## Start predicting ####
#- - - - - - - - - - - - - - - - - - - - - -

par(mfrow=c(4,4))

#- - - - - - - - - - - - - - - - - - - - - -
## GAM ####
gm <- SDMs[["gm_pred"]][[1]]
gm_pred <- raster::predict(Env, gm)

# rescale between 0 and 1
gm_pred <- fct.rescale(gm_pred, x.min = gm_pred@data@min, x.max = gm_pred@data@max)

plot(gm_pred, main = "GAM")

#- - - - - - - - - - - - - - - - - - - - - -
## GLM ####
lm1 <- SDMs[["lm1_pred"]][[1]]
lm1_pred <- raster::predict(Env, lm1)
# rescale between 0 and 1
lm1_pred <- fct.rescale(lm1_pred, x.min = lm1_pred@data@min, x.max = lm1_pred@data@max)

plot(lm1_pred, main = "GLM")

# second GLM
lm_subset <- SDMs[["lm_subset_pred"]][[1]]
lm_subset_pred <- raster::predict(Env, lm_subset)
lm_subset_pred <- fct.rescale(lm_subset_pred, x.min = lm_subset_pred@data@min, x.max = lm_subset_pred@data@max)

plot(lm_subset_pred, main = "GLM subset")

#- - - - - - - - - - - - - - - - - - - - - -
## Lasso ####
lasso_cv <- SDMs[["lasso_pred"]][[1]]

# we have to load again the training data set
# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0("bg.glm", "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# generating the quadratic terms for all continuous variables
# function to creat quadratic terms for lasso and ridge
quad_obj <- make_quadratic(training, cols = covarsNames)

# for the next step, we need the Env as data frame
Env.df <- as.data.frame(Env)

# now we can predict this quadratic object on the training and testing data
# this make two columns for each covariates used in the transformation
testing_quad <- predict(quad_obj, newdata = Env.df)

# Define vector with appropriate covarsNames for sparse matrix.
# More specifically: create names like this: bio1_1, bio1_2, bio2_1, bio2_2...
sparse_covarsNames <- paste0(rep(covarsNames, each=2), "_", 1:2)

# convert the data.frames to sparse matrices
# select all quadratic (and non-quadratic) columns, except the y (occ)
testing_sparse <- sparse.model.matrix( ~. -1, testing_quad[, sparse_covarsNames])

lasso_pred <- predict(lasso_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
#head(lasso_pred)

lasso_pred <- as.numeric(lasso_pred)
names(lasso_pred) <- rownames(Env) #add site names

plot(lasso_pred, main = "Lasso regression")

#- - - - - - - - - - - - - - - - - - - - - -
## Ridge ####
ridge_cv <- SDMs[["ridge_pred"]][[1]]
ridge_pred <- predict(ridge_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
#head(ridge_pred)

ridge_pred <- as.numeric(ridge_pred)
names(ridge_pred) <- rownames(Env) #add site names

plot(ridge_pred, main = "Ridge regression")

#- - - - - - - - - - - - - - - - - - - - - -
## MARS ####
# load pre-calculated predictions
mars_pred <- SDMs[["mars_pred"]][[6]]
#mars_pred

# define background dataset (for testing data)
modelName <- SDMs[["mars_pred"]][[3]]

plot(mars_pred, main = "MARS") #, sub ="threshold = 0.00001 (pre-defined)")

#- - - - - - - - - - - - - - - - - - - - - -
## MaxEnt ####
maxmod <- SDMs[["maxmod_pred"]][[1]]
maxmod_pred <- dismo::predict(maxmod, Env)
maxmod_pred <- fct.rescale(maxmod_pred, x.min = maxmod_pred@data@min, x.max = maxmod_pred@data@max)

# rescale between 0 and 1
maxmod_pred <- fct.rescale(maxmod_pred, x.min = maxmod_pred@data@min, x.max = maxmod_pred@data@max)

plot(maxmod_pred, main = "MaxEnt")
print("Note: If the MaxEnt prediction looks bad, it is maybe because of wrong (or missing) predictors...")

#- - - - - - - - - - - - - - - - - - - - - -
# MaxNet ####
mxnet <- SDMs[["maxnet_pred"]][[1]]
maxnet_pred <- raster::predict(Env, mxnet)
maxnet_pred <- fct.rescale(maxnet_pred, x.min = maxnet_pred@data@min, x.max = maxnet_pred@data@max)
#maxnet_pred <- fct.rescale(maxnet_pred, x.min = maxnet_pred@data@max, x.max = maxnet_pred@data@min)
# Note: if max and min are changed, map fits better to GLM etc. ...

# rescale between 0 and 1
maxnet_pred <- fct.rescale(maxnet_pred, x.min = maxnet_pred@data@min, x.max = maxnet_pred@data@max)

plot(maxnet_pred, main = "MaxNet")

#- - - - - - - - - - - - - - - - - - - - - -
## BRT ####
# extract already calculated predictions
brt_pred <- SDMs[["brt_pred"]][[7]]
#brt_pred

# rescale between 0 and 1
brt_pred <- fct.rescale(brt_pred, x.min = brt_pred@data@min, x.max = brt_pred@data@max)

plot(brt_pred, main = "BRT (GBM)")

## BRT 2 (for ensemble model)
# extract already calculated predictions
brt2_pred <- SDMs[["brt_pred"]][[7]]
#brt_pred

# rescale between 0 and 1
brt2_pred <- fct.rescale(brt2_pred, x.min = brt2_pred@data@min, x.max = brt2_pred@data@max)

plot(brt2_pred, main = "BRT (GBM)")

#- - - - - - - - - - - - - - - - - - - - - -
## XGBoost ####
xgb_pred <- SDMs[["xgb_pred"]][[6]]

names(xgb_pred) <- rownames(Env) #add site names

# rescale between 0 and 1
xgb_pred <- fct.rescale(xgb_pred, x.min = maxmod_pred@data@min, x.max = maxmod_pred@data@max)

plot(xgb_pred, main = "XGBoost")

#- - - - - - - - - - - - - - - - - - - - - -
## RF ####
rf <- SDMs[["rf_pred"]][[1]][,sample(1, 1:ncol(SDMs[["brt_pred"]][[1]]))]
rf_pred <- raster::predict(Env, rf, type = "prob") # prob = continuous prediction
#head(rf_pred)

rf_pred <- as.numeric(rf_pred)
names(rf_pred) <- rownames(Env) #add site names

#- - - - - - - - - - - - - - - - - - - - - -
## RF downsampled ####
rf_downsample <- SDMs[["rf_downsample_pred"]][[1]]
rf_downsample_pred <- predict(rf_downsample, Env, type = "prob")[, "1"] # prob = continuous prediction
#head(rf_downsample_pred)

rf_downsample_pred <- as.numeric(rf_downsample_pred)
names(rf_downsample_pred) <- rownames(Env) #add site names

#- - - - - - - - - - - - - - - - - - - - - -
## SVM ####
svm_e <- SDMs[["svm_pred"]][[1]]
svm_pred <- raster::predict(Env, svm_e, probability = TRUE)
svm_prob <- attr(svm_pred, "probabilities")[,"1"]

# see the first few predictions
#head(svm_prob)

svm_prob <- as.numeric(svm_prob)
names(svm_prob) <- rownames(Env) #add site names

#- - - - - - - - - - - - - - - - - - - - - -
## BIOMOD ####
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

#- - - - - - - - - - - - - - - - - - - - - -
## simple Ensemble model ####
# create Raster stack of the models to be merged
ensm_pred <- raster::stack(gm_pred, maxmod_pred) #lasso_pred, brt2_pred, rf2_pred

# average the predictions
ensm_pred <- calc(ensm_pred, fun = mean, na.rm = T)


#- - - - - - - - - - - - - - - - - - - - - -
## Save maps ####
#- - - - - - - - - - - - - - - - - - - - - -
# stack all predictions
stack_pred <- raster::stack(gm_pred, lm1_pred, lm_subset_pred, lasso_pred, ridge_pred, 
                            mars_pred, maxmod_pred, maxnet_pred, brt_pred, brt2_pred, xgb_pred, svm_pred,
                            rf_pred, rf2_pred, rf_downsample_pred, biomod_pred, ensm_pred)

names(stack_pred) <- c("gm_pred", "lm1_pred", "lm_subset_pred", "lasso_pred", "ridge_pred", 
                 "mars_pred", "maxmod_pred", "maxnet_pred", "brt_pred", "brt2_pred", "xgb_pred", "svm_pred",
                 "rf_pred", "rf2_pred", "rf_downsample_pred", "biomod_pred", "ensm_pred")


# save all predictions
#save(stack_pred, file=paste0(here::here(), "/results/_Maps/SDM_Predictions_", Taxon_name, "_", spID, ".RData"))

## Save best performing model prediction only
# read model performance evaluation table
mod_eval <- read.csv(file=paste0(here::here(), "/results/ModelEvaluation_", Taxon_name, "_", spID, ".csv"))

# select best model based on TSS
best_model <- mod_eval[mod_eval$tss == max(mod_eval$tss), "model"]
best_pred <- stack_pred[[best_model]]

# save best model prediction
save(best_pred, file=paste0(here::here(), "/results/_Maps/SDM_Predictions_", Taxon_name, "_", spID, ".RData"))

plot(best_pred, main = names(best_pred))

# ## -------------------------------------------------------------------------------------- ##
# ## Obtain spatiotemporal projections ----
# ## -------------------------------------------------------------------------------------- ##
# # Goncalas et al. ...
# # https://github.com/joaofgoncalves/GoncalvesAna_et_al_2021/tree/master/RCODE/PostModelAnalyses
# 
# 
# # Models to consider in the ensemble and projection
# modelsToUse <- get_kept_models(myBiomodEM, 1)
# 
# 
# 
# for(projName in projNames){
#   
#   # Obtain spatiotemporal projections
#   myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
#                                     new.env = get(projName),
#                                     proj.name = projName, ## Name of the projection from above variable proj.name
#                                     selected.models = modelsToUse,
#                                     filtered.meth = NULL,
#                                     binary.meth = NULL,
#                                     compress = 'gzip',
#                                     clamping.mask = TRUE,
#                                     output.format = '.grd',
#                                     do.stack = TRUE)
#   
#   
#   # Perform the ensembling of projections
#   myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
#                                            binary.meth = c('TSS','ROC','KAPPA'),
#                                            EM.output = myBiomodEM,
#                                            output.format = '.grd')
#   
#   # Convert all output raster files to GeoTIFF
#   inFolder <- paste(getwd(),"/",sp,"/proj_",projName,sep="")
#   outFolder <- paste(inFolder,"/","GeoTIFF", sep="")
#   dir.create(outFolder)
#   
#   convertToGeoTIFF(inFolder, outFolder)
#   
# } 
# 
# save.image(file=paste(sp,"ModObjects.RData",sep="_"))
