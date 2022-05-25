#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8g")) # expand Java memory
gc()
library(tidyverse)
library(here)

library(CoordinateCleaner)
library(rgbif)
library(taxize)

library(raster)
library(sp)
library(extremevalues)
library(sdmpredictors)

library(biomod2) # also to create pseudo-absences
#library(snowfall) # for parallelization
#library(SSDM)

library(usdm) # for variable inflation factor calculation

library(mgcv) # for GAM
library(gam)  # for GLM (!)
#library(remotes) #to download package from github

## for regularized regressions
# installing the package from github
#remotes::install_github("rvalavi/myspatial")
library(glmnet)
library(myspatial)

library(caret) # for MARS and BRT
library(earth) # for MARS
library(doParallel) # for MARS and XGBoost

library(dismo) # for MaxEnt and BRT
# download maxent.jar 3.3.3k, and place the file in the
# desired folder
# utils::download.file(url = "https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar", 
#     destfile = paste0(system.file("java", package = "dismo"), 
#         "/maxent.jar"), mode = "wb")  ## wb for binary file, otherwise maxent.jar can not execute

library(maxnet) # for MaxNet

library(xgboost) # for XGBoost
library(randomForest) # for RF
library(e1071) # for SVM

library(ROCR) # for Ensemble model
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

## different input data required and package quite newly released...:
#remotes::install_github("peterbat1/fitMaxnet")
#library(fitMaxnet) #for variable importance of maxnet object

# plotting
library(GGally) #for correlations with ggpairs
library(gridExtra)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0("Species_list_", Taxon_name, ".csv"))

corMatPearson <- as.matrix(read.csv(file=paste0("corMatPearson_predictors.csv")))
dimnames(corMatPearson)[[1]] <- dimnames(corMatPearson)[[2]]
# based on Valavi et al. 2021: Pearson 0.8
env_exclude <- caret::findCorrelation(corMatPearson, cutoff = 0.8, names=TRUE)
covarsNames <- dimnames(corMatPearson)[[1]][!(dimnames(corMatPearson)[[1]] %in% env_exclude)]
covarsNames <- covarsNames[covarsNames != "x" & covarsNames != "y"]
# exclude based on VIF
env_vif <- read.csv(file=paste0("VIF_predictors.csv"))
env_exclude <- env_vif %>% filter(is.na(VIF)) %>% dplyr::select(Variables) %>% as.character()
covarsNames <- covarsNames[!(covarsNames %in% env_exclude)]
# excluded:
print("=== We excluded the following variables based on VIF and Pearson correlation: ===")
setdiff(env_vif$Variables, covarsNames)

# final predictor variables
print("=== And we kept the following, final predictor variables: ===")
covarsNames

#- - - - - - - - - - - - - - - - - - - - -
# note: we will load the datasets before each individual model

# load environmental variables (for projections)
Env_norm <- raster::stack(paste0("EnvPredictor_2km_normalized.grd"))
#Env_norm <- stack(Env_norm)

# as dataframe
load(paste0("EnvPredictor_2km_df_normalized.RData")) #Env_norm_df

# define formula for GLM (and biomod)
form <- paste0("occ ~ ", paste0(paste0("s(", covarsNames, ")"), collapse=" + "))


# Calculate the number of cores
no.cores <-  parallel::detectCores()/2 

#- - - - - - - - - - - - - - - - - - - - -
## GAM ####
#- - - - - - - - - - - - - - - - - - - - -
registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
         .export = c("Env_norm", "Env_norm_df", "form"),
         .packages = c("tidyverse", "mgcv", "biomod2", "caret", "gam")) %dopar% { try({

modelName <- "bg.glm"

# identify and load all relevant data files
temp.files <- list.files(path = paste0(here::here(), "/results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# calculate Infinitely Weighted Logistic Regression (IWLR) for GLM and GAM 
# models suggested by Fithian and Hastie (2013)
iwp <- (10^6)^(1-training$occ)

# calculating the case weights (equal weights)
# the order of weights should be the same as presences and backgrounds in the 
# training data
prNum <- as.numeric(table(training$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)

# general settings
tmp <- proc.time()[3]
set.seed(32639)

# setting of family: according to Valavi et al. 2020, but see also following:
# ResearchGate dicussion from 2019: Jincheng Zhou
# "Gaussian family : for continuous decimal data with normal distribution, like weight, length, et al.
# Poisson or quasipoisson family: for positive integer or small natural number like count, individual number, frequency.
# Binomial or quasibinomial family: binary data like 0 and 1, or proportion like survival number vs death number, positive frequency vs negative frequency, winning times vs the number of failtures, et al...
# Gamma family : usually describe time data like the time or duration of the occurrence of the event.
# If your model is overdispered, you should make a correction by choosing quasi-binomial or quasi-poisson"

# run model
gm <- mgcv::gam(formula = as.formula(form),
                data = training,
                family = binomial(link = "logit"),
                weights = wt,
                method = "REML")

# predict (only for model performance)
temp_validation <- mgcv::predict.gam(gm, validation[,colnames(validation) %in% covarsNames], type="response")
#head(temp_validation)

# get running time
temp_model_time <- proc.time()[3] - tmp
temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

temp_runs <- 1

# start time for prediction
tmp <- proc.time()[3]

# predict to whole Europe
temp_prediction <- mgcv::predict.gam(gm, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames], type="response")
temp_prediction <- as.numeric(temp_prediction)
names(temp_prediction) <- rownames(Env_norm_df) #add site names
temp_prediction <- as.data.frame(temp_prediction)

temp_prediction$x <- Env_norm_df$x
temp_prediction$y <- Env_norm_df$y

colnames(temp_prediction)[1] <- "layer"

# get running time
temp_predict_time <- proc.time()[3] - tmp
temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

# varImp
temp_varImp <- caret::varImp(gm, scale=T) #scaled between 0 and 100% 
# check if it looks right
#print(temp_varImp) 
# replace Inf values with 0
temp_varImp[is.infinite(temp_varImp[,1]),1] <- 0   
temp_varImp$Predictor <- rownames(temp_varImp)
colnames(temp_varImp)[1] <- "gm"

gm_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp)
save(gm_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_gm_", spID, ".RData"))

rm(gm, gm_list, modelName, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp)

})}
stopImplicitCluster()


registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
         .export = c("Env_norm", "Env_norm_df", "form"),
         .packages = c("tidyverse", "mgcv", "biomod2", "caret", "gam")) %dopar% { try({

#- - - - - - - - - - - - - - - - - - - - -
## GLM ####
#- - - - - - - - - - - - - - - - - - - - -
modelName <- "bg.glm"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# calculating the weights
# the order of weights should be the same as presences and backgrounds in the 
# training data
prNum <- as.numeric(table(training$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)

tmp <- proc.time()[3]

# the base glm model with linear terms
lm1 <- glm(occ ~., data = training, weights = wt, 
           family = binomial(link = "logit"))
#summary(lm1)

temp_model_time <- proc.time()[3] - tmp
temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

temp_validation <- stats::predict.glm(lm1, validation[,colnames(validation) %in% covarsNames], type="response")
#head(temp_validation)

temp_runs <- 1

# set start time for prediction
tmp <- proc.time()[3]

# predict to whole environment
temp_prediction <- stats::predict.glm(lm1, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames], type="response")

temp_prediction <- as.numeric(temp_prediction)
names(temp_prediction) <- rownames(Env_norm_df) #add site names
temp_prediction <- as.data.frame(temp_prediction)

temp_prediction$x <- Env_norm_df$x
temp_prediction$y <- Env_norm_df$y

colnames(temp_prediction)[1] <- "layer"

# get running time
temp_predict_time <- proc.time()[3] - tmp
temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

# varImp
temp_varImp <- caret::varImp(lm1, scale=T) #scaled between 0 and 100% 
# check if it looks right
#print(temp_varImp) 
# replace Inf values with 0
temp_varImp[is.infinite(temp_varImp[,1]),1] <- 0   
temp_varImp$Predictor <- rownames(temp_varImp)
colnames(temp_varImp)[1] <- "lm1"

lm1_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp)
save(lm1_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_lm1_", spID, ".RData"))

rm(lm1_list, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp)


# load library again to make sure that function works correctly
library(gam)

# model scope for subset selection
mod_scope <- gam::gam.scope(frame = training, form=FALSE, #form=T: formular, else character vector
                            response = 1) #column with response variable

# Note: alternatively, you can define one by your own
# mod_scope <- list("cti" = ~1 + cti + poly(cti, 2),
#                   "disturb" = ~1 + disturb + poly(disturb, 2),
#                   "mi" = ~1 + mi + poly(mi, 2),
#                   "rainann" = ~1 + rainann + poly(rainann, 2),
#                   "raindq" = ~1 + raindq + poly(raindq, 2),
#                   "rugged" = ~1 + rugged + poly(rugged, 2),
#                   "soildepth" = ~1 + soildepth + poly(soildepth, 2),
#                   "soilfert" = ~1 + soilfert + poly(soilfert, 2),
#                   "solrad" = ~1 + solrad + poly(solrad, 2),
#                   "tempann" = ~1 + tempann + poly(tempann, 2),
#                   "topo" = ~1 + topo + poly(topo, 2),
#                   "vegsys" = ~1 + vegsys)

tmp <- proc.time()[3]
set.seed(32639)

lm_subset <- gam::step.Gam(object = lm1,
                           scope = mod_scope,
                           direction = "both", # up and down-moving of terms
                           data = training, # this is optional
                           parallel = FALSE,   # optional
                           trace = FALSE) #if information is printed or not

temp_validation <- stats::predict.glm(lm_subset, validation[,colnames(validation) %in% covarsNames], type="response")
#head(temp_validation)

temp_model_time <- proc.time()[3] - tmp
temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

#summary(lm_subset)

temp_runs <- 1

# start time for prediction
tmp <- proc.time()[3]

temp_prediction <- stats::predict.glm(lm_subset, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames], type="response")

temp_prediction <- as.numeric(temp_prediction)
names(temp_prediction) <- rownames(Env_norm_df) #add site names
temp_prediction <- as.data.frame(temp_prediction)

temp_prediction$x <- Env_norm_df$x
temp_prediction$y <- Env_norm_df$y

colnames(temp_prediction)[1] <- "layer"

# get running time
temp_predict_time <- proc.time()[3] - tmp
temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

# varImp
temp_varImp <- caret::varImp(lm_subset, scale=T) #scaled between 0 and 100% 
# check if it looks right
#print(temp_varImp) 
# replace Inf values with 0
temp_varImp[is.infinite(temp_varImp[,1]),1] <- 0   
temp_varImp$Predictor <- rownames(temp_varImp)
colnames(temp_varImp)[1] <- "lm_subset"

lm_subset_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp)
save(lm_subset_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_lm_subset_", spID, ".RData"))

rm(lm1, lm_subset, lm_subset_list, modelName, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp)

})}
stopImplicitCluster()


registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
         .export = c("Env_norm", "Env_norm_df", "form"),
         .packages = c("tidyverse", "mgcv", "biomod2", "caret", "gam", "glmnet", "myspatial")) %dopar% { try({

#- - - - - - - - - - - - - - - - - - - - -
## Regularized regressions: LASSO and RIDGE regression ####
##- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.lasso"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# generating the quadratic terms for all continuous variables
# function to creat quadratic terms for lasso and ridge
quad_obj <- make_quadratic(training, cols = covarsNames)

# now we can predict this quadratic object on the training, testing, and prediction data
# this make two columns for each covariates used in the transformation
training_quad <- predict(quad_obj, newdata = training)
testing_quad <- predict(quad_obj, newdata = validation[,colnames(validation) %in% covarsNames])
predicting_quad <- predict(quad_obj, newdata = Env_norm_df[,colnames(Env_norm_df) %in% covarsNames])

# Define vector with appropriate covarsNames for sparse matrix.
# More specifically: create names like this: bio1_1, bio1_2, bio2_1, bio2_2...
sparse_covarsNames <- paste0(rep(covarsNames, each=2), "_", 1:2)

# convert the data.frames to sparse matrices
# select all quadratic (and non-quadratic) columns, except the y (occ)
new_vars <- names(training_quad)[names(training_quad) != "occ"]
training_sparse <- Matrix::sparse.model.matrix(~. -1, training_quad[, new_vars])
testing_sparse <- Matrix::sparse.model.matrix( ~. -1, testing_quad[, new_vars])
predicting_sparse <- Matrix::sparse.model.matrix( ~. -1, predicting_quad[, sparse_covarsNames])

# calculating the case weights
prNum <- as.numeric(table(training_quad$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training_quad$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)

#- - - - - - - - - - - - - - - - - - - - -
## fitting lasso
# Note: The alpha parameter in this model ranges from 0 to 1, where selecting an 
# alpha of 0 leads to ridge regression and 1 to lasso and anything in between 
# is a combination of both called elastic-net.

# The lambda parameter controls regularization ? it is the penalty applied 
# to the model?s coefficients. To select the best lambda, internal cross-
# validation was used (in cv.glmnet function).
tmp <- proc.time()[3]
set.seed(32639)
lasso_cv <- glmnet::cv.glmnet(x = training_sparse,
                              y = training_quad$occ,
                              family = "binomial",
                              alpha = 1, # fitting lasso
                              weights = wt,
                              nfolds = 10) # number of folds for cross-validation

# predicting with lasso model
temp_validation <- predict(lasso_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
#head(temp_validation)

temp_model_time <- proc.time()[3] - tmp
temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

# the cross-validation result
#plot(lasso_cv)
# Note: dashed lines show options for best Lambda parameter
# left: Lambda with minimum deviance (lambda.min)
# right: best Lambda within one standard deviation of left-dashed line (l...1se)

#print("One of the two Lambda (dashed lines) needs to be selected for prediction.")

tmp <- proc.time()[3]

# predict for whole environment
#temp_prediction <- glmnet:::predict.glmnet(object = lasso_cv, newx = predicting_sparse, type = "response", s = "lambda.min")[,1]
temp_prediction <- predict(lasso_cv, predicting_sparse, type = "response", s = "lambda.min")[,1]
#head(temp_prediction)
temp_prediction <- data.frame("site" = names(temp_prediction), "layer" = as.numeric(temp_prediction)) %>%
  full_join(Env_norm_df %>% mutate("site" = rownames(Env_norm_df)) %>% dplyr::select(x,y,site), by = "site") %>%
  dplyr::select(-site)

temp_runs <- 1

# get running time
temp_predict_time <- proc.time()[3] - tmp
temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

# varImp
coefs <- coef(lasso_cv) #!!! simple way...
# remove first element (Intercept)
temp_varImp <- data.frame("Predictor" = coefs@Dimnames[[1]][coefs@i + 1], "lasso" = coefs@x)[-1,] 
colnames(temp_varImp)[2] <- "lasso"
temp_varImp$Predictor <- stringr::str_split_fixed(temp_varImp$Predictor, "_[:digit:]", 2)[,1]
temp_varImp <- temp_varImp %>% group_by(Predictor) %>% summarize_all(mean, na.rm=T) %>% as.data.frame()
    
# check if it looks right
#print(temp.vi)

lasso_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp)
save(lasso_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_lasso_", spID, ".RData"))

rm(lasso_cv, lasso_list, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp, coefs)

#- - - - - - - - - - - - - - - - - - - - -
## fitting ridge resgression (alpha=0) while identify the right lambda
tmp <- proc.time()[3]
set.seed(32639)
ridge_cv <- glmnet::cv.glmnet(x = training_sparse,
                              y = training_quad$occ,
                              family = "binomial",
                              alpha = 0, # fitting ridge
                              weights = wt,
                              nfolds = 10) # number of folds for cross-validation
# predict ridge
temp_validation <- predict(ridge_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
#head(temp_validation)

temp_model_time <- proc.time()[3] - tmp
temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

#plot(ridge_cv)

temp_runs <- 1

tmp <- proc.time()[3]

# predict for whole environment
temp_prediction <- predict(ridge_cv, predicting_sparse, type = "response", s = "lambda.min")[,1]
#head(ridge_prediction)
temp_prediction <- data.frame("site" = names(temp_prediction), "layer" = as.numeric(temp_prediction)) %>%
  full_join(Env_norm_df %>% mutate("site" = rownames(Env_norm_df)) %>% dplyr::select(x,y,site), by = "site") %>%
  dplyr::select(-site)

# get running time
temp_predict_time <- proc.time()[3] - tmp
temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

coefs <- coef(ridge_cv) #!!! simple way...
# remove first element (Intercept)
temp_varImp <- data.frame("Predictor" = coefs@Dimnames[[1]][coefs@i + 1], "ridge" = coefs@x)[-1,] 
colnames(temp_varImp)[2] <- "ridge"
temp_varImp$Predictor <- stringr::str_split_fixed(temp_varImp$Predictor, "_[:digit:]", 2)[,1]
temp_varImp <- temp_varImp %>% group_by(Predictor) %>% summarize_all(mean, na.rm=T) %>% as.data.frame()

ridge_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp)
save(ridge_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_ridge_", spID, ".RData"))

rm(ridge_cv, ridge_list, modelName, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction)
rm(testing_sparse, training_sparse, new_vars, predicting_sparse, sparse_covarsNames, testing_quad, training_quad, predicting_quad, quad_obj, temp_varImp, coefs)
gc()

})}
stopImplicitCluster()


registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
         .export = c("Env_norm", "Env_norm_df", "form"),
         .packages = c("tidyverse","biomod2", "caret", "earth", "doParallel")) %dopar% { try({

#- - - - - - - - - - - - - - - - - - - - -
## MARS ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.mars"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name),
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
#lapply(temp.files, load, .GlobalEnv)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/2

mars_pred_list <- list()

for(no.runs in 1:no.loop.runs){
  
  temp.files.subset <- list.files(path = paste0("./results/",Taxon_name), 
                                  pattern = paste0(modelName, no.runs, "_[[:graph:]]*", spID), full.name = T)
  
  lapply(temp.files.subset, load, .GlobalEnv)
  
  # change the response to factor variables with non-numeric levels
  training$occ <- as.factor(training$occ)
  levels(training$occ) <- c("C0", "C1")
  
  # settings
  mytuneGrid <- expand.grid(nprune = 1+length(covarsNames), # max. number of terms including intercept
                            degree = 2) # no interaction = 1 = builds additive model; maximum degree of interaction (Friedman's mi)
  mytrControl <- trainControl(method = "cv", # cross validation
                              number = 5, # 5-fold cross-validation
                              classProbs = TRUE, #calculates class probabilities for held-out samples
                              summaryFunction = twoClassSummary, #calculate alternate performance summaries: sensitivity, specificity and area under the ROC curve
                              allowParallel = TRUE) # allow parallel
  tmp <- proc.time()[3]
  
  # parallelize
  cluster <- makeCluster(2) # you can use all cores of your machine instead e.g. 8
  registerDoParallel(cluster)
  
  set.seed(32639)
  mars_fit <- caret::train(form = occ ~ .,
                           data = training,
                           method = "earth", # = multivariate adaptive regression spline (MARS)
                           metric = "ROC", # select best model based on ROC
                           trControl = mytrControl,
                           tuneGrid = mytuneGrid,
                           thresh = 0.00001) # forward stepping threshold, default is 0.001
  stopCluster(cluster)
  registerDoSEQ()
  
  temp_validation <- caret::predict.train(mars_fit, validation[,colnames(validation) %in% covarsNames], type="prob")
  temp_validation <- data.frame("sites"=rownames(validation[,colnames(validation) %in% covarsNames]), "occ"=temp_validation[,"C1"])
  #head(temp_validation)
  
  temp_model_time <- proc.time()[3] - tmp
  temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))
  #plot(mars_fit)
  
  # caluclate variable importance (for later)
  mars_varImp <- caret::varImp(mars_fit, scale=T) #scaled between 0 and 100% 
  mars_varImp <- data.frame("importance" = mars_varImp$importance, "Predictor"=rownames(mars_varImp$importance))
  
  tmp <- proc.time()[3]

  # create raster layer of predictions for whole environmental space
  temp_prediction <- caret::predict.train(mars_fit, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames], type="prob")[,"C1"]
  temp_prediction <- as.numeric(temp_prediction)
  # add names of grid cell (only for those that have no NA in any layer)
  names(temp_prediction) <- rownames(Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),])
  temp_prediction <- as.data.frame(temp_prediction)
  temp_prediction$x <- Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),]$x
  temp_prediction$y <- Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),]$y
  temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
    rename("layer" = temp_prediction)
  
  # get running time
  temp_predict_time <- proc.time()[3] - tmp
  temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

  mars_pred_list[[no.runs]] <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, validation=temp_validation, prediction=temp_prediction, varImp=mars_varImp)

  print(paste0("Run no. ", no.runs, " of ", no.loop.runs ," is ready now."))
  rm(temp_prediction, mars_varImp, mars_fit, temp_validation, temp_model_time, temp_predict_time)
}

# merge all MARS predictions
temp_validation <- do.call(rbind, lapply(mars_pred_list, "[[", 4))
temp_validation <- temp_validation %>% group_by(sites) %>% summarise(across(everything(), mean))

temp_model_time <- c(sum(as.numeric(sapply(mars_pred_list, "[[", 2)[1,]), na.rm=T), mars_pred_list[[1]][[2]][2])
temp_predict_time <- c(sum(as.numeric(sapply(mars_pred_list, "[[", 3)[1,]), na.rm=T), mars_pred_list[[1]][[3]][2])

temp_varImp <- do.call(rbind, lapply(mars_pred_list, "[[", 6)) %>% 
  group_by(Predictor) %>% summarise_all(mean, na.rm=T)
temp_varImp <- as.data.frame(temp_varImp)
colnames(temp_varImp)[2] <- "mars"
#print(temp_varImp)

#temp_models <- sapply(mars_pred_list, "[[", 1)

temp_runs <- length(mars_pred_list)

## extract predicted probabilities (for later)
temp_prediction <- do.call(rbind, lapply(mars_pred_list, "[[", 5)) %>% 
  group_by(x,y) %>% summarise_all(mean, na.rm=T)
temp_uncertainty <- do.call(rbind, lapply(mars_pred_list, "[[", 5)) %>% 
  group_by(x,y) %>% summarise_all(sd, na.rm=T)

mars_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp, uncertainty=temp_uncertainty)
save(mars_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_mars_", spID, ".RData"))

rm(mars_list, modelName, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_uncertainty, temp_varImp)
gc()

})}
stopImplicitCluster()


for(spID in unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID)){ try({

#- - - - - - - - - - - - - - - - - - - - -
## MaxEnt and MaxNet ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.maxent"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# "We used five different regularization multipliers (0.5, 1, 2, 3 and 4)
# in combination with different features (L, LQ, H, LQH, LQHP) to find the 
# best parameters that maximizes the average AUC-ROC in CV."

# function for simultaneous tuning of maxent regularization multiplier and features
maxent_param <- function(data, y = "occ", k = 5, folds = NULL, filepath){
  require(dismo)
  require(caret)
  require(precrec)
  if(is.null(folds)){
    # generate balanced CV folds
    folds <- caret::createFolds(y = as.factor(data$occ), k = k)
  }
  names(data)[which(names(data) == y)] <- "occ"
  
  # regularization multipliers
  ms <- c(0.5, 1, 2, 3, 4)
  grid <- expand.grid(
    regmult = paste0("betamultiplier=", ms),
    
    # features
    features = list(
      c("noautofeature", "nothreshold"), # LQHP
      c("noautofeature", "nothreshold", "noproduct"), # LQH
      c("noautofeature", "nothreshold", "nohinge", "noproduct"), # LQ
      c("noautofeature", "nothreshold", "nolinear", "noquadratic", "noproduct"), # H
      c("noautofeature", "nothreshold", "noquadratic", "nohinge", "noproduct")), # L
    stringsAsFactors = FALSE
  )
  AUCs <- c()
  for(n in seq_along(grid[,1])){
    full_pred <- data.frame()
    for(i in seq_len(length(folds))){
      trainSet <- unlist(folds[-i])
      testSet <- unlist(folds[i])
      if(inherits(try(
        maxmod <- dismo::maxent(x = data[trainSet, covarsNames],
                                p = data$occ[trainSet],
                                removeDuplicates = FALSE,
                                path = filepath,
                                args = as.character(unlist(grid[n, ])))
      ), "try-error")){
        next
      }
      modpred <- dismo::predict(maxmod, data[testSet, covarsNames], type = c("cloglog")) #, args = "outputformat=cloglog"
      pred_df <- data.frame(score = modpred, label = data$occ[testSet])
      full_pred <- rbind(full_pred, pred_df)
    }
    AUCs[n] <- precrec::auc(precrec::evalmod(scores = full_pred$score,
                                             labels = full_pred$label))[1,4]
  }
  best_param <- as.character(unlist(grid[which.max(AUCs), ]))
  return(best_param)
}

## now use the function to tune MaxEnt
# number of folds
nfolds <- ifelse(sum(as.numeric(training$occ)) < 10, 2, 5)

tmp <- proc.time()[3]
set.seed(32639)

# tune maxent parameters
param_optim <- maxent_param(data = training,
                            k = nfolds,
                            filepath = paste0("results/", Taxon_name, "/maxent_files"))

## fit a maxent model with the tuned parameters
maxent <- dismo::maxent(x = training[, covarsNames],
                        p = training$occ,
                        removeDuplicates = FALSE, #remove occurrences that fall into same grid cell (not necessary)
                        #path = paste0("results/", Taxon_name, "/maxent_files"), #wanna save files?
                        args = param_optim)

#temp_validation <- dismo::predict(maxent, validation[,colnames(validation) %in% covarsNames], type="response")
temp_validation <- dismo::predict(maxent, validation[,colnames(validation) %in% covarsNames], type = c("cloglog"))
temp_validation <- as.numeric(temp_validation)
names(temp_validation) <- rownames(validation[,colnames(validation) %in% covarsNames]) #add site names
#head(temp_validation)

temp_model_time <- proc.time()[3] - tmp
temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

tmp <- proc.time()[3]
# create raster layer of predictions for whole environmental space
#temp_prediction <- raster::predict(Env_norm, maxent)
#temp_prediction <- data.frame(raster::rasterToPoints(temp_prediction))
gc()
#temp_prediction <- dismo::predict(maxent, Env_norm_df %>% dplyr::select(-x, -y)) # Java out of memory
temp_prediction <- dismo::predict(maxent, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames], type = c("cloglog"))
temp_prediction <- as.numeric(temp_prediction)
names(temp_prediction) <- rownames(Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),]) #add site names
temp_prediction <- as.data.frame(temp_prediction)
temp_prediction$x <- Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),]$x
temp_prediction$y <- Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),]$y
colnames(temp_prediction)[1] <- "layer"

temp_runs <- 1

temp_predict_time <- proc.time()[3] - tmp
temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

# varImp
temp_varImp <- as.data.frame(temp.results[str_detect(rownames(maxent@results),"permutation.importance"),])
# Note: permutation importance = determine the importance of predictors 
# calculated by permuting values of each predictor &  resulting reduction
# in training AUC: large reduction = model is influenced by that predictor. 
    
# extract predictor names
temp_varImp$Predictor <- stringr::str_split_fixed(rownames(temp_varImp), "[.]perm", 2)[,1]
colnames(temp_varImp) <- c("maxent", "Predictor")  
# plot variable importance
#plot(maxent)   
# plot response curve
#response(maxent)  

maxent_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp)
save(maxent_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_maxent_", spID, ".RData"))

rm(maxent, maxent_list, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp)
gc()

## MaxNet
presences <- training$occ # presence (1s) and background (0s) points
covariates <- training[, covarsNames] # predictor covariates

tmp <- proc.time()[3]
set.seed(32639)

maxnet <- maxnet::maxnet(p = presences,
                        data = covariates,
                        regmult = as.numeric(stringr::str_split(param_optim[1], "=")[[1]][2]), # regularization multiplier, a constant, taken from maxmod
                        maxnet::maxnet.formula(presences, covariates, classes = "default"))

# predicting with MaxNet
temp_validation <- predict(maxnet, validation[,colnames(validation) %in% covarsNames], type = c("cloglog"))[, 1]
#head(temp_validation)

temp_model_time <- proc.time()[3] - tmp
temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

tmp <- proc.time()[3]
# create raster layer of predictions for whole environmental space
gc()
#temp_prediction <- raster::predict(Env_norm, maxnet)
#temp_prediction <- data.frame(rasterToPoints(temp_prediction))
temp_prediction <- predict(maxnet, Env_norm_df[,colnames(Env_norm_df) %in% names(covariates)], type = "cloglog")#[, 1]
temp_prediction <- data.frame("sites"=rownames(temp_prediction), "layer"=as.numeric(temp_prediction))
temp_env <- Env_norm_df %>% mutate("sites"=rownames(Env_norm_df))
temp_prediction <- temp_prediction %>% dplyr::full_join(temp_env[,c("sites", "x", "y")], by=c("sites")) %>% dplyr::select(-sites) 
colnames(temp_prediction)[1] <- "layer"

temp_predict_time <- proc.time()[3] - tmp
temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

temp_runs <- 1

# varImp
# use package designed to calculate variable importance of maxnet object
#fitMaxnet::varImportance(maxnet)
# instead: extract beta coefficients
temp_varImp <- maxnet$betas
temp_varImp <- as.data.frame(temp_varImp[str_detect(names(temp_varImp),"hinge")]) 
temp_varImp$Predictor <- stringr::str_split_fixed(rownames(temp_varImp), "[()]", 3)[,2]
temp_varImp <- temp_varImp %>% group_by(Predictor) %>% summarize_all(mean, na.rm=T) %>% as.data.frame()
colnames(temp_varImp) <- c("Predictor", "maxnet")
# plot response curves
#plot(maxnet)

maxnet_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp)
save(maxnet_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_maxnet_", spID, ".RData"))

rm(temp_env, maxnet, maxnet_list, modelName, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, presences, covariates, temp_varImp)
gc()

})}
rm(training, validation)

registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
         .export = c("Env_norm", "Env_norm_df", "form"),
         .packages = c("tidyverse", "biomod2", "caret", "dismo")) %dopar% { try({

#- - - - - - - - - - - - - - - - - - - - -
## BRT (GBM) ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/2

brt_pred_list <- list()
temp_varImp <- vector(mode="list", length=no.loop.runs)

for(no.runs in 1:no.loop.runs){
  
  temp.files.subset <- list.files(path = paste0("./results/",Taxon_name), 
                                  pattern = paste0(modelName, no.runs, "_[[:graph:]]*", spID), full.name = T)
  
  lapply(temp.files.subset, load, .GlobalEnv)
  
  # calculating the case weights
  prNum <- as.numeric(table(training$occ)["1"]) # number of presences
  bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
  wt <- ifelse(training$occ == 1, 1, prNum / bgNum)
  set.seed(32639)
  
  # parameter settings
  brt <- NULL
  ntrees <- 50
  tcomplexity <- ifelse(prNum < 50, 1, 5)
  lrate <- 0.001
  m <- 0
  
  tmp <- proc.time()[3]

  # start modeling
  while(is.null(brt)){
    m <- m + 1
    if(m < 11){
      ntrees <- 50
      lrate <- 0.001
    } else if(m < 21){
      lrate <- 0.0001
    } else if(m < 31){
      ntrees <- 25
      lrate <- 0.0001
    } else if(m < 41){
      ntrees <- 25
      lrate <- 0.0001
      tcomplexity <- ifelse(prNum < 50, 1, 3)
    } else{
      break
    }
    
    if(inherits(try(
      brt <- dismo::gbm.step(data = training,
                             gbm.x = 2:ncol(training), # column indices for covariates
                             gbm.y = 1, # column index for response
                             family = "bernoulli", # = binomial
                             tree.complexity = tcomplexity,
                             learning.rate = lrate, # weight for individual trees
                             bag.fraction = 0.75, # proportion of observations used in selecting variables
                             max.trees = 10000,
                             n.trees = ntrees,
                             n.folds = 5, # 5-fold cross-validation
                             site.weights = wt,
				     plot.main = FALSE, plot.folds = FALSE, #stop plotting
                             silent = TRUE) # avoid printing the cv results
    ), "try-error")){
      cat("Couldn't fit model", spID, "in the iteration", m, "\n")
    } 
    
  }
  if(is.null(brt)){
    next
  }
  
  # Note: model tuning with ~50,000 background points may take ~1h.
  
  temp_model_time <- proc.time()[3] - tmp
  temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

  # the optimal number of trees (intersect of green and red line in plot)
  #brt$gbm.call$best.trees
  
  # get interactions of predictors
  temp.find.int <- dismo::gbm.interactions(brt)
  
  # predicting with the best trees
  temp_validation <- dismo::predict(brt, validation[,colnames(validation) %in% covarsNames], n.trees = brt$gbm.call$best.trees, type = "response")
  temp_validation <- data.frame("sites"=rownames(validation[,colnames(validation) %in% covarsNames]), "occ"=temp_validation)
  #head(temp_validation)

  tmp <- proc.time()[3]
  # create raster layer of predictions for whole environmental space
  temp_prediction <- dismo::predict(brt, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames], n.trees = brt$gbm.call$best.trees, type = "response")
  temp_prediction <- as.numeric(temp_prediction)
  # add names of grid cell (only for those that have no NA in any layer)
  names(temp_prediction) <- rownames(Env_norm_df)
  temp_prediction <- as.data.frame(temp_prediction)
  temp_prediction$x <- Env_norm_df$x
  temp_prediction$y <- Env_norm_df$y
  temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
    rename("layer" = temp_prediction)
  
  temp_predict_time <- proc.time()[3] - tmp
  temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

  # varImp
  temp_temp_varImp <- brt$contributions
  # check if it looks right
  #print(temp_temp_varImp)
  temp_temp_varImp <- temp_temp_varImp %>% mutate("Predictor" =rownames(temp_temp_varImp))
  temp_temp_varImp <- temp_temp_varImp[-1]
  colnames(temp_temp_varImp)[1] <- "brt"
  temp_varImp[[no.runs]] <- temp_temp_varImp

  brt_pred_list[[no.runs]] <- list(time_model=temp_model_time, time_predict=temp_predict_time, validation=temp_validation, prediction=temp_prediction, settings=data.frame("n.trees"=ntrees, "l.rate"=lrate, "t.complexity"=tcomplexity), 
                                   interactions=temp.find.int$interactions)
  rm(temp_temp_varImp, brt, temp_prediction, temp_validation, temp.find.int, temp_model_time, temp_predict_time)
}

# merge all BRT predictions
temp_validation <- do.call(rbind, lapply(brt_pred_list, "[[", 3))
temp_validation <- temp_validation %>% group_by(sites) %>% summarise(across(everything(), mean))

# varImp
temp_varImp <- do.call(rbind, temp_varImp) %>% 
      group_by(Predictor) %>% summarise_all(mean, na.rm=T) %>%
      as.data.frame()

# sum of computing time
temp_model_time <- c(sum(as.numeric(sapply(brt_pred_list, "[[", 1)[1,]), na.rm=T), brt_pred_list[[1]][[1]][2])
temp_predict_time <- c(sum(as.numeric(sapply(brt_pred_list, "[[", 2)[1,]), na.rm=T), brt_pred_list[[1]][[2]][2])

temp_settings <- rowMeans(as.data.frame(unlist(sapply(brt_pred_list, "[[", 5))), na.rm=T)
temp_settings <- data.frame("parameter"=names(brt_pred_list[[1]][[5]]), "setting"=temp_settings)

temp_runs <- length(brt_pred_list)

# calculate mean of tuned number of trees (if necessary)
if(ncol(sapply(brt_pred_list, "[[", 6))!=1) {
  temp_interactions <- rowMeans(sapply(brt_pred_list, "[[", 6), na.rm=T)
}else{
  temp_interactions <- sapply(brt_pred_list, "[[", 6)
}

temp_interactions <- as.data.frame(matrix(temp_interactions, ncol=length(covarsNames)))
rownames(temp_interactions) <- covarsNames
colnames(temp_interactions) <- covarsNames

## extract predicted probabilities (for later)
temp_prediction <- do.call(rbind, lapply(brt_pred_list, "[[", 4)) %>% 
  group_by(x,y) %>% summarise_all(mean, na.rm=T)
temp_uncertainty <- do.call(rbind, lapply(brt_pred_list, "[[", 4)) %>%
  group_by(x,y) %>% summarise_all(function(x) sd(x, na.rm=T))

#temp.models <- sapply(brt_pred_list, "[[", 1)

temp_runs <- 1

brt_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, uncertainty=temp_uncertainty, settings=temp_settings, interactions=temp_interactions, varImp=temp_varImp)
save(brt_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_brt_", spID, ".RData"))

rm(brt_pred_list, brt_list, modelName, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_uncertainty, temp_interactions, temp_settings, temp_varImp)

#- - - - - - - - - - - - - - - - 
## Model BRT2 for ensemble modelling with consistent background data (bg.glm)
modelName <- "bg.glm"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# calculating the case weights
prNum <- as.numeric(table(training$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)
set.seed(32639)

# settings
brt2 <- NULL
ntrees <- 50
tcomplexity <- ifelse(prNum < 50, 1, 5)
lrate <- 0.001
m <- 0

tmp <- proc.time()[3]

# start modeling! We use the "try" notation so if a species fails to fit, the loop will continue.
while(is.null(brt2)){
  m <- m + 1
  if(m < 11){
    ntrees <- 50
    lrate <- 0.001
  } else if(m < 21){
    lrate <- 0.0001
  } else if(m < 31){
    ntrees <- 25
    lrate <- 0.0001
  } else if(m < 41){
    ntrees <- 25
    lrate <- 0.0001
    tcomplexity <- ifelse(prNum < 50, 1, 3)
  } else{
    break
  }
 
  
  if(inherits(try(
    
    brt2 <- dismo::gbm.step(data = training,
                            gbm.x = 2:ncol(training), # column indices for covariates
                            gbm.y = 1, # column index for response
                            family = "bernoulli",
                            tree.complexity = tcomplexity,
                            learning.rate = lrate,
                            bag.fraction = 0.75,
                            max.trees = 10000,
                            n.trees = ntrees,
                            n.folds = 5, # 5-fold cross-validation
                            site.weights = wt,
				    plot.main = FALSE, plot.folds = FALSE, #stop plotting
                            silent = TRUE) # avoid printing the cv results
  ), "try-error")){
    cat("Couldn't fit model", spID, "in the iteration", m, "\n")
  } 
  
}
if(is.null(brt2)){
  next
}

# Note: model tuning with ~50,000 background points may take ~1h (or less!).

# the optimal number of trees (intersect of green and red line in plot)
#brt2$gbm.call$best.trees

# get interactions of predictors
temp.find.int <- gbm.interactions(brt2)

# predicting with the best trees
temp_validation <- dismo::predict(brt2, validation[,colnames(validation) %in% covarsNames], n.trees = brt2$gbm.call$best.trees, type = "response")
#head(temp_validation)
temp_validation <- as.numeric(temp_validation)
names(temp_validation) <- rownames(validation[,colnames(validation) %in% covarsNames]) #add site names

temp_model_time <- proc.time()[3] - tmp
  temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

tmp <- proc.time()[3]
# predict to raster (for later)
temp_prediction <- dismo::predict(brt2, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames], n.trees = brt2$gbm.call$best.trees, type = "response") #maybe remove latter part
temp_prediction <- as.numeric(temp_prediction)
# add names of grid cell (only for those that have no NA in any layer)
names(temp_prediction) <- rownames(Env_norm_df)
temp_prediction <- as.data.frame(temp_prediction)
temp_prediction$x <- Env_norm_df$x
temp_prediction$y <- Env_norm_df$y
temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
  rename("layer" = temp_prediction)

temp_predict_time <- proc.time()[3] - tmp
temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

temp_settings <- data.frame("parameter"=c("n.trees", "l.rate", "t.complexity"), "setting"=c(ntrees,lrate,tcomplexity))
temp_runs <- 1

# extract VI
temp_varImp <- brt2$contributions
colnames(temp_varImp) <- c("Predictor", "brt2")
# check if it looks right
#print(temp_varImp)

brt2_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, settings=temp_settings, interactions=temp.find.int$interactions, varImp=temp_varImp)
save(brt2_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_brt2_", spID, ".RData"))

rm(brt2, brt2_list, modelName, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp.find.int, temp_settings, temp_varImp)

})}
stopImplicitCluster()


registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
        .export = c("Env_norm", "Env_norm_df", "form"),
        .packages = c("tidyverse", "caret", "dismo", "xgboost")) %dopar% { try({
          
  print(paste0("Species ", spID, ", model building for XGBoost."))
  
#- - - - - - - - - - - - - - - - - - - - -
## XGBoost ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/2

xgb_pred_list <- list()

for(no.runs in 1:no.loop.runs){
  
  temp.files.subset <- list.files(path = paste0(here::here(), "/results/",Taxon_name), 
                                  pattern = paste0(modelName, no.runs, "_[[:graph:]]*", spID), full.name = T)
  
  lapply(temp.files.subset, load, .GlobalEnv)
  
  # start modeling! We use the "try" notation so if a species fails to fit, the loop will continue.
  xgb_fit <- NULL
  # specify num of folds for cross-validation
  k <- ifelse(sum(training$occ, na.rm=T) < 10, 2, 5)
  
  tmp <- proc.time()[3]

  # change the response to factor variables with non-numeric levels
  training$occ <- as.factor(training$occ)
  levels(training$occ) <- c("C0", "C1") # caret does not accept 0 and 1 as factor levels
  
  # train control for cross-validation for model tuning
  mytrControl <- caret::trainControl(method = "cv",
                              number = k, # 5 fold cross-validation
                              summaryFunction = twoClassSummary,
                              classProbs = TRUE,
                              allowParallel = TRUE)
  # setting the range of parameters for grid search tuning
  mytuneGrid <- expand.grid(
    nrounds = seq(from = 500, to = 15000, by = 500), # number of boosting iterations
    eta = 0.001,       # shrinkage = learning rate, low = robust to overfitting but slower
    max_depth = 5,     # max. tree depth, default = 6
    subsample = 0.75,  # Subsample Percentage: ratio of training instance -> 75% of data used for fitting, the less the faster
    gamma = 0,         # max. loss reduction required to make further partition on leaf node of tree, the larger the more conservative
    colsample_bytree = 0.8, # Subsample Ratio of Columns when constructing each tree
    min_child_weight = 1    # Minimum Sum of Instance Weight (hessian) needed in a child, 1 = default; if partition step results in leaf node with sum of instance weight < min_child_weight then stop partitioning, the larger the more conservative
  )
  
  cluster <- makeCluster(no.cores) # you can use all cores of your machine instead e.g. 8
  registerDoParallel(cluster)
  set.seed(32639)
  if(inherits(try(
    xgb_fit <- caret::train(form = occ ~ .,
                            data = training,
                            method = "xgbTree",
                            trControl = mytrControl,
                            metric ="ROC", 
                            tuneGrid = mytuneGrid,
                            verbose = TRUE#,
                            #eval_metric="auc", #metric="ROC" does not work with early_stop_rounds
                            #early_stopping_rounds = 10, #if performance not improving for 10 rounds, model iteration stops
                            #maximize = TRUE
    )
  ), "try-error")){
    next
  }
  
  stopCluster(cluster)
  registerDoSEQ()
  
  temp_validation <- dismo::predict(xgb_fit, validation[,colnames(validation) %in% covarsNames],  type = "prob")
  temp_validation <- data.frame("sites"=rownames(validation[,colnames(validation) %in% covarsNames]), "occ"=temp_validation[,"C1"])
  #head(temp_validation)
  
  temp_model_time <- proc.time()[3] - tmp
  temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))
  
  #plot(xgb_fit)
  
  #print(xgb_fit$bestTune)
  
  # calculate variable importance (for later)
  xgb_varImp <- caret::varImp(xgb_fit, scale=T) #scaled between 0 and 100% 
  xgb_varImp <- data.frame("importance" = xgb_varImp$importance, "Predictor"=rownames(xgb_varImp$importance))
  
  tmp <- proc.time()[3]
  # predict to raster (for later)
  temp_prediction <- dismo::predict(xgb_fit, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames],  type = "prob")[,"C1"]
  temp_prediction <- as.numeric(temp_prediction)
  # add names of grid cell (only for those that have no NA in any layer)
  names(temp_prediction) <- rownames(Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),])
  temp_prediction <- as.data.frame(temp_prediction)
  temp_prediction$x <- Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),]$x
  temp_prediction$y <- Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),]$y
  temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
    rename("layer" = temp_prediction)
  
  temp_predict_time <- proc.time()[3] - tmp
  temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))
  
  # NULL for xgb_fit (to save memory)
  xgb_pred_list[[no.runs]] <- list(temp_validation, modelName, temp_model_time, xgb_varImp, temp_prediction, temp_predict_time)
  rm(temp_validation, temp_prediction, temp_model_time, temp_predict_time, xgb_varImp, xgb_fit)
}

# merge all XGB predictions
temp_validation <- do.call(rbind, lapply(xgb_pred_list, "[[", 1))
temp_validation <- temp_validation %>% group_by(sites) %>% summarise(across(everything(), mean))

temp_model_time <- c(sum(as.numeric(sapply(xgb_pred_list, "[[", 3)[1,]), na.rm=T), xgb_pred_list[[1]][[3]][2])
temp_predict_time <- c(sum(as.numeric(sapply(xgb_pred_list, "[[", 6)[1,]), na.rm=T), xgb_pred_list[[1]][[6]][2])

#temp_models <- sapply(xgb_pred_list, "[[", 1)

temp_varImp <- do.call(rbind, lapply(xgb_pred_list, "[[", 4)) %>% 
  group_by(Predictor) %>% summarise_all(mean, na.rm=T)
temp_varImp <- as.data.frame(temp_varImp)
colnames(temp_varImp)[2] <- "xgb"
#print(temp_varImp)

temp_runs <- length(xgb_pred_list)

## extract predicted probabilities (for later)
temp_prediction <- do.call(rbind, lapply(xgb_pred_list, "[[", 5)) %>% 
  group_by(x,y) %>% summarise_all(mean, na.rm=T)

xgb_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp)
save(xgb_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_xgb_", spID, ".RData"))

rm(xgb_list, xgb_pred_list, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp)

})}
stopImplicitCluster()


registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
         .export = c("Env_norm", "Env_norm_df", "form"),
         .packages = c("tidyverse","biomod2", "caret", "randomForest")) %dopar% { try({

#- - - - - - - - - - - - - - - - - - - - -
## cforest ####
#- - - - - - - - - - - - - - - - - - - - -

# cforest will not be used here because of long computational time and 
# comparibly low performance (see Valavi et al. 2021).

# modelName <- "bg.rf"
# 
# # identify and load all relevant data files
# temp.files <- list.files(path = paste0("./results/",Taxon_name), 
#                          pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
# lapply(temp.files, load, .GlobalEnv)

#- - - - - - - - - - - - - - - - - - - - -
## RF and RF-downsampled ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/2

rf_pred_list <- list()
rf_downsample_pred_list <- list()
temp_varImp_rf <- vector(mode="list", length=no.loop.runs)
temp_varImp_rf_downsample <- vector(mode="list", length=no.loop.runs)

for(no.runs in 1:no.loop.runs){
  
  # load training (and testing) data
  temp.files.subset <- list.files(path = paste0("./results/",Taxon_name), 
                                  pattern = paste0(modelName, no.runs, "_[[:graph:]]*", spID), full.name = T)
  lapply(temp.files.subset, load, .GlobalEnv)
  
  # convert the response to factor for producing class relative likelihood
  training$occ <- as.factor(training$occ)
  tmp <- proc.time()[3]
  set.seed(32639)
  
  rf <- randomForest::randomForest(formula = occ ~.,
                                   data = training,
                                   importance = T, # assess importance of predictors
                                   ntree = 500) # the default number of trees = 500
  
  # predict with RF
  temp_validation <- dismo::predict(rf, validation[,colnames(validation) %in% covarsNames], type = "prob") # prob = continuous prediction
  temp_validation <- data.frame("sites"=rownames(validation[,colnames(validation) %in% covarsNames]), "occ"=temp_validation[, "1"])
  #head(temp_validation)
  
  temp_model_time <- proc.time()[3] - tmp
  temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))
  
  tmp <- proc.time()[3]

  # predict to raster (for later)
  temp_prediction <- dismo::predict(rf, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames],  type = "prob")[, "1"]
  temp_prediction <- as.numeric(temp_prediction)
  # add names of grid cell (only for those that have no NA in any layer)
  names(temp_prediction) <- rownames(Env_norm_df)
  temp_prediction <- as.data.frame(temp_prediction)
  temp_prediction$x <- Env_norm_df$x
  temp_prediction$y <- Env_norm_df$y
  temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
    rename("layer" = temp_prediction)
  
  # get running time
  temp_predict_time <- proc.time()[3] - tmp
  temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))
  
  # varImp
  temp_temp_varImp <- as.data.frame(rf$importance)
  # check if it looks right
  #print(temp_temp_varImp)
  # extract varImp as mean decreae in Gini index (alternative: decrease in accuracy)
  temp_temp_varImp <- temp_temp_varImp %>% mutate("Predictor" =rownames(temp_temp_varImp))
  temp_temp_varImp <- temp_temp_varImp[,c("Predictor", "MeanDecreaseGini")]
  colnames(temp_temp_varImp)[2] <- "rf"
  temp_varImp_rf[[no.runs]] <- temp_temp_varImp

  rf_pred_list[[no.runs]] <- list(rf, temp_validation, modelName, temp_model_time, temp_prediction, temp_predict_time)
  rm(temp_temp_varImp)  

  #plot(rf, main = "RF")
  
  ## down-sampling RF
  # for this, background and presence data should be the same number
  prNum <- as.numeric(table(training$occ)["1"]) # number of presences
  bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
  
  # the sample size in each class; the same as presence number
  smpsize <- c("0" = prNum, "1" = prNum)
  tmp <- proc.time()[3]
  set.seed(32639)
  
  rf_downsample <- randomForest::randomForest(formula = occ ~.,
                                              data = training,
                                              ntree = 1000,
                                              sampsize = smpsize,
                                              importance = T,
                                              replace = TRUE)
  #plot(rf_downsample, main = "RF down-sampled")
  
  # predict with RF down-sampled
  temp_validation <- predict(rf_downsample, validation[,colnames(validation) %in% covarsNames], type = "prob") # prob = continuous prediction
  temp_validation <- data.frame("sites"=rownames(validation[,colnames(validation) %in% covarsNames]), "occ"=temp_validation[, "1"])
  #head(temp_validation)

  temp_model_time <- proc.time()[3] - tmp
  temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

  tmp <- proc.time()[3]

  # predict to raster (for later)
  temp_prediction <- dismo::predict(rf_downsample, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames],  type = "prob")[,"1"]
  temp_prediction <- as.numeric(temp_prediction)
  # add names of grid cell (only for those that have no NA in any layer)
  names(temp_prediction) <- rownames(Env_norm_df)
  temp_prediction <- as.data.frame(temp_prediction)
  temp_prediction$x <- Env_norm_df$x
  temp_prediction$y <- Env_norm_df$y
  temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
    rename("layer" = temp_prediction)
   
  # get running time
  temp_predict_time <- proc.time()[3] - tmp
  temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

  # varImp
  temp_temp_varImp <- as.data.frame(rf_downsample$importance)
  # check if it looks right
  #print(temp_temp_varImp)
  # extract varImp as mean decreae in Gini index (alternative: decrease in accuracy)
  temp_temp_varImp <- temp_temp_varImp %>% mutate("Predictor" =rownames(temp_temp_varImp))
  temp_temp_varImp <- temp_temp_varImp[,c("Predictor", "MeanDecreaseGini")]
  colnames(temp_temp_varImp)[2] <- "rf_downsample"
  temp_varImp_rf_downsample[[no.runs]] <- temp_temp_varImp

  rf_downsample_pred_list[[no.runs]] <- list(rf_downsample, temp_validation, modelName, temp_model_time, temp_prediction, temp_predict_time)
  rm(rf, rf_downsample, temp_validation, temp_prediction, temp_model_time, temp_predict_time, temp_temp_varImp)
}

# merge all RF predictions
temp_validation <- do.call(rbind, lapply(rf_pred_list, "[[", 2))
temp_validation <- temp_validation %>% group_by(sites) %>% summarise(across(everything(), mean))

temp_model_time <- c(sum(as.numeric(sapply(rf_pred_list, "[[", 4)[1,]), na.rm=T), rf_pred_list[[1]][[4]][2])
#temp_model_time <- sum(sapply(rf_pred_list, "[[", 4)[1,], na.rm=T)
temp_predict_time <- c(sum(as.numeric(sapply(rf_pred_list, "[[", 6)[1,]), na.rm=T), rf_pred_list[[1]][[6]][2])
#temp_predict_time <- sum(sapply(rf_pred_list, "[[", 6)[1,], na.rm=T)

#temp_models <- sapply(rf_pred_list, "[[", 1)

## extract predicted probabilities (for later)
temp_prediction <- do.call(rbind, lapply(rf_pred_list, "[[", 5)) %>% 
  group_by(x,y) %>% summarise_all(mean, na.rm=T)

temp_runs <- length(rf_pred_list)

temp_varImp_rf <- do.call(rbind, temp_varImp_rf) %>% 
      group_by(Predictor) %>% summarise_all(mean, na.rm=T) %>%
      as.data.frame()

rf_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp_rf)
save(rf_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_rf_", spID, ".RData"))

rm(rf_pred_list, rf_list, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp_rf)

# merge all RF_downsample predictions
temp_validation <- do.call(rbind, lapply(rf_downsample_pred_list, "[[", 2))
temp_validation <- temp_validation %>% group_by(sites) %>% summarise(across(everything(), mean))

temp_model_time <- c(sum(as.numeric(sapply(rf_downsample_pred_list, "[[", 4)[1,]), na.rm=T), rf_downsample_pred_list[[1]][[4]][2])
temp_predict_time <- c(sum(as.numeric(sapply(rf_downsample_pred_list, "[[", 6)[1,]), na.rm=T), rf_downsample_pred_list[[1]][[6]][2])

temp.models <- sapply(rf_downsample_pred_list, "[[", 1)

## extract predicted probabilities (for later)
temp_prediction <- do.call(rbind, lapply(rf_downsample_pred_list, "[[", 5)) %>% 
  group_by(x,y) %>% summarise_all(mean, na.rm=T)

temp_runs <- length(rf_downsample_pred_list)

temp_varImp_rf_downsample <- do.call(rbind, temp_varImp_rf_downsample) %>% 
      group_by(Predictor) %>% summarise_all(mean, na.rm=T) %>%
      as.data.frame()

rf_downsample_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp_rf_downsample)
save(rf_downsample_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_rf_downsample_", spID, ".RData"))

rm(rf_downsample_pred_list, modelName, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp_rf_downsample)

#- - - - - - - - - - - - - - - - - - - - - 
## Model RF2 for ensemble modelling using consistent background data (bg.glm)

modelName <- "bg.glm"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)

lapply(temp.files, load, .GlobalEnv)

# convert the response to factor for producing class relative likelihood
training$occ <- as.factor(training$occ)
tmp <- proc.time()[3]
set.seed(32639)
rf2 <- randomForest(formula = occ ~.,
                    data = training,
                    ntree = 500) # the default number of trees

# predict with RF
temp_validation <- dismo::predict(rf2, validation[,colnames(validation) %in% covarsNames], type = "prob")[, "1"] # prob = continuous prediction
#head(temp_validation)
temp_validation <- as.numeric(temp_validation)
names(temp_validation) <- rownames(validation[,colnames(validation) %in% covarsNames]) #add site names

temp_model_time <- proc.time()[3] - tmp
temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

tmp <- proc.time()[3] 
# predict to raster (for later)
temp_prediction <- dismo::predict(rf2, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames],  type = "prob")[, "1"]
temp_prediction <- as.numeric(temp_prediction)
# add names of grid cell (only for those that have no NA in any layer)
names(temp_prediction) <- rownames(Env_norm_df)
temp_prediction <- as.data.frame(temp_prediction)
temp_prediction$x <- Env_norm_df$x
temp_prediction$y <- Env_norm_df$y
temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
  rename("layer" = temp_prediction)

# get running time
temp_predict_time <- proc.time()[3] - tmp
temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

temp_runs <- 1

# variable importance
temp_varImp <- as.data.frame(rf2$importance)
temp_varImp$Predictor <- rownames(temp_varImp)
colnames(temp_varImp)[1] <- "rf2"

rf2_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp)
save(rf2_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_rf2_", spID, ".RData"))

rm(rf2, rf2_list, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp)

})}
stopImplicitCluster()


registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
         .export = c("Env_norm", "Env_norm_df", "form"),
         .packages = c("tidyverse","biomod2", "caret", "e1071")) %dopar% { try({

#- - - - - - - - - - - - - - - - - - - - -
## SVM ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/2

svm_pred_list <- list()
temp_varImp <- vector(mode="list", length=no.loop.runs)

for(no.runs in 1:no.loop.runs){
  
  temp.files.subset <- list.files(path = paste0("./results/",Taxon_name), 
                                  pattern = paste0(modelName, no.runs, "_[[:graph:]]*", spID), full.name = T)
  lapply(temp.files.subset, load, .GlobalEnv)
  
  # change the response to factor
  training$occ <- as.factor(training$occ)
  
  # calculate weights the class weight
  prNum <- as.numeric(table(training$occ)["1"]) # number of presences
  bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
  cwt <- c("0" = prNum / bgNum, "1" = 1)
  
  tmp <- proc.time()[3]
  set.seed(32639)
  
  svm_e <- e1071::svm(formula = occ ~ .,
                      data = training,
                      kernel = "radial", 
                      class.weights = cwt,
                      probability = TRUE) # allow for probability predictions
  
  # predicting on test data
  temp_validation <- e1071:::predict.svm(svm_e, validation[,colnames(validation) %in% covarsNames], probability=TRUE)
  temp_validation <- attr(temp_validation, "probabilities")[,"1"] 
  temp_validation <- data.frame("sites"=rownames(validation[,colnames(validation) %in% covarsNames]), "occ"=temp_validation)
  #head(temp_validation)
    
  temp_model_time <- proc.time()[3] - tmp
  temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))
  
  tmp <- proc.time()[3]
  # predict
  temp_prediction <- e1071:::predict.svm(svm_e, Env_norm_df[,colnames(Env_norm_df) %in% covarsNames], probability=TRUE)
  temp_prediction <- attr(temp_prediction, "probabilities")[,"1"]
  temp_prediction <- as.numeric(temp_prediction)
  # add names of grid cell (only for those that have no NA in any layer)
  names(temp_prediction) <- rownames(Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),])
  temp_prediction <- as.data.frame(temp_prediction)
  temp_prediction$x <- Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),]$x
  temp_prediction$y <- Env_norm_df[complete.cases(Env_norm_df[,colnames(Env_norm_df) %in% covarsNames]),]$y
  temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
    rename("layer" = temp_prediction)

  # get running time
  temp_predict_time <- proc.time()[3] - tmp
  temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

  # varImp
  # calculate variable importance according to 
  # https://stackoverflow.com/questions/34781495/how-to-find-important-factors-in-support-vector-machine
  w <- t(svm_e$coefs) %*% svm_e$SV # weight vectors
  w <- apply(w, 2, function(v){sqrt(sum(v^2))})  # weight
  w <- sort(w, decreasing = T)  
  # structure variable importance
  temp_temp_varImp <- w %>% as.data.frame() %>% 
        mutate(Predictor = names(w))
  colnames(temp_temp_varImp)[1] <- "svm" 
  # check if it looks right
  #print(temp_temp_varImp)
  temp_varImp[[no.runs]] <- temp_temp_varImp

  svm_pred_list[[no.runs]] <- list(svm_e, temp_validation, modelName, temp_model_time, temp_prediction, temp_predict_time)
  rm(svm_e, temp_validation, temp_model_time, temp_prediction, temp_predict_time, temp_temp_varImp)
}

# merge all SVM predictions
temp_validation <- do.call(rbind, lapply(svm_pred_list, "[[", 2))
temp_validation <- temp_validation %>% group_by(sites) %>% summarise(across(everything(), mean))

temp_model_time <- c(sum(as.numeric(sapply(svm_pred_list, "[[", 4)[1,]), na.rm=T), svm_pred_list[[1]][[4]][2])
temp_predict_time <- c(sum(as.numeric(sapply(svm_pred_list, "[[", 6)[1,]), na.rm=T), svm_pred_list[[1]][[6]][2])

#temp_models <- sapply(svm_pred_list, "[[", 1)

## extract predicted probabilities (for later)
temp_prediction <- do.call(rbind, lapply(svm_pred_list, "[[", 5)) %>% 
  group_by(x,y) %>% summarise_all(mean, na.rm=T)

temp_runs <- length(svm_pred_list)

# varImp
temp_varImp <- do.call(rbind, temp_varImp) %>% 
      group_by(Predictor) %>% summarise_all(mean, na.rm=T) %>%
      as.data.frame()

svm_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction, varImp=temp_varImp)
save(svm_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_svm_", spID, ".RData"))

rm(svm_list, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp)

})}
stopImplicitCluster()


registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
         .export = c("Env_norm", "Env_norm_df", "form"),
         .packages = c("tidyverse","biomod2", "caret", "earth")) %dopar% { try({

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
#   MAXENT.Phillips = list( path_to_maxent.jar = paste0(here::here(), "/results/", Taxon_name, "/maxent_files"), # change it to maxent directory
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
  #MAXENT.Phillips = list( path_to_maxent.jar =paste0(here::here(), "/results" )), # change it to maxent directory
  )

# models to predict with
mymodels <- c("GLM","GBM","GAM","CTA","ANN","FDA","MARS","RF","MAXENT.Phillips")

# model fitting
tmp <- proc.time()[3]
setwd(paste0(here::here(), "/results/", Taxon_name))

set.seed(32639)
myBiomodModelOut <- biomod2::BIOMOD_Modeling(myBiomodData,
                                             models = mymodels,
                                             models.options = myBiomodOption,
                                             NbRunEval = 10,   # 3-fold crossvalidation evaluation run
                                             DataSplit = 80, # use subset of the data for training
                                             models.eval.meth = c("ROC"),
                                             SaveObj = FALSE, #save output on hard drive?
                                             rescal.all.models = FALSE, #scale all predictions with binomial GLM?
                                             do.full.models = FALSE, # do evaluation & calibration with whole dataset
                                             modeling.id = paste(myRespName,"_Modeling", sep = ""))

# ensemble modeling using mean probability
myBiomodEM <- biomod2::BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                               chosen.models = "all",  # all algorithms
                                               em.by = "all",    #evaluated over evaluation data if given (it is not, see Prepare_input.R)
                                               # note: evaluation not that important as we will calculate measures on independent data
                                               eval.metric = c("ROC"), # 'all' would takes same as above in BIOMOD_Modelling
                                               eval.metric.quality.threshold = NULL, # since some species's auc are naturally low
                                               prob.mean = TRUE, #estimate mean probabilities across predictions
                                               prob.cv = FALSE,   #estimate coefficient of variation across predictions
                                               prob.ci = FALSE,  #estimate confidence interval around the prob.mean
                                               prob.median = FALSE, #estimate the median of probabilities
                                               committee.averaging = FALSE, #estimate committee averaging across predictions
                                               prob.mean.weight = TRUE, #estimate weighted sum of predictions
                                               prob.mean.weight.decay = "proportional", #the better a model (performance), the higher weight
                                               VarImport = 1)    #number of permutations to estimate variable importance
temp_model_time <- proc.time()[3] - tmp
temp_model_time <- c(round(as.numeric(temp_model_time), 3), units(temp_model_time))

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
temp_predict_time <- c(round(as.numeric(temp_predict_time), 3), units(temp_predict_time))

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
save(biomod_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_biomod_", spID, ".RData"))

rm(biomod_list, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp, myBiomodEM, myBiomodEnProj, myBiomodModelOut, myBiomodProj, myBiomodModelEval, myEnProjDF)

setwd(here::here())  

})}
stopImplicitCluster()


registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
         .export = c("Env_norm", "Env_norm_df", "form"),
         .packages = c("tidyverse","biomod2", "caret", "earth", "doParallel")) %dopar% { try({

#- - - - - - - - - - - - - - - - - - - - -
## simple Ensemble model ####
#- - - - - - - - - - - - - - - - - - - - -
# note: this model should be run after all the component 
# models (MaxEnt, Lasso, GAM, BRT, RF down-sampled)

modelName <- "bg.glm"

temp_validation <- c()
temp_prediction <- data.frame("layer"=NA, "x"=NA, "y"=NA)[0,]
temp_model_time <- c()
temp_predict_time <- c()
temp_models <- c()

for(m in c("gm", "lasso", "brt2", "rf2", "maxent")){ try({
  # load and scale validation data
  temp_list <- get(load(file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_",m, "_", spID, ".RData")))
  temp_vali <- scales::rescale(temp_list[["validation"]], to = c(0,1))

  temp_validation <- cbind(temp_validation, temp_vali)
  
  # sum up computational time
  temp_model_time <- paste0(temp_model_time, temp_list[["time_model"]], " + ") 
  temp_predict_time <- paste0(temp_predict_time, temp_list[["time_model"]], " + ") 

  # average predictions
  temp_pred <- temp_list[["prediction"]]
  temp_pred[,m] <- temp_pred$layer
  temp_prediction <- temp_prediction %>% full_join(temp_pred, by=c("x", "y"))
    
  temp_models <- c(temp_models, m)
  
}, silent=T)}

# average the predictions
temp_validation <- rowMeans(temp_validation, na.rm=T) 

temp_runs <- 1

temp_prediction <- temp_prediction %>% mutate("layer"=rowMeans(temp_prediction %>% dplyr::select(-x, -y), na.rm=T)) %>%
  dplyr::select(x, y, layer)

ensm_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, validation=temp_validation, prediction=temp_prediction)
save(ensm_list, file=paste0(here::here(),"/results/", Taxon_name, "/temp_files/SDM_ensm_", spID, ".RData"))

rm(temp_list, ensm_list, temp_prediction, temp_pred, temp_vali, temp_model_time, temp_predict_time, temp_runs, temp_validation)

#- - - - - - - - - - - - - - - - - - - - -
## Check if everything went well ####
#- - - - - - - - - - - - - - - - - - - - -
})}
stopImplicitCluster()


for(spID in unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID)){ try({

check_files <- list.files(paste0(here::here(),"/results/", Taxon_name, "/temp_files"))
check_files <- check_files[stringr::str_detect(check_files, spID)]

print(check_files)
if(length(check_files)==17){ print("Everything looks well. Please continue :)")
}else{ print("Some algorithms weren't saved; please check which ones are missing.")}


# #- - - - - - - - - - - - - - - - - - - - -
# ## Saving ####
# #- - - - - - - - - - - - - - - - - - - - -
# 
# # save all models to calculate model performance later
# SDMs <- list(gm_pred, lm1_pred, lm_subset_pred, lasso_pred, ridge_pred, 
#              mars_pred, maxmod_pred, maxnet_pred, brt_pred, brt2_pred, xgb_pred, svm_pred,
#              rf_pred, rf2_pred, rf_downsample_pred, biomod_pred, ensm_pred)
# names(SDMs) <- c("gm_pred", "lm1_pred", "lm_subset_pred", "lasso_pred", "ridge_pred", 
#                  "mars_pred", "maxmod_pred", "maxnet_pred", "brt_pred", "brt2_pred", "xgb_pred", "svm_pred",
#                  "rf_pred", "rf2_pred", "rf_downsample_pred", "biomod_pred", "ensm_pred")
# #head(SDMs)
# 
# save(SDMs, file=paste0(here::here(), "/sdm/SDM_Models_", spID, ".RData")) 

})}

stopImplicitCluster()

