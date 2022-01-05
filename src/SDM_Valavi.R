#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(mgcv) # for GAM
library(gam)  # for GLM (!)
library(remotes) #to download package from github

## for regularized regressions
# installing the package from github
remotes::install_github("rvalavi/myspatial")
library(glmnet)
library(myspatial)

library(caret) # for MARS and BRT
library(earth) # for MARS
library(doParallel) # for MARS and XGBoost

library(dismo) # for MaxEnt and BRT
library(maxnet) # for MaxNet

library(xgboost) # for XGBoost
library(randomForest) # for RF
library(e1071) # for SVM

library(precrec) # for performance measure
# devtools::install_github("meeliskull/prg/R_package/prg")
library(prg)

# note: we will load the datasets before each individual model

#- - - - - - - - - - - - - - - - - - - - -
## GAM ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.glm"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
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

# define formula
form <- paste0("occ ~ ", paste0(paste0("s(", covarsNames, ")"), collapse=" + "))

# general settings
tmp <- Sys.time()
set.seed(32639)

# run model
gm <- mgcv::gam(formula = as.formula(form),
                data = training,
                family = binomial(link = "logit"),
                weights = wt,
                method = "REML")

# get running time
Sys.time() - tmp

# check the appropriateness of Ks
# the model parameter k should not be higher than the number of unique values 
# in each covariate
par(mfrow=c(2,2))
gam.check(gm)
# Note: k-index The further below 1 this is, the more likely it is that there 
# is missed pattern left in the residuals;
# and if edf close to k', basic dimension k' has been set too low

par(mfrow=c(1,1))
#plot(gm, pages = 1, rug = TRUE, shade = TRUE)
# Note: plots are useless for binary data

# predict
gm_pred <- predict(gm, testing_env)
head(gm_pred)

gm_pred <- list(gm_pred, modelName)

#- - - - - - - - - - - - - - - - - - - - -
## GLM ####
#- - - - - - - - - - - - - - - - - - - - -

# calculating the weights
# the order of weights should be the same as presences and backgrounds in the 
# training data
prNum <- as.numeric(table(training$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)

# the base glm model with linear terms
lm1 <- glm(occ ~., data = training, weights = wt, 
           family = binomial(link = "logit"))
summary(lm1)

lm1_pred <- predict(lm1, testing_env)
head(lm1_pred)

lm1_pred <- list(lm1_pred, modelName)

# model scope for subset selection
mod_scope <- gam.scope(frame = training, form=F,
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

tmp <- Sys.time()
set.seed(32639)

lm_subset <- step.Gam(object = lm1,
                      scope = mod_scope,
                      direction = "both",
                      data = training, # this is optional
                      trace = FALSE)
Sys.time() - tmp

summary(lm_subset)

lm_subset_pred <- predict(lm_subset, testing_env)
head(lm_subset_pred)

lm_subset_pred <- list(lm_subset_pred, modelName)

#- - - - - - - - - - - - - - - - - - - - -
## Regularized regressions: lasso and ridge regression ####
##- - - - - - - - - - - - - - - - - - - - -

# generating the quadratic terms for all continuous variables
# function to creat quadratic terms for lasso and ridge
quad_obj <- make_quadratic(training, cols = covarsNames)

# now we can predict this quadratic object on the training and testing data
# this make two columns for each covariates used in the transformation
training_quad <- predict(quad_obj, newdata = training)
testing_quad <- predict(quad_obj, newdata = testing_env)

# convert the data.frames to sparse matrices
# select all quadratic (and non-quadratic) columns, except the y (occ)
new_vars <- names(training_quad)[names(training_quad) != "occ"]
training_sparse <- sparse.model.matrix(~. -1, training_quad[, new_vars])
testing_sparse <- sparse.model.matrix( ~. -1, testing_quad[, new_vars])

# calculating the case weights
prNum <- as.numeric(table(training_quad$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training_quad$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)

## fitting lasso
# Note: The alpha parameter in this model ranges from 0 to 1, where selecting an 
# alpha of 0 leads to ridge regression and 1 to lasso and anything in between 
# is a combination of both called elastic-net.

# # plot the regularization path and shrinkage in the coefficients
# plot(lasso, xvar = "lambda", label = TRUE)
# # As you can see by changing (log) Lambda the coefficients shrink (the y-axes) 
# # and the number of covariates included in the model, decreases (x-axes, top) 
# # as the coefficients can be set to zeros in the lasso

# The lambda parameter controls regularization – it is the penalty applied 
# to the model’s coefficients. To select the best lambda, internal cross-
# validation was used (in cv.glmnet function).
tmp <- Sys.time()
set.seed(32639)
lasso_cv <- cv.glmnet(x = training_sparse,
                      y = training_quad$occ,
                      family = "binomial",
                      alpha = 1, # fitting lasso
                      weights = wt,
                      nfolds = 10) # number of folds for cross-validation
Sys.time() - tmp

# the cross-validation result
plot(lasso_cv)
# Note: dashed lines show options for best Lambda parameter
# left: Lambda with minimum deviance (lambda.min)
# right: best Lambda within one standard deviation of left-dashed line (l...1se)

print("One of the two Lambda (dashed lines) needs to be selected for prediction.")

# fitting ridge resgression (alpha=0) while identify the right lambda
tmp <- Sys.time()
set.seed(32639)
ridge_cv <- cv.glmnet(x = training_sparse,
                      y = training_quad$occ,
                      family = "binomial",
                      alpha = 0, # fitting lasso
                      weights = wt,
                      nfolds = 10) # number of folds for cross-validation
Sys.time() - tmp

plot(ridge_cv)

# predicting with lasso model
lasso_pred <- predict(lasso_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
head(lasso_pred)

lasso_pred <- list(lasso_pred, modelName)

# predict ridge
ridge_pred <- predict(ridge_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
head(ridge_pred)

ridge_pred <- list(ridge_pred, modelName)

#- - - - - - - - - - - - - - - - - - - - -
## MARS ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.mars"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# change the response to factor variables with non-numeric levels
training$occ <- as.factor(training$occ)
levels(training$occ) <- c("C0", "C1")

mytuneGrid <- expand.grid(nprune = 2:20,
                          degree = 2) # no interaction = 1
mytrControl <- trainControl(method = "cv",
                            number = 5, # 5-fold cross-validation
                            classProbs = TRUE,
                            summaryFunction = twoClassSummary,
                            allowParallel = TRUE)
tmp <- Sys.time()
cluster <- makeCluster(6) # you can use all cores of your machine instead e.g. 8
registerDoParallel(cluster)
set.seed(32639)
mars_fit <- train(form = occ ~ .,
                  data = training,
                  method = "earth",
                  metric = "ROC",
                  trControl = mytrControl,
                  tuneGrid = mytuneGrid,
                  thresh = 0.00001)
stopCluster(cluster)
registerDoSEQ()
Sys.time() - tmp
plot(mars_fit)

mars_pred <- predict(mars_fit, testing_env)
# transform occurrence back into numeric
mars_pred <- as.character(mars_pred)
mars_pred[mars_pred=="C0"] <- 0
mars_pred[mars_pred=="C1"] <- 1
mars_pred <- as.numeric(mars_pred)
names(mars_pred) <- rownames(testing_env) #add site names
head(mars_pred)

mars_pred <- list(mars_pred, modelName)

# transform occurrence column back to numeric
training$occ <- as.numeric(training$occ)

#- - - - - - - - - - - - - - - - - - - - -
## MaxEnt and MaxNet ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.glm"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# "We used five different regularization multipliers (0.5, 1, 2, 3 and 4)
# in combination with different features (L, LQ, H, LQH, LQHP) to find the 
# best parameters that maximizes the average AU CROC in CV."

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
  covars <- names(data)[which(names(data) != y)]
  # regularization multipliers
  ms <- c(0.5, 1, 2, 3, 4)
  grid <- expand.grid(
    regmult = paste0("betamultiplier=", ms),
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
        maxmod <- dismo::maxent(x = data[trainSet, covars],
                                p = data$occ[trainSet],
                                removeDuplicates = FALSE,
                                path = filepath,
                                args = as.character(unlist(grid[n, ])))
      ), "try-error")){
        next
      }
      modpred <- predict(maxmod, data[testSet, covars], args = "outputformat=cloglog")
      pred_df <- data.frame(score = modpred, label = data$occ[testSet])
      full_pred <- rbind(full_pred, pred_df)
    }
    AUCs[n] <- precrec::auc(precrec::evalmod(scores = full_pred$score,
                                             labels = full_pred$label))[1,4]
  }
  best_param <- as.character(unlist(grid[which.max(AUCs), ]))
  return(best_param)
}

# now use the function to tune MaxEnt
# number of folds
nfolds <- ifelse(sum(as.numeric(training$occ)) < 10, 2, 5)
tmp <- Sys.time()
set.seed(32639)
# tune maxent parameters
param_optim <- maxent_param(data = training,
                            k = nfolds,
                            filepath = paste0("results/", Taxon_name, "/maxent_files"))
# fit a maxent model with the tuned parameters
maxmod <- dismo::maxent(x = training[, covars],
                        p = training$occ,
                        removeDuplicates = FALSE,
                        path = paste0("results/", Taxon_name, "/maxent_files"),
                        rgs = param_optim)
Sys.time() - tmp

maxmod_pred <- predict(maxmod, testing_env)
names(maxmod_pred) <- rownames(testing_env) #add site names
head(maxmod_pred)

maxmod_pred <- list(maxmod_pred, modelName)

## MaxNet
presences <- training$occ # presence (1s) and background (0s) points
covariates <- training[, 2:ncol(training)] # predictor covariates
tmp <- Sys.time()
set.seed(32639)
mxnet <- maxnet::maxnet(p = presences,
                        data = covariates,
                        regmult = 1, # regularization multiplier
                        maxnet::maxnet.formula(presences, covariates, classes = "default"))
Sys.time() - tmp

# predicting with MaxNet
maxnet_pred <- predict(mxnet, testing_env, type = c("cloglog"))[, 1]
head(maxnet_pred)

maxnet_pred <- list(maxnet_pred, modelName)

#- - - - - - - - - - - - - - - - - - - - -
## BRT (GBM) ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# calculating the case weights
prNum <- as.numeric(table(training$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)

tmp <- Sys.time()
set.seed(32639)
brt <- gbm.step(data = training,
                gbm.x = 2:ncol(training), # column indices for covariates
                gbm.y = 1, # column index for response
                family = "bernoulli",
                tree.complexity = ifelse(prNum < 50, 1, 5),
                learning.rate = 0.001,
                bag.fraction = 0.75,
                max.trees = 10000,
                n.trees = 50,
                n.folds = 5, # 5-fold cross-validation
                site.weights = wt,
                silent = TRUE) # avoid printing the cv results

Sys.time() - tmp
# Note: model tuning with ~50,000 background points may take ~1h.

#the optimal number of trees (intersect of green and red line in plot)
brt$gbm.call$best.trees

# predicting with the best trees
brt_pred <- predict(brt, testing_env, n.trees = brt$gbm.call$best.trees, type = "response")
head(brt_pred)

brt_pred <- list(brt_pred, modelName)

#- - - - - - - - - - - - - - - - - - - - -
## XGBoost ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# change the response to factor variables with non-numeric levels
training$occ <- as.factor(training$occ)
levels(training$occ) <- c("C0", "C1") # caret does not accept 0 and 1 as factor levels

# train control for cross-validation for model tuning
mytrControl <- trainControl(method = "cv",
                            number = 5, # 5 fold cross-validation
                            summaryFunction = twoClassSummary,
                            classProbs = TRUE,
                            allowParallel = TRUE)
# setting the range of parameters for grid search tuning
mytuneGrid <- expand.grid(
  nrounds = seq(from = 500, to = 15000, by = 500),
  eta = 0.001,
  max_depth = 5,
  subsample = 0.75,
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1
)
tmp <- Sys.time()
cluster <- makeCluster(6) # you can use all cores of your machine instead e.g. 8
registerDoParallel(cluster)
set.seed(32639)
xgb_fit <- train(form = occ ~ .,
                 data = training,
                 method = "xgbTree",
                 metric = "ROC",
                 trControl = mytrControl,
                 tuneGrid = mytuneGrid,
                 verbose = TRUE)

stopCluster(cluster)
registerDoSEQ()
Sys.time() - tmp

plot(xgb_fit)

print(xgb_fit$bestTune)

xgb_pred <- predict(xgb_fit, testing_env)
head(xgb_pred)

xgb_pred <- list(xgb_pred, modelName)

# transform occurrence column back to numeric
training$occ <- as.numeric(training$occ)

#- - - - - - - - - - - - - - - - - - - - -
## cforest ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# cforest will not be used here because of long computational time and 
# comparibly low performance (see Valavi et al. 2021).


#- - - - - - - - - - - - - - - - - - - - -
## RF and RF-downsampled ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# convert the response to factor for producing class relative likelihood
training$occ <- as.factor(training$occ)
tmp <- Sys.time()
set.seed(32639)
rf <- randomForest(formula = occ ~.,
                   data = training,
                   ntree = 500) # the default number of trees
Sys.time() - tmp

# predict with RF
rf_pred <- predict(rf, testing_env, type = "prob")[, "1"] # prob = continuous prediction
head(rf_pred)

rf_pred <- list(rf_pred, modelName)

plot(rf, main = "RF")

## down-sampling RF
# for this, background and presence data should be the same number
prNum <- as.numeric(table(training$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds

# the sample size in each class; the same as presence number
smpsize <- c("0" = prNum, "1" = prNum)
tmp <- Sys.time()
set.seed(32639)
rf_downsample <- randomForest(formula = occ ~.,
                              data = training,
                              ntree = 1000,
                              sampsize = smpsize,
                              replace = TRUE)
Sys.time() - tmp

plot(rf_downsample, main = "RF down-sampled")

# predict with RF down-sampled
rf_downsample_pred <- predict(rf_downsample, testing_env, type = "prob")[, "1"] # prob = continuous prediction
head(rf_downsample_pred)

rf_downsample_pred <- list(rf_downsample_pred, modelName)

# transform occurrence column back to numeric
training$occ <- as.numeric(training$occ)

#- - - - - - - - - - - - - - - - - - - - -
## SVM ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# change the response to factor
training$occ <- as.factor(training$occ)

# calculate weights the class weight
prNum <- as.numeric(table(training$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
cwt <- c("0" = prNum / bgNum, "1" = 1)
tmp <- Sys.time()
set.seed(32639)
svm_e <- e1071::svm(formula = occ ~ .,
                    data = training,
                    kernel = "radial",
                    lass.weights = cwt,
                    probability = TRUE)
Sys.time() - tmp

# predicting on test data
svm_pred <- predict(svm_e, testing_env, probability = TRUE)
svm_prob <- attr(svm_pred, "probabilities")[,"1"]

# see the first few predictions
head(svm_prob)

svm_pred <- list(svm_prob, modelName)
 
# transform occurrence column back to numeric
training$occ <- as.numeric(training$occ)
                    
#- - - - - - - - - - - - - - - - - - - - -
## biomod ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.glm"

# re-loading the species data (we need the x & y column)
# load background data (pseudo-absences) for each modeling approach
load(file=paste0(here::here(), "/results/", Taxon_name, "/PA_Env_", Taxon_name, "_", spID, ".RData"))
data <- bg.list[[modelName]] %>% rename("occ"=1)

pr <- data[!is.na(data$occ),]
bg <- data[is.na(data$occ),]

# define parameters
training <- rbind(pr, bg)
myRespName <- "occ"
myResp <- as.numeric(training[, myRespName])
myResp[which(myResp == 0)] <- NA
myExpl <- data.frame(training[, covars])
myRespXY <- training[, c("x", "y")]

# create biomod data format
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.name = myRespName,
                                     resp.xy = myRespXY,
                                     PA.nb.absences = 50000,
                                     PA.strategy = 'random',
                                     na.rm = TRUE)

# using the default options
# you can change the mentioned parameters by changes this
myBiomodOption <- BIOMOD_ModelingOptions()

# models to predict with
mymodels <- c("GLM","GBM","GAM","CTA","ANN","FDA","MARS","RF","MAXENT.Phillips")

# model fitting
tmp <- Sys.time()
set.seed(32639)
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = mymodels,
                                    models.options = myBiomodOption,
                                    NbRunEval = 1,
                                    DataSplit = 100, # use all the data for training
                                    models.eval.meth = c("ROC"),
                                    SaveObj = TRUE,
                                    rescal.all.models = FALSE,
                                    do.full.models = TRUE,
                                    modeling.id = paste(myRespName,"NCEAS_Modeling", sep = ""))

# ensemble modeling using mean probability
myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                      chosen.models = 'all',
                                      em.by = 'all',
                                      eval.metric = c("ROC"),
                                      eval.metric.quality.threshold = NULL, # since some species's auc are naturally low
                                      prob.mean = TRUE,
                                      prob.cv = FALSE,
                                      prob.ci = FALSE,
                                      prob.median = FALSE,
                                      committee.averaging = FALSE,
                                      prob.mean.weight = FALSE)
Sys.time() - tmp

# project single models
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                  new.env = as.data.frame(testing_env[, covars]),
                                  proj.name = "nceas_modeling",
                                  selected.models = "all",
                                  binary.meth = "ROC",
                                  compress = TRUE,
                                  clamping.mask = TRUE)

# project ensemble of all models
myBiomodEnProj <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                             EM.output = myBiomodEM,
                                             selected.models = "all")

# extracting the values for ensemble prediction
myEnProjDF <- as.data.frame(get_predictions(myBiomodEnProj))

# see the first few predictions
# the prediction scale of biomod is between 0 and 1000
head(myEnProjDF)

biomod_pred <- myEnProjDF[,1]
names(biomod_pred) <- rownames(testing_env)
biomod_pred <- list(biomod_pred, modelName)
   
                                          
#- - - - - - - - - - - - - - - - - - - - -
## Saving ####
#- - - - - - - - - - - - - - - - - - - - -

# save all models to calculate model performance later
SDMs <- list(gm_pred, lm1_pred, lm_subset_pred, lasso_pred, ridge_pred, 
     mars_pred, maxmod_pred, maxnet_pred, brt_pred, xgb_pred, svm_pred,
     rf_pred, rf_downsample_pred, biomod_pred)
names(SDMs) <- c("gm_pred", "lm1_pred", "lm_subset_pred", "lasso_pred", "ridge_pred", 
                "mars_pred", "maxmod_pred", "maxnet_pred", "brt_pred", "xgb_pred", "svm_pred",
                "rf_pred", "rf_downsample_pred", "biomod_pred")
head(SDMs)

save(SDMs, file=paste0(here::here(), "/results/", Taxon_name, "/Predicted_SDMs_", spID, ".RData"))




