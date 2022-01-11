#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

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
temp.time <- as.numeric(Sys.time() - tmp)

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
#head(gm_pred)

gm_pred <- as.numeric(gm_pred)
names(gm_pred) <- rownames(testing_env) #add site names


gm_pred <- list(gm, gm_pred, modelName, temp.time)

#- - - - - - - - - - - - - - - - - - - - -
## GLM ####
#- - - - - - - - - - - - - - - - - - - - -

# calculating the weights
# the order of weights should be the same as presences and backgrounds in the 
# training data
prNum <- as.numeric(table(training$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)

tmp <- Sys.time()

# the base glm model with linear terms
lm1 <- glm(occ ~., data = training, weights = wt, 
           family = binomial(link = "logit"))
summary(lm1)

temp.time <- as.numeric(Sys.time() - tmp)

lm1_pred <- predict(lm1, testing_env)
#head(lm1_pred)

lm1_pred <- as.numeric(lm1_pred)
names(lm1_pred) <- rownames(testing_env) #add site names

lm1_pred <- list(lm1, lm1_pred, modelName, temp.time)

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
temp.time <- as.numeric(Sys.time() - tmp)

summary(lm_subset)

lm_subset_pred <- predict(lm_subset, testing_env)
#head(lm_subset_pred)

lm_subset_pred <- as.numeric(lm_subset_pred)
names(lm_subset_pred) <- rownames(testing_env) #add site names


lm_subset_pred <- list(lm_subset, lm_subset_pred, modelName, temp.time)

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
temp.time <- as.numeric(Sys.time() - tmp)

# the cross-validation result
#plot(lasso_cv)
# Note: dashed lines show options for best Lambda parameter
# left: Lambda with minimum deviance (lambda.min)
# right: best Lambda within one standard deviation of left-dashed line (l...1se)

print("One of the two Lambda (dashed lines) needs to be selected for prediction.")

# predicting with lasso model
lasso_pred <- predict(lasso_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
#head(lasso_pred)

lasso_pred <- as.numeric(lasso_pred)
names(lasso_pred) <- rownames(testing_env) #add site names


lasso_pred <- list(lasso_cv, lasso_pred, modelName, temp.time)


## fitting ridge resgression (alpha=0) while identify the right lambda
tmp <- Sys.time()
set.seed(32639)
ridge_cv <- cv.glmnet(x = training_sparse,
                      y = training_quad$occ,
                      family = "binomial",
                      alpha = 0, # fitting lasso
                      weights = wt,
                      nfolds = 10) # number of folds for cross-validation
temp.time <- as.numeric(Sys.time() - tmp)

#plot(ridge_cv)

# predict ridge
ridge_pred <- predict(ridge_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
#head(ridge_pred)

ridge_pred <- as.numeric(ridge_pred)
names(ridge_pred) <- rownames(testing_env) #add site names


ridge_pred <- list(ridge_cv, ridge_pred, modelName, temp.time)

#- - - - - - - - - - - - - - - - - - - - -
## MARS ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.mars"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/4

mars_pred_list <- list()

for(no.runs in 1:no.loop.runs){
  
  temp.files.subset <- list.files(path = paste0("./results/",Taxon_name), 
                           pattern = paste0(modelName, no.runs, "_[[:graph:]]*", spID), full.name = T)
  
  lapply(temp.files.subset, load, .GlobalEnv)
  
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
  mars_fit <- caret::train(form = occ ~ .,
                    data = training,
                    method = "earth",
                    metric = "ROC",
                    trControl = mytrControl,
                    tuneGrid = mytuneGrid,
                    thresh = 0.00001)
  stopCluster(cluster)
  registerDoSEQ()
  temp.time <- as.numeric(Sys.time() - tmp)
  #plot(mars_fit)
  
  mars_pred <- predict(mars_fit, testing_env)
  # transform occurrence back into numeric
  mars_pred <- as.character(mars_pred)
  mars_pred[mars_pred=="C0"] <- 0
  mars_pred[mars_pred=="C1"] <- 1
  mars_pred <- as.numeric(mars_pred)
  names(mars_pred) <- rownames(testing_env) #add site names
  #head(mars_pred)
  
  mars_pred_list[[no.runs]] <- list(mars_fit, mars_pred, modelName, temp.time)
}

# average all MARS predictions
mars_pred <- as.data.frame(sapply(mars_pred_list, "[[", 2))
mars_pred <- rowMeans(mars_pred, na.rm=T)
temp.time <- mean(sapply(mars_pred_list, "[[", 4), na.rm=T)

temp.models <- sapply(mars_pred_list, "[[", 1)

mars_pred <- list(temp.models, mars_pred, modelName, temp.time)

mars_pred <- as.numeric(mars_pred)
names(mars_pred) <- rownames(testing_env) #add site names


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
        maxmod <- dismo::maxent(x = data[trainSet, covarsNames],
                                p = data$occ[trainSet],
                                removeDuplicates = FALSE,
                                path = filepath,
                                args = as.character(unlist(grid[n, ])))
      ), "try-error")){
        next
      }
      modpred <- predict(maxmod, data[testSet, covarsNames], args = "outputformat=cloglog")
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
maxmod <- dismo::maxent(x = training[, covarsNames],
                        p = training$occ,
                        removeDuplicates = FALSE,
                        path = paste0("results/", Taxon_name, "/maxent_files"),
                        rgs = param_optim)
temp.time <- as.numeric(Sys.time() - tmp)

maxmod_pred <- predict(maxmod, testing_env)
names(maxmod_pred) <- rownames(testing_env) #add site names
#head(maxmod_pred)

maxmod_pred <- as.numeric(maxmod_pred)
names(maxmod_pred) <- rownames(testing_env) #add site names


maxmod_pred <- list(maxmod, maxmod_pred, modelName, temp.time)

## MaxNet
presences <- training$occ # presence (1s) and background (0s) points
covariates <- training[, 2:ncol(training)] # predictor covariates
tmp <- Sys.time()
set.seed(32639)
mxnet <- maxnet::maxnet(p = presences,
                        data = covariates,
                        regmult = 1, # regularization multiplier
                        maxnet::maxnet.formula(presences, covariates, classes = "default"))
temp.time <- as.numeric(Sys.time() - tmp)

# predicting with MaxNet
maxnet_pred <- predict(mxnet, testing_env, type = c("cloglog"))[, 1]
head(maxnet_pred)

maxnet_pred <- as.numeric(maxnet_pred)
names(maxnet_pred) <- rownames(testing_env) #add site names

maxnet_pred <- list(mxnet, maxnet_pred, modelName, temp.time)

#- - - - - - - - - - - - - - - - - - - - -
## BRT (GBM) ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/4

brt_pred_list <- list()

for(no.runs in 1:no.loop.runs){
  
  temp.files.subset <- list.files(path = paste0("./results/",Taxon_name), 
                                  pattern = paste0(modelName, no.runs, "_[[:graph:]]*", spID), full.name = T)
  
  lapply(temp.files.subset, load, .GlobalEnv)
    
  # calculating the case weights
  prNum <- as.numeric(table(training$occ)["1"]) # number of presences
  bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
  wt <- ifelse(training$occ == 1, 1, prNum / bgNum)
  set.seed(32639)
  
  # settings
  brt <- NULL
  ntrees <- 50
  tcomplexity <- ifelse(prNum < 50, 1, 5)
  lrate <- 0.001
  m <- 0
  
  # start modeling! We use the "try" notation so if a species fails to fit, the loop will continue.
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
    
    tmp <- Sys.time()
  
    if(inherits(try(
      
      brt <- gbm.step(data = training,
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
                      silent = TRUE) # avoid printing the cv results
    ), "try-error")){
      cat("Couldn't fit model", spID, "in the iteration", m, "\n")
    } 
    temp.time <- as.numeric(Sys.time() - tmp)
  }
  if(is.null(brt)){
    next
  }
  
  # Note: model tuning with ~50,000 background points may take ~1h.
  
  #the optimal number of trees (intersect of green and red line in plot)
  brt$gbm.call$best.trees
  
  # get interactions of predictors
  temp.find.int <- gbm.interactions(brt)
  
  # predicting with the best trees
  brt_pred <- predict(brt, testing_env, n.trees = brt$gbm.call$best.trees, type = "response")
  #head(brt_pred)
  
  brt_pred <- as.numeric(brt_pred)
  names(brt_pred) <- rownames(testing_env) #add site names
  
  brt_pred_list[[no.runs]] <- list(brt, brt_pred, modelName, temp.time, 
                                   data.frame("n.trees"=ntrees, "l.rate"=lrate, "t.complexity"=tcomplexity), 
                                   temp.find.int$interactions)
}

# average all BRT predictions
brt_pred <- as.data.frame(sapply(brt_pred_list, "[[", 2))
brt_pred <- rowMeans(brt_pred, na.rm=T)
temp.time <- mean(sapply(brt_pred_list, "[[", 4), na.rm=T)

temp.settings <- rowMeans(as.data.frame(unlist(sapply(brt_pred_list, "[[", 5))), na.rm=T)
temp.settings <- data.frame("parameter"=names(brt_pred_list[[1]][[5]]), "setting"=temp.settings)

temp.interactions <- rowMeans(sapply(brt_pred_list, "[[", 5), na.rm=T)
temp.interactions <- as.data.frame(matrix(temp.interactions, ncol=length(covarsNames)))
rownames(temp.interactions) <- covarsNames
colnames(temp.interactions) <- covarsNames

temp.models <- sapply(brt_pred_list, "[[", 1)

brt_pred <- list(temp.models, brt_pred, modelName, temp.time, temp.settings, temp.interactions)

#- - - - - - - - - - - - - - - - - - - - -
## XGBoost ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/4

xgb_pred_list <- list()

for(no.runs in 1:no.loop.runs){
  
  temp.files.subset <- list.files(path = paste0("./results/",Taxon_name), 
                                  pattern = paste0(modelName, no.runs, "_[[:graph:]]*", spID), full.name = T)
  
  lapply(temp.files.subset, load, .GlobalEnv)
  
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
  xgb_fit <- caret::train(form = occ ~ .,
                   data = training,
                   method = "xgbTree",
                   metric = "ROC",
                   trControl = mytrControl,
                   tuneGrid = mytuneGrid,
                   verbose = TRUE)
  
  stopCluster(cluster)
  registerDoSEQ()
  temp.time <- as.numeric(Sys.time() - tmp)
  
  #plot(xgb_fit)
  
  #print(xgb_fit$bestTune)
  
  xgb_pred <- predict(xgb_fit, testing_env)
  #head(xgb_pred)
  
  # transform occurrence back into numeric
  xgb_pred <- as.character(xgb_pred)
  xgb_pred[xgb_pred=="C0"] <- 0
  xgb_pred[xgb_pred=="C1"] <- 1
  xgb_pred <- as.numeric(xgb_pred)
  names(xgb_pred) <- rownames(testing_env) #add site names
  
  xgb_pred_list[[no.runs]] <- list(xgb_fit, xgb_pred, modelName, temp.time)
}

# average all XGBoost predictions
xgb_pred <- as.data.frame(sapply(xgb_pred_list, "[[", 2))

xgb_pred <- rowMeans(xgb_pred, na.rm=T)

temp.time <- mean(sapply(xgb_pred_list, "[[", 4), na.rm=T)

temp.models <- sapply(xgb_pred_list, "[[", 1)

xgb_pred <- list(temp.models, xgb_pred, modelName, temp.time)

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

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/4

rf_pred_list <- list()
rf_downsampled_pred_list <- list()

for(no.runs in 1:no.loop.runs){
  
  temp.files.subset <- list.files(path = paste0("./results/",Taxon_name), 
                                  pattern = paste0(modelName, no.runs, "_[[:graph:]]*", spID), full.name = T)
  
  lapply(temp.files.subset, load, .GlobalEnv)
  
  # convert the response to factor for producing class relative likelihood
  training$occ <- as.factor(training$occ)
  tmp <- Sys.time()
  set.seed(32639)
  rf <- randomForest(formula = occ ~.,
                     data = training,
                     ntree = 500) # the default number of trees
  temp.time <- as.numeric(Sys.time() - tmp)
  
  # predict with RF
  rf_pred <- predict(rf, testing_env, type = "prob")[, "1"] # prob = continuous prediction
  #head(rf_pred)
  
  rf_pred <- as.numeric(rf_pred)
  names(rf_pred) <- rownames(testing_env) #add site names
  
  rf_pred_list[[no.runs]] <- list(rf, rf_pred, modelName, temp.time)
  
  #plot(rf, main = "RF")
  
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
  temp.time <- as.numeric(Sys.time() - tmp)
  
  #plot(rf_downsample, main = "RF down-sampled")
  
  # predict with RF down-sampled
  rf_downsample_pred <- predict(rf_downsample, testing_env, type = "prob")[, "1"] # prob = continuous prediction
  #head(rf_downsample_pred)
  
  rf_downsample_pred <- as.numeric(rf_downsample_pred)
  names(rf_downsample_pred) <- rownames(testing_env) #add site names
  
  rf_downsample_pred_list[[no.runs]] <- list(rf_downsample, rf_downsample_pred, modelName, temp.time)
}

# average all RF predictions
rf_pred <- as.data.frame(sapply(rf_pred_list, "[[", 2))
rf_pred <- rowMeans(rf_pred, na.rm=T)
temp.time <- mean(sapply(rf_pred_list, "[[", 4), na.rm=T)

temp.models <- sapply(rf_pred_list, "[[", 1)

rf_pred <- list(temp.models, rf_pred, modelName, temp.time)

# average all RF_downsampled predictions
rf_downsample_pred <- as.data.frame(sapply(rf_downsample_pred_list, "[[", 2))
rf_downsample_pred <- rowMeans(rf_downsample_pred, na.rm=T)
temp.time <- mean(sapply(rf_downsample_pred_list, "[[", 4), na.rm=T)

temp.models <- sapply(rf_downsample_pred_list, "[[", 1)

rf_downsample_pred <- list(temp.models, rf_downsample_pred, modelName, temp.time)

# transform occurrence column back to numeric
training$occ <- as.numeric(training$occ)

#- - - - - - - - - - - - - - - - - - - - -
## SVM ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/4

svm_pred_list <- list()

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
  tmp <- Sys.time()
  set.seed(32639)
  svm_e <- e1071::svm(formula = occ ~ .,
                      data = training,
                      kernel = "radial",
                      lass.weights = cwt,
                      probability = TRUE)
  temp.time <- as.numeric(Sys.time() - tmp)
  
  # predicting on test data
  svm_pred <- predict(svm_e, testing_env, probability = TRUE)
  svm_prob <- attr(svm_pred, "probabilities")[,"1"]
  
  # see the first few predictions
  #head(svm_prob)
  
  svm_prob <- as.numeric(svm_prob)
  names(svm_prob) <- rownames(testing_env) #add site names
  
  svm_pred_list[[no.runs]] <- list(svm_e, svm_prob, modelName, temp.time)
}

# average all SVM predictions
svm_pred <- as.data.frame(sapply(svm_pred_list, "[[", 2))
svm_pred <- rowMeans(svm_pred, na.rm=T)
temp.time <- mean(sapply(svm_pred_list, "[[", 4), na.rm=T)

temp.models <- sapply(svm_pred_list, "[[", 1)

svm_pred <- list(temp.models, svm_pred, modelName, temp.time)

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
myExpl <- data.frame(training[, covarsNames])
myRespXY <- training[, c("x", "y")]

tmp <- Sys.time()

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
myBiomodOption <- BIOMOD_ModelingOptions(
  GAM = list( algo = 'GAM_mgcv',
               type = 's_smoother',
               k = -1,
               interaction.level = 0,
               myFormula = as.formula(form),
               family = binomial(link = 'logit'),
               method = 'GCV.Cp',
               optimizer = c('outer','newton'),
               select = FALSE,
               knots = NULL,
               paraPen = NULL,
               control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07, maxit = 200, trace = FALSE,
                              mgcv.tol = 1e-07, mgcv.half = 15, rank.tol = 1.49011611938477e-08,
                              nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0),
                              optim = list(factr=1e+07),
                              newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0), outerPIsteps = 0,
                              idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE, scale.est = "fletcher",
                              edge.correct = FALSE)),
   MAXENT.Phillips = list( path_to_maxent.jar = paste0(here::here(), "/results/", Taxon_name, "/maxent_files"), # change it to maxent directory
                           memory_allocated = 512,
                           background_data_dir = 'default',
                           maximumbackground = 50000,
                           maximumiterations = 200,
                           visible = FALSE,
                           linear = TRUE,
                           quadratic = TRUE,
                           product = TRUE,
                           threshold = TRUE,
                           hinge = TRUE,
                           lq2lqptthreshold = 80,
                           l2lqthreshold = 10,
                           hingethreshold = 15,
                           beta_threshold = -1,
                           beta_categorical = -1,
                           beta_lqp = -1,
                           beta_hinge = -1,
                           betamultiplier = 1,
                           defaultprevalence = 0.5)
)

# models to predict with
mymodels <- c("GLM","GBM","GAM","CTA","ANN","FDA","MARS","RF","MAXENT.Phillips")

# model fitting

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
                                    modeling.id = paste(myRespName,"_Modeling", sep = ""))

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
temp.time <- as.numeric(Sys.time() - tmp)

# project single models
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                  new.env = as.data.frame(testing_env[, covarsNames]),
                                  proj.name = "modeling",
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
#head(myEnProjDF)

biomod_pred <- myEnProjDF[,1]
names(biomod_pred) <- rownames(testing_env)
biomod_pred <- list(myBiomodEM, biomod_pred, modelName, temp.time)
 
  
#- - - - - - - - - - - - - - - - - - - - -
## simple Ensemble model ####
#- - - - - - - - - - - - - - - - - - - - -
# note: this model should be run after all the component 
# models (MaxEnt, Lasso, GAM, BRT, RF down-sampled)

modelName <- "bg.glm"

gm <- scales::rescale(gm_pred[[2]], to = c(0,1))
lasso <- scales::rescale(lasso_pred[[2]], to = c(0,1))
brt <- scales::rescale(brt_pred[[2]], to = c(0,1))
rf <- scales::rescale(rf_pred[[2]], to = c(0,1))
maxt <- scales::rescale(maxmod_pred[[2]], to = c(0,1))

ensm_pred <- rowMeans(cbind(gm, lasso, maxt), na.rm=T) #,brt, rf  !they have different background data

temp.time <- sum(gm_pred[[4]], lasso_pred[[4]], maxmod_pred[[4]]) #,brt_pred[[3]], rf_pred[[3]]

ensm_pred <- list("No model output available. Predictions have to be averaged again if necessary.", ensm_pred, modelName, temp.time)

#- - - - - - - - - - - - - - - - - - - - -
## Saving ####
#- - - - - - - - - - - - - - - - - - - - -

# save all models to calculate model performance later
SDMs <- list(gm_pred, lm1_pred, lm_subset_pred, lasso_pred, ridge_pred, 
     mars_pred, maxmod_pred, maxnet_pred, brt_pred, xgb_pred, svm_pred,
     rf_pred, rf_downsample_pred, biomod_pred, ensm_pred)
names(SDMs) <- c("gm_pred", "lm1_pred", "lm_subset_pred", "lasso_pred", "ridge_pred", 
                "mars_pred", "maxmod_pred", "maxnet_pred", "brt_pred", "xgb_pred", "svm_pred",
                "rf_pred", "rf_downsample_pred", "biomod_pred", "ensm_pred")
#head(SDMs)

save(SDMs, file=paste0(here::here(), "/results/", Taxon_name, "/Predicted_SDMs_", spID, ".RData"))




