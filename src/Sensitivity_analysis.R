#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Sensitivity analysis              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# Sensitivity analysis: We will check, if SDMs for species with many occurrence 
# records (i.e., present in >200 grid cells) are similar to the SDMs for the 
# same species made with only a subset of available records (i.e., with only 50, 
# 10, and 5 records).

# We will test this for RF, GAM and biomod only.
# This script is based on the following scripts:
# "Create_backgroundData.R"
# "SDMs_Valavi.R"
# ...

#- - - - - - - - - - - - - - - - - - - - - -
## Create background data ####
#- - - - - - - - - - - - - - - - - - - - - -
# environmental (explanatory) variables as raster file
myExpl <- stack(paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))

# response variable (i.e., species occurrences) in wide format
mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_1km_", Taxon_name, ".csv"))
#mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))

# get species with more than equal to 200 records:
temp_species <- speciesNames[speciesNames$NumCells_1km >= 200,]$SpeciesID

# subset species' records
mySpeciesOcc <- mySpeciesOcc[,c("x", "y", temp_species)]

## parallelize
# Calculate the number of cores
no.cores <- detectCores()/2; no.cores

# Initiate cluster used in foreach function
registerDoParallel(no.cores)

#- - - - - - - - - - - - - - - - - - - - - - - 
## For loop through all selected species ####
foreach(myRespName = temp_species, .export = c("mySpeciesOcc"), 
        .packages = c("biomod2", "tidyverse")) %dopar% {
          
          # define response variable index
          myResp <- as.numeric(mySpeciesOcc[,myRespName])
          
          # get NAs id
          na.id <- which(is.na(myResp))
          
          # remove NAs to enforce PA sampling to be done on explanatory rasters
          myResp <- myResp[-na.id]
          myRespCoord <- mySpeciesOcc[-na.id,c('x','y')]
          
          # create summary table for model settings
          model.settings <- data.frame(SpeciesID=myRespName, model="X", strategy="test", 
                                       No.runs=1, No.points=1, min.distance=1, run.time=0)[0,]
          
          ## GLM, GAM ####
          # randomly performs consistently well, excepted when presences are climatically 
          #   biased for which ‘2°far’ is the best method
          # 10,000 PA or a minimum of 10 runs with 1,000 PA with an equal weight for 
          #   presences and absences
          temp.runs <- 1
          temp.number <- 50000
          temp.strategy <- "random"
          
          tmp <- Sys.time()
          bg.glm <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                                  expl.var = myExpl,
                                                  resp.xy = myRespCoord,
                                                  resp.name = myRespName,
                                                  PA.nb.rep = temp.runs,
                                                  PA.nb.absences = temp.number,
                                                  PA.strategy = temp.strategy)
          temp.time <- Sys.time() - tmp
          
          if(temp.runs==1){
            bg.glm <- cbind(bg.glm@data.species, bg.glm@coord, bg.glm@data.env.var)
          }else{
            bg.glm <- cbind(bg.glm@PA, bg.glm@coord, bg.glm@data.env.var)
          }
          bg.glm$SpeciesID <- myRespName
          
          str(bg.glm)
          print("Note: TRUE in PA means presence.")
          
          # add model settings into summary table
          temp.dat <- data.frame(SpeciesID=myRespName, model="GLM.GAM", strategy=temp.strategy,
                                 No.runs=temp.runs, No.points=temp.number, 
                                 min.distance=NA, run.time=temp.time)
          model.settings <- rbind(model.settings, temp.dat)
          
          rm(temp.strategy, temp.runs, temp.number, temp.min.dist, temp.time)
          
          ## CTA, BRT, RF ####	
          # ‘2°far’ performs consistently better with few presences, ‘SRE’ performs 
          #   better with a large number of presences 	
          # same as number of presences, 10 runs when less than 1000 PA with an equal 
          # weight for presences and absences
          if(length(myResp) < 500){
            temp.runs <- 10
            temp.number <- 500
            temp.strategy <- "disk"
            temp.min.dist <- 222222
            
          }else{
            if(length(myResp) < 1000){
              temp.runs <- 10
              temp.number <- length(myResp)
              temp.strategy <- "sre"
              temp.min.dist <- NULL
            }else{
              temp.runs <- 1
              temp.number <- length(myResp)
              temp.strategy <- "sre"
              temp.min.dist <- NULL
            }
          }
          
          tmp <- Sys.time()
          bg.rf <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = myExpl,
                                        resp.xy = myRespCoord,
                                        resp.name = myRespName,
                                        PA.nb.rep = temp.runs,
                                        PA.nb.absences = temp.number,
                                        PA.strategy = temp.strategy, 
                                        PA.dist.min = temp.min.dist)
          temp.time <- Sys.time() - tmp
          
          if(temp.runs==1){
            bg.rf <- cbind(bg.rf@data.species, bg.rf@coord, bg.rf@data.env.var)
            bg.rf$SpeciesID <- myRespName
          }else{
            bg.rf <- cbind(bg.rf@PA, bg.rf@coord, bg.rf@data.env.var)
            bg.rf$SpeciesID <- myRespName
          }
          
          # add model settings into summary table
          temp.dat <- data.frame(SpeciesID=myRespName, model="RF.BRT.CTA",strategy=temp.strategy, 
                                 No.runs=temp.runs, No.points=temp.number, 
                                 min.distance=NA, run.time=temp.time)
          model.settings <- rbind(model.settings, temp.dat)
          
          rm(temp.strategy, temp.runs, temp.number, temp.min.dist, temp.time)
          
          #- - - - - - - - - - - - - - - - - - - 
          ## BIOMOD ####
          temp.runs <- 1
          temp.number <- 50000 * 0.8 # take only 80% of the data as the other ones will
          # be splitted in "Prepare_input.r" in 80% trainign and 20% testing...
          temp.strategy <- "random"
          
          tmp <- Sys.time()
          bg.biomod <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                                     expl.var = myExpl,
                                                     resp.xy = myRespCoord,
                                                     resp.name = myRespName,
                                                     PA.nb.rep = temp.runs,
                                                     PA.nb.absences = temp.number,
                                                     PA.strategy = temp.strategy)
          temp.time <- Sys.time() - tmp
          
          # NOTE: testing data and validation data will be used from bg.glm
          # therefore, we need to save some parameters for later
          bg.biomod <- list("bg.biomod"=bg.biomod, "myResp"=myResp, "myRespCoord"=myRespCoord, "myRespName"=myRespName,
                            "temp.runs"=temp.runs, "temp.number"=temp.number, "temp.strategy"=temp.strategy)
          
          # add model settings into summary table
          temp.dat <- data.frame(SpeciesID=myRespName, model="BIOMOD", strategy=temp.strategy,
                                 No.runs=temp.runs, No.points=temp.number, 
                                 min.distance=NA, run.time=temp.time)
          model.settings <- rbind(model.settings, temp.dat)
          
          rm(temp.strategy, temp.runs, temp.number, temp.min.dist, temp.time)
          
          print(model.settings)
          print("The strategy used for GLM/GAM can be used for other models than listed.")
          
          
          # save model setting summary for later
          save(model.settings, file=paste0(here::here(), "/results/",  Taxon_name, "/_Sensitivity/BackgroundData_modelSettings_", Taxon_name, "_", myRespName, ".RData"))
          
          # save background data
          bg.list <- list(bg.glm, bg.rf, bg.biomod)
          names(bg.list) <- c("bg.glm", "bg.rf", "bg.biomod")
          save(bg.list, file=paste0(here::here(), "/results/", Taxon_name, "/_Sensitivity/SensAna_PA_Env_", Taxon_name, "_", myRespName, ".RData"))
          # write.csv(bg.glm, file=paste0(here::here(), "/results/", Taxon_name, "/_Sensitivity/SensAna_PA_Env_GLM_", Taxon_name, "_", myRespName, ".csv"), row.names = F)
          # write.csv(bg.mars, file=paste0(here::here(), "/results/", Taxon_name, "/_Sensitivity/SensAna_PA_Env_MARS_", Taxon_name, "_", myRespName, ".csv"), row.names = F)
          # write.csv(bg.mda, file=paste0(here::here(), "/results/", Taxon_name, "/_Sensitivity/SensAna_PA_Env_MDA_", Taxon_name, "_", myRespName, ".csv"), row.names = F)
          # write.csv(bg.rf, file=paste0(here::here(), "/results/", Taxon_name, "/_Sensitivity/SensAna_PA_Env_RF_", Taxon_name, "_", myRespName, ".csv"), row.names = F)
        }

# at the end
stopImplicitCluster()

#- - - - - - - - - - - - - - - - - - - - -
## Do SDMs ####
#- - - - - - - - - - - - - - - - - - - - -
# note: we will load the datasets before each individual model

# define formula for GLM (and biomod)
form <- paste0("occ ~ ", paste0(paste0("s(", covarsNames, ")"), collapse=" + "))

# load environmental variables (for projections)
myExpl <- stack(paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))
myExpl <- crop(myExpl, extent_Europe) # crop to Europe
myExpl <- stack(myExpl)

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

# general settings
tmp <- Sys.time()
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

# get running time
temp_time <- Sys.time() - tmp
temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))

# predict (only for model performance)
gm_pred <- predict(gm, validation_env)
#head(gm_pred)

gm_pred <- as.numeric(gm_pred)
names(gm_pred) <- rownames(validation_env) #add site names

temp_runs <- 1

gm_pred <- list(gm, gm_pred, modelName, temp_time, temp_runs)
rm(gm, temp_time)

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

tmp <- Sys.time()

# the base glm model with linear terms
lm1 <- glm(occ ~., data = training, weights = wt, 
           family = binomial(link = "logit"))
summary(lm1)

temp_time <- Sys.time() - tmp
temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))

lm1_pred <- predict(lm1, validation_env)
#head(lm1_pred)

lm1_pred <- as.numeric(lm1_pred)
names(lm1_pred) <- rownames(validation_env) #add site names

temp_runs <- 1

lm1_pred <- list(lm1, lm1_pred, modelName, temp_time, temp_runs)
rm(temp_time)

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

tmp <- Sys.time()
set.seed(32639)

lm_subset <- gam::step.Gam(object = lm1,
                           scope = mod_scope,
                           direction = "both", # up and down-moving of terms
                           data = training, # this is optional
                           parallel = FALSE,   # optional
                           trace = FALSE) #if information is printed or not
temp_time <- Sys.time() - tmp
temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))

summary(lm_subset)

lm_subset_pred <- predict(lm_subset, validation_env)
#head(lm_subset_pred)

lm_subset_pred <- as.numeric(lm_subset_pred)
names(lm_subset_pred) <- rownames(validation_env) #add site names

temp_runs <- 1

lm_subset_pred <- list(lm_subset, lm_subset_pred, modelName, temp_time, temp_runs)
rm(lm_subset, lm1, temp_time, mod_scope)

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

# for the next step, we need the myExpl as data frame
myExpl.df <- raster::rasterToPoints(myExpl)

# now we can predict this quadratic object on the training, testing, and prediction data
# this make two columns for each covariates used in the transformation
training_quad <- predict(quad_obj, newdata = training)
testing_quad <- predict(quad_obj, newdata = validation_env)
predicting_quad <- predict(quad_obj, newdata = as.data.frame(myExpl.df))

# Define vector with appropriate covarsNames for sparse matrix.
# More specifically: create names like this: bio1_1, bio1_2, bio2_1, bio2_2...
sparse_covarsNames <- paste0(rep(covarsNames, each=2), "_", 1:2)

# convert the data.frames to sparse matrices
# select all quadratic (and non-quadratic) columns, except the y (occ)
new_vars <- names(training_quad)[names(training_quad) != "occ"]
training_sparse <- sparse.model.matrix(~. -1, training_quad[, new_vars])
testing_sparse <- sparse.model.matrix( ~. -1, testing_quad[, new_vars])
predicting_sparse <- sparse.model.matrix( ~. -1, predicting_quad[, sparse_covarsNames])

# calculating the case weights
prNum <- as.numeric(table(training_quad$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training_quad$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)

#- - - - - - - - - - - - - - - - - - - - -
## fitting lasso
# Note: The alpha parameter in this model ranges from 0 to 1, where selecting an 
# alpha of 0 leads to ridge regression and 1 to lasso and anything in between 
# is a combination of both called elastic-net.

# The lambda parameter controls regularization – it is the penalty applied 
# to the model’s coefficients. To select the best lambda, internal cross-
# validation was used (in cv.glmnet function).
tmp <- Sys.time()
set.seed(32639)
lasso_cv <- glmnet::cv.glmnet(x = training_sparse,
                              y = training_quad$occ,
                              family = "binomial",
                              alpha = 1, # fitting lasso
                              weights = wt,
                              nfolds = 10) # number of folds for cross-validation
temp_time <- Sys.time() - tmp
temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))

# the cross-validation result
#plot(lasso_cv)
# Note: dashed lines show options for best Lambda parameter
# left: Lambda with minimum deviance (lambda.min)
# right: best Lambda within one standard deviation of left-dashed line (l...1se)

#print("One of the two Lambda (dashed lines) needs to be selected for prediction.")

# predicting with lasso model
lasso_pred <- predict(lasso_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
#head(lasso_pred)

lasso_pred <- as.numeric(lasso_pred)
names(lasso_pred) <- rownames(validation_env) #add site names

# predict for whole environment
#lasso_prediction <- glmnet:::predict.glmnet(object = lasso_cv, newx = predicting_sparse, type = "response", s = "lambda.min")[,1]
lasso_prediction <- raster::predict(lasso_cv, predicting_sparse, type = "response", s = "lambda.min")[,1]
#head(lasso_prediction)
lasso_prediction <- data.frame("site" = names(lasso_prediction), "prediction" = as.numeric(lasso_prediction)) %>%
  full_join(as.data.frame(myExpl.df) %>% mutate("site" = rownames(as.data.frame(myExpl.df))), by = "site")
lasso_prediction <- raster::rasterFromXYZ(lasso_prediction[,c("x", "y", "prediction")])

temp_runs <- 1

lasso_pred <- list(lasso_cv, lasso_pred, modelName, temp_time, temp_runs, lasso_prediction)
rm(lasso_cv, temp_time)

#- - - - - - - - - - - - - - - - - - - - -
## fitting ridge resgression (alpha=0) while identify the right lambda
tmp <- Sys.time()
set.seed(32639)
ridge_cv <- glmnet::cv.glmnet(x = training_sparse,
                              y = training_quad$occ,
                              family = "binomial",
                              alpha = 0, # fitting ridge
                              weights = wt,
                              nfolds = 10) # number of folds for cross-validation
temp_time <- Sys.time() - tmp
temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))

#plot(ridge_cv)

# predict ridge
ridge_pred <- predict(ridge_cv, testing_sparse, type = "response", s = "lambda.min")[,1]
#head(ridge_pred)

ridge_pred <- as.numeric(ridge_pred)
names(ridge_pred) <- rownames(validation_env) #add site names

temp_runs <- 1

# predict for whole environment
ridge_prediction <- raster::predict(ridge_cv, predicting_sparse, type = "response", s = "lambda.min")[,1]
#head(ridge_prediction)
ridge_prediction <- data.frame("site" = names(ridge_prediction), "prediction" = as.numeric(ridge_prediction)) %>%
  full_join(as.data.frame(myExpl.df) %>% mutate("site" = rownames(as.data.frame(myExpl.df))), by = "site")
ridge_prediction <- raster::rasterFromXYZ(ridge_prediction[,c("x", "y", "prediction")])

ridge_pred <- list(ridge_cv, ridge_pred, modelName, temp_time, temp_runs, ridge_prediction)
rm(ridge_cv, temp_time, training_sparse, tratining_quad, testing_sparse)

#- - - - - - - - - - - - - - - - - - - - -
## MARS ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.mars"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/3

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
  tmp <- Sys.time()
  
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
  
  temp_time <- Sys.time() - tmp
  temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))
  #plot(mars_fit)
  
  mars_pred <- predict(mars_fit, validation_env)
  # transform occurrence back into numeric
  mars_pred <- as.character(mars_pred)
  mars_pred[mars_pred=="C0"] <- 0
  mars_pred[mars_pred=="C1"] <- 1
  mars_pred <- as.numeric(mars_pred)
  names(mars_pred) <- rownames(validation_env) #add site names
  #head(mars_pred)
  
  # caluclate variable importance (for later)
  mars_varImp <- caret::varImp(mars_fit, scale=T) #scaled between 0 and 100% 
  mars_varImp <- data.frame("importance" = mars_varImp$importance, "Predictor"=rownames(mars_varImp$importance))
  
  # create raster layer of predictions for whole environmental space
  mars_raster <- raster::predict(myExpl, mars_fit)
  
  mars_pred_list[[no.runs]] <- list(mars_fit, mars_pred, modelName, temp_time, mars_varImp, mars_raster)
}

# average all MARS predictions
mars_pred <- as.data.frame(sapply(mars_pred_list, "[[", 2))
mars_pred <- rowMeans(mars_pred, na.rm=T)
temp_time <- c(mean(as.numeric(sapply(mars_pred_list, "[[", 4)[1,]), na.rm=T), mars_pred_list[[1]][[4]][2])
temp.varImp <- do.call(rbind, lapply(mars_pred_list, "[[", 5)) %>% 
  group_by(Predictor) %>% summarise_all(mean, na.rm=T)
temp.models <- sapply(mars_pred_list, "[[", 1)

temp_runs <- length(mars_pred_list)

## extract predicted probabilities (for later)
temp.prediction <- raster::stack(lapply(mars_pred_list, "[[", 6), raster)
# convert to 0 and 1 (from factors)
values(temp.prediction)[values(temp.prediction)=="1"] <- 0
values(temp.prediction)[values(temp.prediction)=="2"] <- 1
# take mean across layers (i.e., across runs)
temp.prediction <- raster::calc(temp.prediction, fun = mean)

# make predictions for validation numeric
mars_pred <- as.numeric(mars_pred)
names(mars_pred) <- rownames(validation_env) #add site names

mars_pred <- list(temp.models, mars_pred, modelName, temp_time, temp_runs, temp.varImp, temp.prediction)
rm(mars_fit, temp_time, mars_varImp, mars_pred_list, temp.prediction)

# transform occurrence column back to numeric
training$occ <- as.numeric(training$occ)

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
      modpred <- predict(maxmod, data[testSet, covarsNames]) #, args = "outputformat=cloglog"
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

tmp <- Sys.time()
set.seed(32639)

# tune maxent parameters
param_optim <- maxent_param(data = training,
                            k = nfolds,
                            filepath = paste0("results/", Taxon_name, "/maxent_files"))

## fit a maxent model with the tuned parameters
maxmod <- dismo::maxent(x = training[, covarsNames],
                        p = training$occ,
                        removeDuplicates = FALSE, #remove occurrences that fall into same grid cell (not necessary)
                        #path = paste0("results/", Taxon_name, "/maxent_files"), #wanna save files?
                        args = param_optim)
temp_time <- Sys.time() - tmp
temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))

maxmod_pred <- dismo::predict(maxmod, validation_env)
maxmod_pred <- as.numeric(maxmod_pred)
names(maxmod_pred) <- rownames(validation_env) #add site names
#head(maxmod_pred)

# create raster layer of predictions for whole environmental space
maxmod_raster <- raster::predict(myExpl, maxmod)

temp_runs <- 1

maxmod_pred <- list(maxmod, maxmod_pred, modelName, temp_time, temp_runs, maxmod_raster)
rm(maxmod, temp_time)

## MaxNet
presences <- training$occ # presence (1s) and background (0s) points
covariates <- training[, covarsNames] # predictor covariates

tmp <- Sys.time()
set.seed(32639)

maxnet <- maxnet::maxnet(p = presences,
                        data = covariates,
                        regmult = as.numeric(stringr::str_split(param_optim[1], "=")[[1]][2]), # regularization multiplier, a constant, taken from maxmod
                        maxnet::maxnet.formula(presences, covariates, classes = "default"))
temp_time <- Sys.time() - tmp
temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))

# predicting with MaxNet
maxnet_pred <- predict(maxnet, validation_env, type = c("cloglog"))[, 1]
maxnet_pred <- as.numeric(maxnet_pred)
head(maxnet_pred)

names(maxnet_pred) <- rownames(validation_env) #add site names

# create raster layer of predictions for whole environmental space
maxnet_raster <- raster::predict(myExpl, maxnet)

temp_runs <- 1

maxnet_pred <- list(maxnet, maxnet_pred, modelName, temp_time, temp_runs, maxnet_raster)
rm(maxnet, temp_time)

#- - - - - - - - - - - - - - - - - - - - -
## BRT (GBM) ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/3

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
  
  # parameter settings
  brt <- NULL
  ntrees <- 50
  tcomplexity <- ifelse(prNum < 50, 1, 5)
  lrate <- 0.001
  m <- 0
  
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
    
    tmp <- Sys.time()
    
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
                             silent = TRUE) # avoid printing the cv results
    ), "try-error")){
      cat("Couldn't fit model", spID, "in the iteration", m, "\n")
    } 
    temp_time <- Sys.time() - tmp
    temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))
  }
  if(is.null(brt)){
    next
  }
  
  # Note: model tuning with ~50,000 background points may take ~1h.
  
  # the optimal number of trees (intersect of green and red line in plot)
  #brt$gbm.call$best.trees
  
  # get interactions of predictors
  temp.find.int <- dismo::gbm.interactions(brt)
  
  # predicting with the best trees
  brt_pred <- predict(brt, validation_env, n.trees = brt$gbm.call$best.trees, type = "response")
  #head(brt_pred)
  
  # create raster layer of predictions for whole environmental space
  brt_raster <- raster::predict(myExpl, brt)
  
  brt_pred <- as.numeric(brt_pred)
  names(brt_pred) <- rownames(validation_env) #add site names
  
  brt_pred_list[[no.runs]] <- list(brt, brt_pred, modelName, temp_time, 
                                   data.frame("n.trees"=ntrees, "l.rate"=lrate, "t.complexity"=tcomplexity), 
                                   temp.find.int$interactions, brt_raster)
}

# average all BRT predictions
brt_pred <- as.data.frame(sapply(brt_pred_list, "[[", 2))
brt_pred <- rowMeans(brt_pred, na.rm=T)
temp_time <- c(mean(as.numeric(sapply(brt_pred_list, "[[", 4)[1,]), na.rm=T), brt_pred_list[[1]][[4]][2])

temp.settings <- rowMeans(as.data.frame(unlist(sapply(brt_pred_list, "[[", 5))), na.rm=T)
temp.settings <- data.frame("parameter"=names(brt_pred_list[[1]][[5]]), "setting"=temp.settings)

temp_runs <- length(brt_pred_list)

# calculate mean of tuned number of trees (if necessary)
if(ncol(sapply(brt_pred_list, "[[", 6))!=1) {
  temp.interactions <- rowMeans(as.numeric(sapply(brt_pred_list, "[[", 6)), na.rm=T)
}else{
  temp.interactions <- sapply(brt_pred_list, "[[", 6)
}

## extract predicted probabilities (for later)
temp.prediction <- raster::stack(lapply(brt_pred_list, "[[", 7), raster)
# take mean across layers (i.e., across runs)
temp.prediction <- raster::calc(temp.prediction, fun = mean)

temp.interactions <- as.data.frame(matrix(temp.interactions, ncol=length(covarsNames)))
rownames(temp.interactions) <- covarsNames
colnames(temp.interactions) <- covarsNames

temp.models <- sapply(brt_pred_list, "[[", 1)

brt_pred <- list(temp.models, brt_pred, modelName, temp_time, temp_runs, temp.settings, temp.interactions, temp.prediction)
rm(brt, temp_time, brt_pred_list, temp.prediction)

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
  
  tmp <- Sys.time()
  
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
                            silent = TRUE) # avoid printing the cv results
  ), "try-error")){
    cat("Couldn't fit model", spID, "in the iteration", m, "\n")
  } 
  temp_time <- Sys.time() - tmp
  temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))
}
if(is.null(brt2)){
  next
}

# Note: model tuning with ~50,000 background points may take ~1h (or less!).

# the optimal number of trees (intersect of green and red line in plot)
#brt2$gbm.call$best.trees

# get interactions of predictors
temp.find.int <- gbm.interactions(brt2)

# predict to raster (for later)
brt2_raster <- raster::predict(myExpl, brt2, n.trees = brt2$gbm.call$best.trees, type = "response") #maybe remove latter part

# predicting with the best trees
brt2_pred <- predict(brt2, validation_env, n.trees = brt2$gbm.call$best.trees, type = "response")
#head(brt_pred)

temp.settings <- data.frame("parameter"=c("n.trees", "l.rate", "t.complexity"), "setting"=c(ntrees,lrate,tcomplexity))

brt2_pred <- as.numeric(brt2_pred)
names(brt2_pred) <- rownames(validation_env) #add site names

temp_runs <- 1

brt2_pred <- list(brt2, brt2_pred, modelName, temp_time, temp_runs, temp.settings, temp.find.int$interactions, brt2_raster)
rm(brt2, temp_time)

#- - - - - - - - - - - - - - - - - - - - -
## XGBoost ####
#- - - - - - - - - - - - - - - - - - - - -

modelName <- "bg.rf"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
lapply(temp.files, load, .GlobalEnv)

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/3

xgb_pred_list <- list()

for(no.runs in 1:no.loop.runs){
  
  temp.files.subset <- list.files(path = paste0(here::here(), "/results/",Taxon_name), 
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
    nrounds = seq(from = 500, to = 15000, by = 500), # number of boosting iterations
    eta = 0.001,       # shrinkage = learning rate, low = robust to overfitting but slower
    max_depth = 5,     # max. tree depth, default = 6
    subsample = 0.75,  # Subsample Percentage: ratio of training instance -> 75% of data used for fitting, the less the faster
    gamma = 0,         # max. loss reduction required to make further partition on leaf node of tree, the larger the more conservative
    colsample_bytree = 0.8, # Subsample Ratio of Columns when constructing each tree
    min_child_weight = 1    # Minimum Sum of Instance Weight (hessian) needed in a child, 1 = default; if partition step results in leaf node with sum of instance weight < min_child_weight then stop partitioning, the larger the more conservative
  )
  tmp <- Sys.time()
  
  cluster <- makeCluster(2) # you can use all cores of your machine instead e.g. 8
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
  
  temp_time <- Sys.time() - tmp
  temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))
  
  #plot(xgb_fit)
  
  #print(xgb_fit$bestTune)
  
  xgb_pred <- predict(xgb_fit, validation_env)
  #head(xgb_pred)
  
  # transform occurrence back into numeric
  xgb_pred <- as.character(xgb_pred)
  xgb_pred[xgb_pred=="C0"] <- 0
  xgb_pred[xgb_pred=="C1"] <- 1
  xgb_pred <- as.numeric(xgb_pred)
  names(xgb_pred) <- rownames(validation_env) #add site names
  
  # calculate variable importance (for later)
  xgb_varImp <- caret::varImp(xgb_fit, scale=T) #scaled between 0 and 100% 
  xgb_varImp <- data.frame("importance" = xgb_varImp$importance, "Predictor"=rownames(xgb_varImp$importance))
  
  # predict to raster (for later)
  xgb_raster <- raster::predict(myExpl, xgb_fit, type="prob")
   
  # NULL for xgb_fit (to save memory)
  xgb_pred_list[[no.runs]] <- list(NULL, xgb_pred, modelName, temp_time, xgb_varImp, xgb_raster)
}

# average all XGBoost predictions
xgb_pred <- as.data.frame(sapply(xgb_pred_list, "[[", 2))
xgb_pred <- rowMeans(xgb_pred, na.rm=T)
temp_time <- c(mean(as.numeric(sapply(xgb_pred_list, "[[", 4)[1,]), na.rm=T), xgb_pred_list[[1]][[4]][2])
temp.models <- sapply(xgb_pred_list, "[[", 1)
temp.varImp <- do.call(rbind, lapply(xgb_pred_list, "[[", 5)) %>% 
  group_by(Predictor) %>% summarise_all(mean, na.rm=T)

temp_runs <- length(xgb_pred_list)

## extract predicted probabilities (for later)
temp.prediction <- raster::stack(lapply(xgb_pred_list, "[[", 6), raster)
# take mean across layers (i.e., across runs)
temp.prediction <- raster::calc(temp.prediction, fun = mean)

xgb_pred <- list(temp.models, xgb_pred, modelName, temp_time, temp_runs, temp.varImp, temp.prediction)
rm(xgb_fit, temp_time, xgb_pred_list, temp.varImp, xgb_varImp, xgb_raster)

# transform occurrence column back to numeric
training$occ <- as.numeric(training$occ)

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
no.loop.runs <- length(temp.files)/3

rf_pred_list <- list()
rf_downsample_pred_list <- list()

for(no.runs in 1:no.loop.runs){
  
  # load training (and testing) data
  temp.files.subset <- list.files(path = paste0("./results/",Taxon_name), 
                                  pattern = paste0(modelName, no.runs, "_[[:graph:]]*", spID), full.name = T)
  lapply(temp.files.subset, load, .GlobalEnv)
  
  # convert the response to factor for producing class relative likelihood
  training$occ <- as.factor(training$occ)
  tmp <- Sys.time()
  set.seed(32639)
  
  rf <- randomForest::randomForest(formula = occ ~.,
                                   data = training,
                                   importance = T, # assess importance of predictors
                                   ntree = 500) # the default number of trees = 500
  temp_time <- Sys.time() - tmp
  temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))
  
  # predict with RF
  rf_pred <- predict(rf, validation_env, type = "prob")[, "1"] # prob = continuous prediction
  #head(rf_pred)
  
  rf_pred <- as.numeric(rf_pred)
  names(rf_pred) <- rownames(validation_env) #add site names
  
  # predict to raster (for later)
  rf_raster <- raster::predict(myExpl, rf, type="prob")
  
  rf_pred_list[[no.runs]] <- list(rf, rf_pred, modelName, temp_time, rf_raster)
  
  #plot(rf, main = "RF")
  
  ## down-sampling RF
  # for this, background and presence data should be the same number
  prNum <- as.numeric(table(training$occ)["1"]) # number of presences
  bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
  
  # the sample size in each class; the same as presence number
  smpsize <- c("0" = prNum, "1" = prNum)
  tmp <- Sys.time()
  set.seed(32639)
  
  rf_downsample <- randomForest::randomForest(formula = occ ~.,
                                              data = training,
                                              ntree = 1000,
                                              sampsize = smpsize,
                                              importance = T,
                                              replace = TRUE)
  temp_time <- Sys.time() - tmp
  temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))
  
  #plot(rf_downsample, main = "RF down-sampled")
  
  # predict with RF down-sampled
  rf_downsample_pred <- predict(rf_downsample, validation_env, type = "prob")[, "1"] # prob = continuous prediction
  #head(rf_downsample_pred)
  
  rf_downsample_pred <- as.numeric(rf_downsample_pred)
  names(rf_downsample_pred) <- rownames(validation_env) #add site names
  
  # predict to raster (for later)
  rf_downsample_raster <- raster::predict(myExpl, rf_downsample, type="prob")
  
  rf_downsample_pred_list[[no.runs]] <- list(rf_downsample, rf_downsample_pred, modelName, temp_time, rf_downsample_raster)
}

# average all RF predictions
rf_pred <- as.data.frame(sapply(rf_pred_list, "[[", 2))
rf_pred <- rowMeans(rf_pred, na.rm=T)
temp_time <- c(mean(as.numeric(sapply(rf_pred_list, "[[", 4)[1,]), na.rm=T), rf_pred_list[[1]][[4]][2])temp_time <- mean(sapply(rf_pred_list, "[[", 4)[1,], na.rm=T)

temp.models <- sapply(rf_pred_list, "[[", 1)

## extract predicted probabilities (for later)
temp.prediction <- raster::stack(lapply(rf_pred_list, "[[", 5), raster)
# take mean across layers (i.e., across runs)
temp.prediction <- raster::calc(temp.prediction, fun = mean)

temp_runs <- length(rf_pred_list)

rf_pred <- list(temp.models, rf_pred, modelName, temp_time, temp_runs, temp.prediction)
rm(rf, temp_time, rf_pred_list)

# average all RF_downsampled predictions
rf_downsample_pred <- as.data.frame(sapply(rf_downsample_pred_list, "[[", 2))
rf_downsample_pred <- rowMeans(rf_downsample_pred, na.rm=T)
temp_time <- c(mean(as.numeric(sapply(rf_downsample_pred_list, "[[", 4)[1,]), na.rm=T), rf_downsample_pred_list[[1]][[4]][2])

temp.models <- sapply(rf_downsample_pred_list, "[[", 1)

## extract predicted probabilities (for later)
temp.prediction <- raster::stack(lapply(rf_downsample_pred_list, "[[", 5), raster)
# take mean across layers (i.e., across runs)
temp.prediction <- raster::calc(temp.prediction, fun = mean)

temp_runs <- length(rf_downsample_pred_list)

rf_downsample_pred <- list(temp.models, rf_downsample_pred, modelName, temp_time, temp_runs, temp.prediction)
rm(rf_downsample, temp_time, rf_downsample_pred_list)

#- - - - - - - - - - - - - - - - - - - - - 
## Model RF2 for ensemble modelling using consistent background data (bg.glm)

modelName <- "bg.glm"

# identify and load all relevant data files
temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                         pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)

lapply(temp.files, load, .GlobalEnv)

# convert the response to factor for producing class relative likelihood
training$occ <- as.factor(training$occ)
tmp <- Sys.time()
set.seed(32639)
rf2 <- randomForest(formula = occ ~.,
                    data = training,
                    ntree = 500) # the default number of trees
temp_time <- Sys.time() - tmp
temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))

# predict with RF
rf2_pred <- predict(rf2, validation_env, type = "prob")[, "1"] # prob = continuous prediction
#head(rf_pred)

rf2_pred <- as.numeric(rf2_pred)
names(rf2_pred) <- rownames(validation_env) #add site names

# predict to raster (for later)
temp_raster <- raster::predict(myExpl, rf2, type="prob")

temp_runs <- 1

rf2_pred <- list(rf2, rf2_pred, modelName, temp_time, temp_runs, temp_raster)
rm(rf2, temp_time, temp_raster)

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

# how often do we have to run the loop? depending on number of background data simulated
no.loop.runs <- length(temp.files)/3

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
                      class.weights = cwt,
                      probability = TRUE) # allow for probability predictions
  temp_time <- Sys.time() - tmp
  temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))
  
  # predicting on test data
  svm_pred <- e1071:::predict.svm(svm_e, validation_env, probability=TRUE)
  svm_pred <- attr(svm_pred, "probabilities")[,"1"] 
  
  # see the first few predictions
  #head(svm_pred)
  
  temp.names <- names(svm_pred)
  svm_pred <- as.numeric(svm_pred)
  names(svm_pred) <- temp.names #add site names
  
  ## predict for whole environmental space (for later)
  # prepare raster
  Env_df <- raster::rasterToPoints(myExpl, spatial=T)
  
  # predict
  svm_raster <- e1071:::predict.svm(svm_e, Env_df@data[,covarsNames], probability=TRUE)
  Env_df$pred <- attr(svm_raster, "probabilities")[,"1"] 
  
  # make to raster
  svm_raster <- raster::rasterize(Env_df, Env, field="pred")

  svm_pred_list[[no.runs]] <- list(svm_e, svm_pred, modelName, temp_time, svm_raster)
}

# average all SVM predictions
svm_pred <- as.data.frame(sapply(svm_pred_list, "[[", 2))
svm_pred <- rowMeans(svm_pred, na.rm=T)
temp_time <- c(mean(as.numeric(sapply(svm_pred_list, "[[", 4)[1,]), na.rm=T), svm_pred_list[[1]][[4]][2])

temp.models <- sapply(svm_pred_list, "[[", 1)

## extract predicted probabilities (for later)
temp.prediction <- raster::stack(lapply(svm_pred_list, "[[", 5), raster)
# take mean across layers (i.e., across runs)
temp.prediction <- raster::calc(temp.prediction, fun = mean)

temp_runs <- length(svm_pred_list)

svm_pred <- list(temp.models, svm_pred, modelName, temp_time, temp_runs, temp.prediction)
rm(svm_e, temp_time, svm_prob, svm_pred_list)

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
  MAXENT.Phillips = list( path_to_maxent.jar =paste0(here::here(), "/results" )), # change it to maxent directory
  )

# models to predict with
mymodels <- c("GLM","GBM","GAM","CTA","ANN","FDA","MARS","RF","MAXENT.Phillips")

# model fitting
tmp <- Sys.time()
setwd(paste0(here::here(), "/results/", Taxon_name))

set.seed(32639)
myBiomodModelOut <- biomod2::BIOMOD_Modeling(myBiomodData,
                                             models = mymodels,
                                             models.options = myBiomodOption,
                                             NbRunEval = 3,   # 3-fold crossvalidation evaluation run
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
                                               VarImport = 3)    #number of permutations to estimate variable importance
temp_time <- Sys.time() - tmp
temp_time <- c(round(as.numeric(temp_time), 3), units(temp_time))

## NOTE: because biomod output can hardly be stored in list file, we will do calculations based on model output now
# project single models (also needed for ensemble model)
myBiomodProj <- biomod2::BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                           new.env = myExpl,        #column/variable names have to perfectly match with training
                                           proj.name = "modeling",  #name of the new folder being created
                                           selected.models = "all", #use all models
                                           binary.meth = "ROC",     #binary transformation according to criteria
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

# extracting the values for ensemble prediction
myEnProjDF <- as.data.frame(get_predictions(myBiomodEM)[,2]) #for weighted probability mean

# see the first few predictions
# note: the prediction scale of biomod is between 0 and 1000
#head(myEnProjDF)

biomod_pred <- myEnProjDF[,1]
names(biomod_pred) <- rownames(validation_env) #add site names

# Get model evaluation values for later
myBiomodModelEval <- as.data.frame(biomod2::get_evaluations(myBiomodEM),
                                   col.names = colnames(myBiomodEM))

# Calculate variable importance across all PA sets, eval runs and algorithms
# and extract only the one for weighed mean predictions (for later)
temp.varImp <- biomod2::get_variables_importance(myBiomodEM)[, , 2]

# save predictions as raster file
temp_raster <- myBiomodEnProj@proj@val[[2]]

temp_runs <- 1

biomod_pred <- list(myBiomodModelEval, biomod_pred, modelName, temp_time, temp_runs, temp.varImp, temp_raster)
#rm(temp.varImp, myBiomodEM, myBiomodEnProj, myBiomodModelOut, myBiomodProj, myBiomodModelEval, myEnProjDF)

setwd(here::here())  

#- - - - - - - - - - - - - - - - - - - - -
## simple Ensemble model ####
#- - - - - - - - - - - - - - - - - - - - -
# note: this model should be run after all the component 
# models (MaxEnt, Lasso, GAM, BRT, RF down-sampled)

modelName <- "bg.glm"

# scale the predictions between 0 and 1
gm <- scales::rescale(gm_pred[[2]], to = c(0,1))
lasso <- scales::rescale(lasso_pred[[2]], to = c(0,1))
brt <- scales::rescale(brt2_pred[[2]], to = c(0,1))
rf <- scales::rescale(rf2_pred[[2]], to = c(0,1))
maxt <- scales::rescale(maxmod_pred[[2]], to = c(0,1))

# average the predictions
ensm_pred <- rowMeans(cbind(gm, lasso, maxt, brt, rf), na.rm=T) #, brt, rf  !brt and rf have different background data

# sum up computational time
temp_time <- sum(gm_pred[[4]], lasso_pred[[4]], maxmod_pred[[4]], brt2_pred[[4]], rf2_pred[[4]]) 

# model output
ensm_pred <- list("No model output available. Predictions have to be averaged again if necessary.", ensm_pred, modelName, temp_time)

#- - - - - - - - - - - - - - - - - - - - -
## Saving ####
#- - - - - - - - - - - - - - - - - - - - -

# save all models to calculate model performance later
SDMs <- list(gm_pred, lm1_pred, lm_subset_pred, lasso_pred, ridge_pred, 
             mars_pred, maxmod_pred, maxnet_pred, brt_pred, brt2_pred, xgb_pred, svm_pred,
             rf_pred, rf2_pred, rf_downsample_pred, biomod_pred, ensm_pred)
names(SDMs) <- c("gm_pred", "lm1_pred", "lm_subset_pred", "lasso_pred", "ridge_pred", 
                 "mars_pred", "maxmod_pred", "maxnet_pred", "brt_pred", "brt2_pred", "xgb_pred", "svm_pred",
                 "rf_pred", "rf2_pred", "rf_downsample_pred", "biomod_pred", "ensm_pred")
#head(SDMs)


save(SDMs, file=paste0(here::here(), "/sdm/SDM_Models_", spID, ".RData")) 



