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

# read model performance evaluation table (for threshold MaxEnt & saving best model)
mod_eval <- read.csv(file=paste0(here::here(), "/results/ModelEvaluation_", Taxon_name, "_", spID, ".csv"))

# define function to rescale raster (for predicted occurrence) between 0 and 1
fct.rescale <- function(x, x.min, x.max, new.min = 0, new.max = 1) {
  if(is.null(x.min)) {x.min = min(x)}
  if(is.null(x.max)) {x.max = max(x)}
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

#- - - - - - - - - - - - - - - - - - - - - -
## Start predicting ####
#- - - - - - - - - - - - - - - - - - - - - -

# use par() only if using plot() function
#par(mfrow=c(4,4))

#- - - - - - - - - - - - - - - - - - - - - -
## GAM ####
gm <- SDMs[["gm_pred"]][[1]]
gm_pred <- raster::predict(Env, gm)

# rescale between 0 and 1
gm_pred <- fct.rescale(gm_pred, x.min = gm_pred@data@min, x.max = gm_pred@data@max)

#plot(gm_pred, main = "GAM")
# gm <- ggplot(data=data.frame(rasterToPoints(gm_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("GAM")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## GLM ####
lm1 <- SDMs[["lm1_pred"]][[1]]
lm1_pred <- raster::predict(Env, lm1)
# rescale between 0 and 1
lm1_pred <- fct.rescale(lm1_pred, x.min = lm1_pred@data@min, x.max = lm1_pred@data@max)

#plot(lm1_pred, main = "GLM")
# lm1 <- ggplot(data=data.frame(rasterToPoints(lm1_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("GLM")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

# second GLM
lm_subset <- SDMs[["lm_subset_pred"]][[1]]
lm_subset_pred <- raster::predict(Env, lm_subset)
lm_subset_pred <- fct.rescale(lm_subset_pred, x.min = lm_subset_pred@data@min, x.max = lm_subset_pred@data@max)

#plot(lm_subset_pred, main = "GLM subset")
# lm_subset <- ggplot(data=data.frame(rasterToPoints(lm_subset_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("GLM subset")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## Lasso ####
lasso_pred <- SDMs[["lasso_pred"]][[6]]

#plot(lasso_pred, main = "Lasso regression")
# lasso <- ggplot(data=data.frame(rasterToPoints(lasso_pred)), aes(x=x, y=y, fill=prediction))+
#   geom_tile()+
#   ggtitle("Lasso regression")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## Ridge ####
ridge_pred <- SDMs[["ridge_pred"]][[6]]

#plot(ridge_pred, main = "Ridge regression")
# ridge <- ggplot(data=data.frame(rasterToPoints(ridge_pred)), aes(x=x, y=y, fill=prediction))+
#   geom_tile()+
#   ggtitle("Ridge regression")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## MARS ####
# load pre-calculated predictions
mars_pred <- SDMs[["mars_pred"]][[7]]
#mars_pred

# define background dataset (for testing data)
modelName <- SDMs[["mars_pred"]][[3]]

#plot(mars_pred, main = "MARS") #, sub ="threshold = 0.00001 (pre-defined)")
# mars <- ggplot(data=data.frame(rasterToPoints(mars_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("MARS")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## MaxEnt ####
maxent_pred <- SDMs[["maxmod_pred"]][[6]]
#maxent_pred <- fct.rescale(maxmod_pred, x.min = maxent_pred@data@min, x.max = maxent_pred@data@max)

# get threshold from model evaluation table
maxent_thresh <- mod_eval[mod_eval$species==spID & mod_eval$model=="maxmod_pred", "thres.maxTSS"]

#plot(maxent_pred, main = "MaxEnt")
# maxent <- ggplot(data=data.frame(rasterToPoints(maxent_pred>=maxent_thresh)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("MaxEnt", subtitle=paste0("Threshold = ", maxent_thresh))+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))
# print("Note: If the MaxEnt prediction looks bad, it is maybe because of wrong (or missing) predictors...")

#- - - - - - - - - - - - - - - - - - - - - -
# MaxNet ####
maxnet_pred <- SDMs[["maxnet_pred"]][[6]]
#maxnet_pred <- fct.rescale(maxnet_pred, x.min = maxnet_pred@data@min, x.max = maxnet_pred@data@max)
# Note: if max and min are changed, map fits better to GLM etc. ...

# get threshold from model evaluation table
maxnet_thresh <- mod_eval[mod_eval$species==spID & mod_eval$model=="maxnet_pred", "thres.maxTSS"]

#plot(maxnet_pred, main = "MaxNet")
# maxnet <- ggplot(data=data.frame(rasterToPoints(maxnet_pred>=maxnet_thresh)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("MaxNet", subtitle=paste0("Threshold = ", maxnet_thresh))+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## BRT ####
# extract already calculated predictions
brt_pred <- SDMs[["brt_pred"]][[8]]
#brt_pred

# rescale between 0 and 1
brt_pred <- fct.rescale(brt_pred, x.min = brt_pred@data@min, x.max = brt_pred@data@max)

#plot(brt_pred, main = "BRT (GBM)")
# brt <- ggplot(data=data.frame(rasterToPoints(brt_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("BRT (GBM)")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

## BRT 2 (for ensemble model)
# extract already calculated predictions
brt2_pred <- SDMs[["brt2_pred"]][[8]]
#brt_pred

# rescale between 0 and 1
brt2_pred <- fct.rescale(brt2_pred, x.min = brt2_pred@data@min, x.max = brt2_pred@data@max)

#plot(brt2_pred, main = "BRT (GBM)")
# brt2 <- ggplot(data=data.frame(rasterToPoints(brt2_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("BRT 2 (GBM)")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## XGBoost ####
xgb_pred <- SDMs[["xgb_pred"]][[7]]

names(xgb_pred) <- rownames(Env) #add site names

## scaling not necessary if predict(..., type="prob"), see SDM script
# # rescale between 0 and 1
# xgb_pred <- fct.rescale(xgb_pred, x.min = xgb_pred@data@min, x.max = xgb_pred@data@max)

#plot(xgb_pred, main = "XGBoost")
# xgb <- ggplot(data=data.frame(rasterToPoints(xgb_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("XGBoost")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## RF ####
rf_pred <- SDMs[["rf_pred"]][[6]]
#rf_pred

## scaling not necessary if predict(..., type="prob"), see SDM script
# # rescale between 0 and 1
# rf_pred <- fct.rescale(rf_pred, x.min = rf_pred@data@min, x.max = rf_pred@data@max)

#plot(rf_pred, main = "RF")
# rf <- ggplot(data=data.frame(rasterToPoints(rf_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("RF")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

## RF 2
rf2_pred <- SDMs[["rf2_pred"]][[6]]
#rf_pred

## scaling not necessary if predict(..., type="prob"), see SDM script
# # rescale between 0 and 1
# rf_pred <- fct.rescale(rf_pred, x.min = rf_pred@data@min, x.max = rf_pred@data@max)

#plot(rf_pred, main = "RF")
# rf2 <- ggplot(data=data.frame(rasterToPoints(rf2_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("RF for Ensemble")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## RF downsampled ####
rf_downsample_pred <- SDMs[["rf_downsample_pred"]][[6]]
#rf_downsample_pred

## scaling not necessary if predict(..., type="prob"), see SDM script
# # rescale between 0 and 1
# rf_downsample_pred <- fct.rescale(rf_downsample_pred, x.min = rf_downsample_pred@data@min, x.max = rf_downsample_pred@data@max)

#plot(rf_downsample_pred, main = "RF downsampled")
# rf_downsample <- ggplot(data=data.frame(rasterToPoints(rf_downsample_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("RF downsampled")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## SVM ####
svm_pred <- SDMs[["svm_pred"]][[6]]

## scaling not necessary if predict(..., probability=TRUE), see SDM script
# # rescale between 0 and 1
# rf_downsample_pred <- fct.rescale(rf_downsample_pred, x.min = rf_downsample_pred@data@min, x.max = rf_downsample_pred@data@max)

#plot(svm_pred, main = "SVM")
# svm <- ggplot(data=data.frame(rasterToPoints(svm_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("SVM")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## BIOMOD ####
biomod_pred <- SDMs[["biomod_pred"]][[7]]

# rescale between 0 and 1
biomod_pred <- fct.rescale(biomod_pred, x.min = biomod_pred@data@min, x.max = biomod_pred@data@max)
names(biomod_pred) <- "layer"

#plot(biomod_pred, main = "Biomod")
# biomod <- ggplot(data=data.frame(rasterToPoints(biomod_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("Biomod")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## simple Ensemble model ####
# create Raster stack of the models to be merged
ensm_pred <- raster::stack(gm_pred, maxmod_pred) #lasso_pred, brt2_pred, rf2_pred

# average the predictions
ensm_pred <- calc(ensm_pred, fun = mean, na.rm = T)

#plot(ensm_pred, main = "Biomod")
# ensm <- ggplot(data=data.frame(rasterToPoints(ensm_pred)), aes(x=x, y=y, fill=layer))+
#   geom_tile()+
#   ggtitle("Ensemble")+
#   scale_fill_viridis_c(limits = c(0,1))+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position = c(0.1,0.4))

#- - - - - - - - - - - - - - - - - - - - - -
## Plot all maps ####
#- - - - - - - - - - - - - - - - - - - - - -
modelNames <- c("lm1", "lm_subset", "gm", 
                # "lasso", "ridge",
                "mars","maxent","maxnet", "xgb", 
                "rf", "rf_downsample", "svm", "biomod", "ensm")

model_list <- list(lm1_pred, lm_subset_pred, gm_pred, 
                   # "lasso", "ridge",
                   mars_pred, maxent_pred, 
                   maxnet_pred,  xgb_pred, rf_pred, rf_downsample_pred, 
                   svm_pred, biomod_pred, ensm_pred)
names(model_list) <- c("lm1_pred", "lm_subset_pred", "gm_pred", 
                       # "lasso", "ridge",
                       "mars_pred", 
                       "maxent_pred", "maxnet_pred", "xgb_pred", 
                       "rf_pred", "rf_downsample_pred", "svm_pred", "biomod_pred","ensm_pred")

plots <- lapply(1:length(modelNames), function(x) {
  temp.pred <- model_list[[paste0(modelNames[x], "_pred")]]
  
  # exception for maxent/maxnet: plot only above threshold
  if(modelNames[x] == "maxent"){
    # get threshold from model evaluation table
    temp_thresh <- mod_eval[mod_eval$species==spID & mod_eval$model=="maxmod_pred", "thres.maxTSS"]
    
    ggplot(data=data.frame(rasterToPoints(temp.pred >= temp_thresh)), aes(x=x, y=y, fill=layer))+
      geom_raster()+
      ggtitle(modelNames[x], subtitle=paste0("Threshold = ", temp_thresh))+
      scale_fill_viridis_c(limits = c(0,1))+
      theme_bw()+
      theme(axis.title = element_blank(), legend.title = element_blank(),
            legend.position = c(0.1,0.4))
  } else {
    if( modelNames[x] == "maxnet"){
      # get threshold from model evaluation table
      temp_thresh <- mod_eval[mod_eval$species==spID & mod_eval$model=="maxnet_pred", "thres.maxTSS"]
      
      ggplot(data=data.frame(rasterToPoints(temp.pred >= temp_thresh)), aes(x=x, y=y, fill=layer))+
        geom_raster()+
        ggtitle(modelNames[x], subtitle=paste0("Threshold = ", temp_thresh))+
        scale_fill_viridis_c(limits = c(0,1))+
        theme_bw()+
        theme(axis.title = element_blank(), legend.title = element_blank(),
              legend.position = c(0.1,0.4))
  } else {
  
  ggplot(data=data.frame(rasterToPoints(temp.pred)), aes(x=x, y=y, fill=layer))+
    geom_raster()+
    ggtitle(modelNames[x])+
    scale_fill_viridis_c(limits = c(0,1))+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.1,0.4))
}}})

require(gridExtra)
#pdf(file=paste0(here::here(), "/figures/DistributionMaps_", Taxon_name, "_", spID, ".pdf"))
do.call(grid.arrange, plots)
dev.off()

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

#- - - - - - - - - - - - - - - - - - - - - -
## Save best performing model prediction only
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
