#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#             Predict SDMs                  #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# Note: we will load the datasets and models before predicting each individual model

# load models
#load(file=paste0(here::here(), "/sdm/SDM_Models_", spID, ".RData"))

# load environmental space
#Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km_normalized.grd"))

# as dataframe
#load(paste0(here::here(),"/results/EnvPredictor_2km_df_normalized.RData")) #Env_norm_df

for(spID in unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID)){ try({
print(paste0("Species: ", spID))

# read model performance evaluation table (for threshold MaxEnt & saving best model)
mod_eval <- read.csv(file=paste0(here::here(), "/results/", Taxon_name, "/ModelEvaluation_", Taxon_name, "_", spID, ".csv"))

# # define function to rescale raster (for predicted occurrence) between 0 and 1
# fct.rescale <- function(x, x.min, x.max, new.min = 0, new.max = 1) {
#   if(is.null(x.min)) {x.min = min(x, na.rm=T)}
#   if(is.null(x.max)) {x.max = max(x, na.rm=T)}
#   new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
# }

#- - - - - - - - - - - - - - - - - - - - - -
## Plot all maps ####
#- - - - - - - - - - - - - - - - - - - - - -
modelNames <- c("lm1", "lm_subset", "gm", 
                 "lasso", "ridge",
                "mars","maxent","maxnet", "xgb", "brt", "brt2",
                "rf", "rf2", "rf_downsample", "svm", "biomod", "ensm")

plots <- lapply(c(1:length(modelNames)), function(m) {try({
  temp_pred <- get(load(file=paste0(here::here(), "/results/", Taxon_name, "/temp_files/SDM_", modelNames[m],"_", spID, ".RData")))[["prediction"]]
  print(m)
  ggplot(data=temp_pred, aes(x=x, y=y, fill=layer))+
      geom_tile()+
      ggtitle(modelNames[m])+
      scale_fill_viridis_c(limits = c(0,1))+
      theme_bw()+
      theme(axis.title = element_blank(), legend.title = element_blank(),
           legend.position = c(0.1,0.4))
})})

require(gridExtra)
#pdf(file=paste0(here::here(), "/figures/DistributionMaps_", Taxon_name, "_", spID, ".pdf"))
png(file=paste0(here::here(), "/figures/DistributionMaps_", Taxon_name, "_", spID, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots)
dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Save maps ####
#- - - - - - - - - - - - - - - - - - - - - -
# stack all predictions
model_list <- lapply(c(1:length(modelNames)), function(m) {try({
  get(load(file=paste0(here::here(), "/results/", Taxon_name, "/temp_files/SDM_", modelNames[m],"_", spID, ".RData")))[["prediction"]]
})})
names(model_list) <- modelNames

save(model_list, file=paste0(here::here(), "/results/_Maps/SDM_Predictions_", Taxon_name, "_", spID, ".RData"))

#- - - - - - - - - - - - - - - - - - - - - -
## Save best performing model prediction only
# select best model based on TSS
if(is.na(max(mod_eval$tss))){
  best_model <- "ensm"
}else{
  best_model <- mod_eval[mod_eval$tss+mod_eval$roc == max(mod_eval$tss+mod_eval$roc, na.rm=T) & !is.na(mod_eval$model), "model"]}

best_pred <- model_list[[best_model]]

# save best model prediction
save(best_pred, file=paste0(here::here(), "/results/_Maps/SDM_bestPrediction_", Taxon_name, "_", spID, ".RData"))

})}

ggplot(data=best_pred, aes(x=x, y=y, fill=layer))+
  geom_tile()+
  ggtitle(best_model)+
  scale_fill_viridis_c(limits = c(0,1))+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))

