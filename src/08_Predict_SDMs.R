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

require(gridExtra)

for(spID in unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID)){ try({
print(paste0("Species: ", spID))

# read model performance evaluation table (for threshold MaxEnt & saving best model)
mod_eval <- read.csv(file=paste0(here::here(), "/results/ModelEvaluation_", Taxon_name, ".csv"))
mod_eval <- mod_eval[mod_eval$species==spID,]

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

temp_files <- list.files(paste0(here::here(), "/results/", Taxon_name, "/temp_files"))
temp_files <- temp_files[stringr::str_detect(temp_files, spID)]

if(length(temp_files)>0){
plots <- lapply(c(1:length(temp_files)), function(m) {try({
  temp_pred <- get(load(file=paste0(here::here(), "/results/", Taxon_name, "/temp_files/", temp_files[m])))[["prediction"]]
  #print(m)
  temp_model <- gsub(paste0("SDM_(.+)_", spID, ".RData"), "\\1", temp_files[m])
  ggplot(data=temp_pred, aes(x=x, y=y, fill=layer))+
      geom_tile()+
      ggtitle(temp_model)+
      scale_fill_viridis_c(limits = c(0,1))+
      theme_bw()+
      theme(axis.title = element_blank(), legend.title = element_blank(),
           legend.position = c(0.1,0.4))
})})

#pdf(file=paste0(here::here(), "/figures/DistributionMaps_", Taxon_name, "_", spID, ".pdf"))
png(file=paste0(here::here(), "/figures/DistributionMaps_", Taxon_name, "_", spID, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots)
dev.off()

while (!is.null(dev.list()))  dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Save maps ####
#- - - - - - - - - - - - - - - - - - - - - -
# stack all predictions
model_list <- lapply(c(1:length(temp_files)), function(m) {try({
  get(load(file=paste0(here::here(), "/results/", Taxon_name, "/temp_files/", temp_files[m])))[["prediction"]]
})})
names(model_list) <- gsub(paste0("SDM_(.+)", ".RData"), "\\1", temp_files)

save(model_list, file=paste0(here::here(), "/results/_Maps/SDM_Predictions_", Taxon_name, "_", spID, ".RData"))

#- - - - - - - - - - - - - - - - - - - - - -
## Save best performing model prediction only
# select best model based on TSS
best_model <- mod_eval[mod_eval$species==spID & !is.na(mod_eval$species) & mod_eval$best==1 & !is.na(mod_eval$best),"model"]
best_pred <- model_list[[paste0(best_model, "_", spID)]]

# save best model prediction
save(best_pred, file=paste0(here::here(), "/results/_Maps/SDM_bestPrediction_", Taxon_name, "_", spID, ".RData"))

}

})}

ggplot(data=best_pred, aes(x=x, y=y, fill=layer))+
  geom_tile()+
  ggtitle(paste0(best_model, " of species ", spID))+
  scale_fill_viridis_c(limits = c(0,1))+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))

