#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#    Identify top predictors with MaxEnt    #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8g")) # expand Java memory

gc()
library(tidyverse)
library(here)

library(raster)

library(biomod2) # also to create pseudo-absences

library(dismo) # for MaxEnt and BRT
# download maxent.jar 3.3.3k, and place the file in the
# desired folder
# utils::download.file(url = "https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar", 
#     destfile = paste0(system.file("java", package = "dismo"), 
#         "/maxent.jar"), mode = "wb")  ## wb for binary file, otherwise maxent.jar can not execute

library(rJava)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Crassiclitellata"
speciesNames <- read.csv(file=paste0(here::here(), "/results/Species_list_", Taxon_name, ".csv"))

corMatPearson <- as.matrix(read.csv(file=paste0(here::here(), "/results/corMatPearson_predictors.csv")))
dimnames(corMatPearson)[[1]] <- dimnames(corMatPearson)[[2]]
# based on Valavi et al. 2021: Pearson 0.8
env_exclude <- caret::findCorrelation(corMatPearson, cutoff = 0.8, names=TRUE)
covarsNames <- dimnames(corMatPearson)[[1]][!(dimnames(corMatPearson)[[1]] %in% env_exclude)]
covarsNames <- covarsNames[covarsNames != "x" & covarsNames != "y"]
# exclude based on VIF
env_vif <- read.csv(file=paste0(here::here(), "/results/VIF_predictors.csv"))
env_exclude <- env_vif %>% filter(is.na(VIF)) %>% dplyr::select(Variables) %>% as.character()
covarsNames <- covarsNames[!(covarsNames %in% env_exclude)]
# excluded:
print("=== We excluded the following variables based on VIF and Pearson correlation: ===")
setdiff(env_vif$Variables, covarsNames)

# final predictor variables
print("=== And we kept the following, final predictor variables: ===")
covarsNames

speciesSub <- unique(speciesNames[speciesNames$NumCells_5km >= 10,]$SpeciesID)

#- - - - - - - - - - - - - - - - - - - - -
# note: we will load the datasets before each individual model

# load environmental variables (for projections)
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_5km_normalized.grd"))
#Env_norm <- stack(Env_norm)

# as dataframe
load(paste0(here::here(), "/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df

# define formula for GLM (and biomod)
form <- paste0("occ ~ ", paste0(paste0("s(", covarsNames, ")"), collapse=" + "))


# Calculate the number of cores
no.cores <-  parallel::detectCores()/2 

#- - - - - - - - - - - - - - - - - - - - -
## Prepare model input ####

mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_5km_", Taxon_name, ".csv"))
mySpeciesOcc <- mySpeciesOcc %>% dplyr::select(-year) %>% unique() 

for(spID in speciesSub) { try({
 
  myResp <- as.numeric(mySpeciesOcc[,spID])
  
  # get NAs id
  na.id <- which(is.na(myResp))
  # remove NAs to enforce PA sampling to be done on explanatory rasters
  myResp <- myResp[-na.id]
  
  myRespCoord <- mySpeciesOcc[-na.id,c('x','y')]
  
  myData <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                                expl.var = Env_norm,
                                                resp.xy = myRespCoord,
                                                resp.name = spID,
                                                PA.nb.rep = 1,
                                                PA.nb.absences = 10000,
                                                PA.strategy = "random")
  
  myData <- cbind(myData@data.species, myData@coord, myData@data.env.var)
  myData$SpeciesID <- spID
  myData <- myData %>% rename("occ" = "myData@data.species")
  myData[is.na(myData$occ),"occ"] <- 0
  
  random.rows <- sample(1:nrow(myData), nrow(myData))
  
  validation <- myData[random.rows[1:round(0.2*nrow(myData))], 
                         c("x","y", "SpeciesID", "occ", covarsNames[covarsNames %in% colnames(myData)])]
  
  training <- myData[random.rows[round(0.2*nrow(myData)):nrow(myData)],]
  
  # subset uncorrelated covariates
  training <- training[, c("occ", covarsNames[covarsNames %in% colnames(myData)])]
  
  # save all datasets
  save(training, file=paste0(here::here(), "/results/", Taxon_name, "/_TopPredictor/MaxentData_train_", Taxon_name,"_", spID, ".RData"))
  save(validation, file=paste0(here::here(), "/results/", Taxon_name, "/_TopPredictor/MaxentData_valid_", Taxon_name,"_", spID, ".RData"))
  
})}


#- - - - - - - - - - - - - - - - - - - - -
## Modeling ####

modelName <- "MaxentData"
files <- list.files(path = paste0(here::here(), "/results/",Taxon_name, "/_TopPredictor"), 
                    pattern = paste0(modelName), full.name = T)

for(spID in speciesSub){ try({
  
  # identify and load all relevant data files
  temp.files <- files[stringr::str_detect(files, spID)]
  lapply(temp.files, load, .GlobalEnv)
  
  # "We used five different regularization multipliers (0.5, 1, 2, 3 and 4)
  # in combination with different features (L, LQ, H, LQH, LQHP) to find the 
  # best parameters that maximizes the average AUC-ROC in CV."
  
  # function for simultaneous tuning of MaxEnt regularization multiplier and features
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
  
  # tune MaxEnt parameters
  param_optim <- maxent_param(data = training,
                              k = nfolds,
                              filepath = paste0(here::here(), "/results/", Taxon_name, "/_TopPredictor/maxent_files"))
  
  ## fit a maxent model with the tuned parameters
  maxent <- dismo::maxent(x = training[, covarsNames],
                          p = training$occ,
                          removeDuplicates = FALSE, #remove occurrences that fall into same grid cell (not necessary)
                          path = paste0("results/", Taxon_name, "/maxent_files/", spID), #wanna save files?
                          args = param_optim)
  
  #temp_validation <- dismo::predict(maxent, validation[,colnames(validation) %in% covarsNames], type="response")
  temp_validation <- dismo::predict(maxent, validation[,colnames(validation) %in% covarsNames], type = c("cloglog"))
  temp_validation <- as.numeric(temp_validation)
  names(temp_validation) <- rownames(validation[,colnames(validation) %in% covarsNames]) #add site names
  #head(temp_validation)
  
  temp_model_time <- proc.time()[3] - tmp
  
  
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
  
  # varImp
  temp_varImp <- as.data.frame(maxent@results[str_detect(rownames(maxent@results),"permutation.importance"),])
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
  save(maxent_list, file=paste0(here::here(), "/results/", Taxon_name, "/_TopPredictor/SDM_maxent_", spID, ".RData"))
  
  rm(maxent, maxent_list, param_optim, temp_model_time, temp_predict_time, temp_runs, temp_validation, temp_prediction, temp_varImp)
  gc()
  
  rm(training, validation)

})}


#- - - - - - - - - - - - - - - - - - - - - -
## Calculate variable importance (VI) ####

# create result data frame
var_imp <- data.frame("Predictor"= c("Aridity", "MAP", "MAP_Seas", "MAT", 
                                              "MAT_Seas", "Snow", "Agriculture", "Dist_Urban",
                                              "Forest_Coni", "Forest_Deci", "NDVI", 
                                              "Pastures", "Pop_Dens", "Shrubland", "Aspect",
                                              "Dist_Coast", "Dist_River", "Elev", 
                                              "Slope", "CEC", "Clay.Silt", "Cu", "Hg",
                                              "Moisture", "N", "P", "pH", "SOC", "SoilT"),
                      "maxent"=NA, "Species"=NA)

for(spID in unique(speciesNames[speciesNames$NumCells_5km >= 10,]$SpeciesID)){ try({
  
  print("=====================================")
  print(spID)
  
  temp_sdm <- get(load(file=paste0(here::here(), "/results/", Taxon_name, "/_TopPredictor/SDM_maxent_", spID, ".RData")))

  temp_vi <- temp_sdm[["varImp"]]
  
  # harmonizes predictor column structure
  temp_vi$Predictor <- as.character(temp_vi$Predictor)

  # round variable importance to 3 decimals
  temp_vi[,"maxent"] <- round(temp_vi[,"maxent"], 3)
      
  # add species name
  temp_vi$Species <- spID
  
  var_imp <- rbind(var_imp, temp_vi)

  rm(temp_vi, temp_sdm, temp_varImp)
    
  #var_imp
          
})}

rownames(var_imp) <- 1:nrow(var_imp)
var_imp <- var_imp[!is.na(var_imp$Species),] %>% unique()

var_imp
str(var_imp)

## Save ####
write_csv(var_imp, file=paste0(here::here(), "/results/Variable_importance_MaxEnt_", Taxon_name, ".csv"))
 

#- - - - - - - - - - - - - - - - - -
## Vizualize and get top 10 ####
var_imp <- read.csv(file=paste0(here::here(), "/results/Variable_importance_MaxEnt_", Taxon_name, ".csv"))
var_imp

# load predictor table to get classification of variables
# load the predictor table containing the individual file names
pred_tab <- readr::read_csv(file=paste0(here::here(), "/doc/Env_Predictors_table.csv"))

# transform to long format and add variable categories
var_imp <- var_imp %>%
  left_join(pred_tab %>% dplyr::select(Predictor, Category), by="Predictor")

# add category for clay.silt
var_imp[var_imp$Predictor=="Clay.Silt","Category"] <- "Soil"

# plot VIF
plotVarImp <- ggplot(data=var_imp, aes(x=maxent, y=reorder(Predictor, maxent), fill=Category))+
  geom_boxplot(cex=0.2, outlier.size=0.2)+
  xlab("Variable importance (Permutation importance)")+
  ylab("Predictor")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 5))
plotVarImp

pdf(paste0(here::here(), "/figures/VariableImportance_MaxEnt_", Taxon_name, ".pdf")); plotVarImp; dev.off()

# plot barplot with top 10
plotTopVI <- var_imp %>% dplyr::select(maxent, Predictor, Category) %>%
  group_by(Predictor, Category) %>% summarize_all(mean, na.rm=T) %>% arrange(desc(maxent)) %>%
  ggplot(aes(x=maxent, y=reorder(Predictor, maxent), fill=Category)) + 
  geom_bar(stat="identity") + geom_line(y=length(covarsNames)-9.5)+
  geom_text(aes(label=round(maxent,3)), position=position_dodge(width=0.5), vjust=0.5, hjust=1.1, cex=3)+
  theme_bw()
plotTopVI

pdf(paste0(here::here(), "/figures/VariableImportance_MaxEnt_top10_", Taxon_name, ".pdf")); plotTopVI; dev.off()

# plot barplot with top 10 (based on 75% quartile)
plotTopVI <- var_imp %>% dplyr::select(maxent, Predictor, Category) %>%
  group_by(Predictor, Category) %>% 
	summarize(q75 = quantile(maxent, probs = .75)) %>% 
	arrange(desc(q75)) %>%
  ggplot(aes(x=q75, y=reorder(Predictor, q75), fill=Category)) + 
  xlim(0, 20)+
  geom_bar(stat="identity") + geom_line(y=length(covarsNames)-9.5)+
  geom_text(aes(label=round(q75,3)), position=position_dodge(width=0.5), vjust=0.5, hjust=-0.1, cex=3)+
  theme_bw()
plotTopVI

pdf(paste0(here::here(), "/figures/VariableImportance_MaxEnt_top10_q75_", Taxon_name, ".pdf")); plotTopVI; dev.off()

# plot maps
temp_files <- list.files(paste0(here::here(), "/results/", Taxon_name, "/_TopPredictor"))
temp_files <- temp_files[stringr::str_detect(temp_files, "SDM_maxent_[:graph:]*.RData")]

plots <- lapply(c(1:length(temp_files)), function(m) {try({
  print(temp_files[m]); print(m)
  temp_pred <- get(load(file=paste0(here::here(), "/results/", Taxon_name, "/_TopPredictor/", temp_files[m])))[["prediction"]]
  #print(m)
  spID <- substr(temp_files[m], 12, 21)
  ggplot(data=temp_pred, aes(x=x, y=y, fill=layer))+
    geom_tile()+
    ggtitle(paste0("MaxEnt for ", spID))+
    scale_fill_viridis_c(limits = c(0,1))+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = "right")
})})

lapply(plots, class) # if error, remove that species

#pdf(file=paste0(here::here(), "/figures/DistributionMaps_", Taxon_name, "_", spID, ".pdf"))
png(file=paste0(here::here(), "/figures/DistributionMaps_", Taxon_name, "_MaxEnt.png"),width=3000, height=3000)
do.call(gridExtra::grid.arrange, plots)
dev.off()


#- - - - - - - - - - - - - - - - - -
## Estimate richness ####

# create empty data frame
species_stack <- Env_norm_df %>% dplyr::select(x, y)

# for loop through all species
for(spID in unique(speciesNames[speciesNames$NumCells_5km >= 10,]$SpeciesID)){ try({
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Load probability maps ####
  temp_pred <- get(load(file=paste0(here::here(), "/results/", Taxon_name, "/_TopPredictor/SDM_maxent_", spID, ".RData")))[["prediction"]] #maxent_list
  
  print(paste0(spID, " successfully loaded."))
  
  ## Transform to binary maps 
  
  # extract threshold to define presence/absence
  temp_thresh <- 0.8
  
  # change to binary
  temp_pred[temp_pred$layer>=temp_thresh & !is.na(temp_pred$layer), "layer"] <- 1
  temp_pred[temp_pred$layer<temp_thresh & !is.na(temp_pred$layer), "layer"] <- 0
  
  temp_pred[,paste0(spID,"_", "MaxEnt")] <- temp_pred$layer
  temp_pred <- temp_pred[,c("x","y",paste0(spID,"_", "MaxEnt"))]
  
  ## Stack species binary maps 
  
  # add species dataframe to stacked dataframe
  species_stack <- species_stack %>% full_join(temp_pred, by=c("x","y"))
  
  print(paste0("Added binary prediction of ", spID, " to the species stack."))
  
  rm(temp_thresh, temp_pred)
}, silent=T)}  

head(species_stack)


#- - - - - - - - - - - - - - - - - - - - - -
## Calculate richness ####
species_stack$Richness <- rowSums(species_stack %>% dplyr::select(-x, -y), na.rm=T)

#- - - - - - - - - - - - - - - - - - - - - -
## Save species stack ####
save(species_stack, file=paste0(here::here(), "/results/_Maps/SDM_stack_MaxEnt_binary0.8_", Taxon_name, ".RData"))


#- - - - - - - - - - - - - - - - - - - - - -
## View individual binary maps and species stack ####

species_stack$Richness[species_stack$Richness == 0] <- NA

# species richness
png(file=paste0(here::here(), "/figures/SpeciesRichness_MaxEnt0.8_", Taxon_name, ".png"),width=1000, height=1000)
ggplot(data=species_stack, aes(x=x, y=y, fill=Richness))+
  geom_tile()+
  ggtitle("Species richness (number of species)",
          subtitle=paste0("n = ", ncol(species_stack)-3))+
  scale_fill_viridis_c(na.value = "grey")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(0.1,0.4))
dev.off()


#- - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - -
## Combine results of 10 identification runs  ####
#- - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - -


temp_files <- list.files(paste0(here::here(), "/results/", Taxon_name, "/_TopPredictor"), recursive=T, full.names=T)
temp_files <- temp_files[str_detect(temp_files,"202208") & str_detect(temp_files, ".csv")]
temp_files

if(length(temp_files)==10) print("Now that you have the variable importance estimates from 10 runs, we can identify the top 10 predictors")
if(length(temp_files)!=10) print("####### PLEASE stop here, you did not have variable importance estimates from 10 runs !!! ########")

var_imp <- data.frame("Predictor"= c("Aridity", "MAP", "MAP_Seas", "MAT", 
                                     "MAT_Seas", "Snow", "Agriculture", "Dist_Urban",
                                     "Forest_Coni", "Forest_Deci", "NDVI", 
                                     "Pastures", "Pop_Dens", "Shrubland", "Aspect",
                                     "Dist_Coast", "Dist_River", "Elev", 
                                     "Slope", "CEC", "Clay.Silt", "Cu", "Hg",
                                     "Moisture", "N", "P", "pH", "SOC", "SoilT"),
                      "maxent"=NA, "Species"=NA, "Run"=NA)

for(files in 1:length(temp_files)){ try({
  
  print("=====================================")
  print(temp_files[files])
  
  temp_vi <- read.csv(file=temp_files[files])
  
  # add species name
  temp_vi$Run <- files
  
  var_imp <- rbind(var_imp, temp_vi)
  
  rm(temp_vi)
  
  #var_imp
  
})}

var_imp <- var_imp[!is.na(var_imp$Species),] %>% unique()
rownames(var_imp) <- 1:nrow(var_imp)

head(var_imp)
str(var_imp)

## Save ####
write_csv(var_imp, file=paste0(here::here(), "/results/Variable_importance_MaxEnt_10runs_", Taxon_name, ".csv"))
var_imp <- read.csv(file=paste0(here::here(), "/results/Variable_importance_MaxEnt_10runs_", Taxon_name, ".csv"))

## Plotting ####
# load predictor table to get classification of variables
# load the predictor table containing the individual file names
pred_tab <- readr::read_csv(file=paste0(here::here(), "/doc/Env_Predictors_table.csv"))

# transform to long format and add variable categories
var_imp <- var_imp %>%
  left_join(pred_tab %>% dplyr::select(Predictor, Category), by="Predictor")

# add category for clay.silt
var_imp[var_imp$Predictor=="Clay.Silt","Category"] <- "Soil"

# plot VIF
plotVarImp <- ggplot(data=var_imp , aes(x=maxent, y=reorder(Predictor, maxent), fill=Category))+
  geom_boxplot(cex=0.2, outlier.size=0.2)+
  xlab("Variable importance (Permutation importance)")+
  ylab("Predictor")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 5))
plotVarImp

pdf(paste0(here::here(), "/figures/VariableImportance_MaxEnt_10runs_", Taxon_name, ".pdf")); plotVarImp; dev.off()

# plot barplot with top 10
plotTopVI <- var_imp %>% dplyr::select(maxent, Predictor, Category) %>%
  group_by(Predictor, Category) %>% summarize_all(mean, na.rm=T) %>% arrange(desc(maxent)) %>%
  ggplot(aes(x=maxent, y=reorder(Predictor, maxent), fill=Category)) + 
  geom_bar(stat="identity") + geom_hline(yintercept=length(covarsNames)-9.5, lty=2)+
  geom_text(aes(label=round(maxent,3)), position=position_dodge(width=0.5), vjust=0.5, hjust=1.1, cex=3)+
  theme_bw()
plotTopVI

pdf(paste0(here::here(), "/figures/VariableImportance_MaxEnt_top10_10runs_", Taxon_name, ".pdf")); plotTopVI; dev.off()

# plot barplot with top 10 (based on 75% quartile)
plotTopVI <- var_imp %>% dplyr::select(maxent, Predictor, Category) %>%
  group_by(Predictor, Category) %>% 
  summarize(q75 = quantile(maxent, probs = .75)) %>% 
  arrange(desc(q75)) %>%
  ggplot(aes(x=q75, y=reorder(Predictor, q75), fill=Category)) + 
  xlim(0, 20)+
  geom_bar(stat="identity") + geom_hline(yintercept=length(covarsNames)-9.5, lty=2)+
  geom_text(aes(label=round(q75,3)), position=position_dodge(width=0.5), vjust=0.5, hjust=-0.1, cex=3)+
  theme_bw()
plotTopVI

pdf(paste0(here::here(), "/figures/VariableImportance_MaxEnt_top10_q75_10runs", Taxon_name, ".pdf")); plotTopVI; dev.off()

## count number of times that predictors are in top10 per species ####
top10 <- c()

for(spID in unique(var_imp$Species)){
  for(run in 1:10){
    temp_top10 <- var_imp %>% filter(Species==spID & Run==run) %>% 
      filter(maxent > 0 & !is.na(maxent)) %>%
      top_n(n = 10, wt = maxent) %>% dplyr::select(Predictor, Species)
    top10 <- rbind(top10, temp_top10)
  }
}

top10 <- top10 %>% count(Predictor, Species)

# load the predictor table containing the individual file names
pred_tab <- readr::read_csv(file=paste0(here::here(), "/doc/Env_Predictors_table.csv"))

# transform to long format and add variable categories
top10 <- top10 %>%
  left_join(pred_tab %>% dplyr::select(Predictor, Category), by="Predictor")

# add category for clay.silt
top10[top10$Predictor=="Clay.Silt","Category"] <- "Soil"


# plot barplot with top 10 (based on top10 counts)
plotTop10 <- top10 %>% dplyr::select(n, Predictor, Category) %>%
  group_by(Predictor, Category) %>%  summarize(sum=sum(n)) %>%
  arrange(desc(sum)) %>%
  ggplot(aes(x=sum, y=reorder(Predictor, sum), fill=Category)) + 
  xlim(0, 400)+
  geom_bar(stat="identity") + geom_hline(yintercept=length(covarsNames)-9.5, lty=2)+
  geom_text(aes(label=sum), position=position_dodge(width=0.5), vjust=0.5, hjust=-0.1, cex=3)+
  theme_bw()
plotTop10

pdf(paste0(here::here(), "/figures/VariableImportance_MaxEnt_top10_count10_10runs", Taxon_name, ".pdf")); plotTop10; dev.off()

## same but per species ####
plotTop10 <- top10 %>% dplyr::select(n, Predictor, Category, Species) %>%
  ggplot(aes(x=n, y=reorder(Predictor, desc(Category)), fill=Category)) + 
  xlim(0, 10.5)+
  geom_bar(stat="identity") + geom_hline(yintercept=length(covarsNames)-9.5, lty=2)+
  geom_text(aes(label=n), position=position_dodge(width=0.5), vjust=0.5, hjust=-0.1, cex=2)+
  facet_wrap(vars(Species))+
  theme_bw()+theme(axis.text = element_text(size=3))
plotTop10

pdf(paste0(here::here(), "/figures/VariableImportance_MaxEnt_top10_count10_species_10runs", Taxon_name, ".pdf"), width=20, height=15); plotTop10; dev.off()

# plot barplot with top 10 (based on top10 counts) for species n>=100 records
plotTop10 <- top10 %>% filter(Species %in% speciesNames[speciesNames$NumCells_2km >=100,"SpeciesID"]) %>%
  dplyr::select(n, Predictor, Category) %>%
  group_by(Predictor, Category) %>%  summarize(sum=sum(n)) %>%
  arrange(desc(sum)) %>%
  ggplot(aes(y=sum, x=reorder(Predictor, sum), fill=Category)) + 
  ylim(0, 200)+
  geom_segment(aes(x=reorder(Predictor, sum), xend=reorder(Predictor, sum), y=0, yend=sum), color="black") +
  geom_point(aes(color=Category), size=4, alpha=1) +
  geom_vline(xintercept=length(covarsNames)-10.5, lty=2)+
  coord_flip() +
  theme_bw()
plotTop10

pdf(paste0(here::here(), "/figures/VariableImportance_MaxEnt_top10_count10_10runs_n100_", Taxon_name, ".pdf")); plotTop10; dev.off()

# plot barplot with top 10 (based on top10 counts) for species 10<n<100 records
plotTop10 <- top10 %>% filter(Species %in% unique(speciesNames[speciesNames$NumCells_2km >=10 & speciesNames$NumCells_2km <100,"SpeciesID"])) %>%
  dplyr::select(n, Predictor, Category) %>%
  group_by(Predictor, Category) %>%  summarize(sum=sum(n)) %>%
  arrange(desc(sum)) %>%
  ggplot(aes(y=sum, x=reorder(Predictor, sum), fill=Category)) + 
  #ylim(0, 320)+
  geom_segment(aes(x=reorder(Predictor, sum), xend=reorder(Predictor, sum), y=0, yend=sum), color="black") +
  geom_point(aes(color=Category), size=4, alpha=1) +
  geom_vline(xintercept=length(covarsNames)-9.5, lty=2)+
  coord_flip()+
 theme_bw()
plotTop10

pdf(paste0(here::here(), "/figures/VariableImportance_MaxEnt_top10_count10_10runs_n10-99_", Taxon_name, ".pdf")); plotTop10; dev.off()


## stacked barplot all species
var_imp$Predictor <- factor(var_imp$Predictor, levels=c("MAP", "MAP_Seas", "MAT",
                                                        "Aspect", "Dist_Coast", "Elev", "Dist_River", "Slope",
                                                        "Agriculture","Dist_Urban", "Forest_Coni", "Forest_Deci", "NDVI", "Pastures", "Pop_Dens", "Shrubland",
                                                        "CEC", "Clay.Silt", "Cu", "Hg","Moisture", "N", "P", "pH", "SOC"))
plotAllVI <- ggplot(var_imp %>% dplyr::select(-Run) %>% group_by(Species, Predictor, Category) %>% summarize_all(mean) %>%
         full_join(var_imp %>% dplyr::select(-Run) %>%  filter(Category=="Climate") %>%
                     group_by(Species, Predictor, Category) %>% summarize_all(mean) %>% ungroup() %>%
                     dplyr::select(Species, Category, maxent)  %>%
                     group_by(Species, Category) %>% summarize_all(sum) %>%
                     group_by(Species) %>% top_n(1, maxent) %>% unique() %>% 
                     dplyr::select(-Category) %>% rename("SumClimate"=maxent), by="Species"),
       aes(fill=Predictor, alpha=Predictor, y=maxent, x=reorder(Species, SumClimate))) + 
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  xlab("Species")+
  scale_y_continuous(expand = c(0, 0))+
  scale_alpha_manual(values=c("MAP"=0.75, "MAP_Seas"=0.55, "MAT"=0.35, 
                              "Aspect"=0.85, "Dist_Coast"=0.65, "Elev"=0.55,"Dist_Coast"=0.45, "Elev"=0.35, "Dist_River"=0.25, "Slope"=0.15,
                              "Agriculture"=0.85,"Dist_Urban"=0.65, "Forest_Coni"=0.45, "Forest_Deci"=0.25, "NDVI"=0.75, "Pastures"=0.55,  "Pop_Dens"=0.35, "Shrubland"=0.15, 
                              "CEC"=0.75,"Clay.Silt"=0.65, "Cu"=0.55, "Hg"=0.45,"Moisture"=0.75, "N"=0.65, "P"=0.55, "pH"=0.45, "SOC"=0.35))+
  scale_fill_manual(values=c("MAP"="#F8766D", "MAP_Seas"="#F8766D", "MAT"="#F8766D", 
                             "Aspect"="#00BFC4", "Dist_Coast"="#00BFC4", "Elev"="#00BFC4","Dist_Coast"="#00BFC4", "Elev"="#00BFC4", "Dist_River"="#00BFC4", "Slope"="#00BFC4",
                             "Agriculture"="#7CAE00","Dist_Urban"="#7CAE00", "Forest_Coni"="#7CAE00", "Forest_Deci"="#7CAE00", "NDVI"="#698B22", "Pastures"="#698B22",  "Pop_Dens"="#698B22", "Shrubland"="#698B22", 
                             "CEC"="#C77CFF","Clay.Silt"="#C77CFF", "Cu"="#C77CFF", "Hg"="#C77CFF","Moisture"="#BF3EFF", "N"="#BF3EFF", "P"="#BF3EFF", "pH"="#BF3EFF", "SOC"="#BF3EFF"))+
  theme_bw()+
  theme(legend.position = "bottom")

png(paste0(here::here(), "/figures/VariableImportance_maxent_species_", Taxon_name, ".png"), height=800, width=600); plotAllVI; dev.off()


# select based on at least 15 times always used
top10 %>% filter(n==10) %>% dplyr::select(Predictor) %>% count(Predictor) %>% filter(n>=15)


