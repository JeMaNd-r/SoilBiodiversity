#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Variable importance              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

## load models ####
#load(file=paste0("./sdm/SDM_Models_", spID, ".RData")) #SDMs

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate variable importance (VI) ####
# create result data frame
var_imp_template <- data.frame("Predictor"= c("Aridity", "MAP", "MAP_Seas", "MAT", 
                                     "MAT_Seas", "Snow", "Agriculture", "Dist_Urban",
                                     "Forest_Coni", "Forest_Deci", "NDVI", 
                                     "Pastures", "Pop_Dens", "Shrubland", "Aspect",
                                     "Dist_Coast", "Dist_River", "Elev", 
                                     "Slope", "CEC", "Clay.Silt", "Cu", "Hg",
                                     "Moisture", "N", "P", "pH", "SOC", "SoilT"))

sdm_names <- c("gm", "lm1", "lm_subset", "lasso", "ridge", "mars", "maxent", "maxnet", "brt", 
               "brt2", "xgb", "svm", "rf", "rf2", "rf_downsample", "biomod", "ensm")

# Calculate the number of cores
no.cores <-  parallel::detectCores()/2 

registerDoParallel(no.cores)
foreach(spID = unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID), 
         .export = c("var_imp_temp"),
         .packages = c("tidyverse", "biomod2", "caret")) %dopar% { try({

var_imp <- var_imp_template

# for loop through all models
for(i in 1:length(sdm_names)){ try({
  temp.model <- sdm_names[i]
  number.models <- 1
 
  print("=====================================")
  print(temp.model)

  temp_sdm <- get(load(file=paste0("./results/", Taxon_name, "/temp_files/SDM_",temp.model,"_", spID, ".RData")))
  #if(exists("temp_sdm")==T){

  temp_vi <- temp_sdm[["varImp"]]
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Ensemble ##
  # see below
  if(temp.model == "ensm_pred"){ next }
   
  #- - - - - - - - - - - - - - - - - - - - - -
  ## ALL ## 
  #- - - - - - - - - - - - - - - - - - - - - -
  # harmonizes predictor column structure
  temp_vi$Predictor <- as.character(temp_vi$Predictor)
  
  # # scale variable importance values between 0 and 10
  # temp.vi[,temp.model] <- scales::rescale(abs(temp_vi[,temp.model]), c(0,10))
   
  # scale to get relative importance: % of sum of importances
  temp_sum <- sum(abs(temp_vi[,temp.model]), na.rm=T)
  temp_vi[,temp.model] <- abs(temp_vi[,temp.model]) / temp_sum
  
  # round variable importance to 3 decimals
  temp_vi[,temp.model] <- round(temp_vi[,temp.model], 3)
   
  var_imp <- dplyr::full_join(var_imp, temp_vi)
  
  # replace NA with 0
  var_imp[,ncol(var_imp)][is.na(var_imp[,ncol(var_imp)])] <- 0
  
  rm(temp_vi, temp_sdm)
  
}, silent=F)}

print("Error in ensm_pred is not an issue")

#- - - - - - - - - - - - - - - - - - - - - -
## Ensemble ####
# average variable importances of individual models
temp_varImp <- var_imp_template
try(temp_varImp <- temp_varImp %>% full_join(var_imp$gm), silent=T)
try(temp_varImp <- temp_varImp %>% full_join(var_imp$lasso), silent=T)
try(temp_varImp <- temp_varImp %>% full_join(var_imp$brt2), silent=T)
try(temp_varImp <- temp_varImp %>% full_join(var_imp$rf2), silent=T)
try(temp_varImp <- temp_varImp %>% full_join(var_imp$maxent), silent=T)

try(var_imp$ensemble <- round(rowMeans(temp_varImp, na.rm=T),3))
rm(temp_varImp)

var_imp

## Save ####
write_csv(var_imp, file=paste0("./results/", Taxon_name, "/Variable_importance_", spID, ".csv"))

})}

stopImplicitCluster()

#- - - - - - - - - - - - - - - - - - - - - -
## Merge all varImp ####

var_imp <- var_imp_template
var_imp$Species <- NA
var_imp$Model <- NA
var_imp$varImp <- NA

temp_files <- list.files(paste0("./results/", Taxon_name))
temp_files <- temp_files[stringr::str_detect(temp_files, "Variable_importance")]
temp_files

for(file in temp_files){
  print(file)
  temp_varImp <- read.csv(file=paste0("./results/", Taxon_name, "/", file))
  temp_varImp$varImp <- rowMeans(temp_varImp %>% dplyr::select(where(is.numeric)), na.rm=T)
  temp_varImp$Species <- substr(file, 21, 30)
  var_imp <- full_join(var_imp, temp_varImp)
}

var_imp <- var_imp[!is.na(var_imp$Species),] %>% unique()
str(var_imp)

write.csv(var_imp, file=paste0("./results/Variable_importance_", Taxon_name,".csv"), row.names = F)


#- - - - - - - - - - - - - - - - - - - - - -
## LASSO and RIDGE varImp - long and accurrate way... unfinished ####
# # define background dataset (for testing data)
# modelName <- SDMs[[i]][[3]]
# 
# # identify and load all relevant data files
# temp.files <- list.files(path = paste0("./results/",Taxon_name), 
#                          pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
# lapply(temp.files, load, .GlobalEnv)
# 
# quad_obj <- make_quadratic(training, cols = covarsNames)
# training_quad <- predict(quad_obj, newdata = training)
# testing_quad <- predict(quad_obj, newdata = validation_env)
# new_vars <- names(training_quad)[names(training_quad) != "occ"]
# testing_sparse <- sparse.model.matrix( ~. -1, testing_quad[, new_vars])
# 
# # load prediction
# prediction <-  prediction(SDMs[[temp.model]][[1]], validation_pa)
# 
# perf <- ROCR::performance(prediction, measure = "auc")
# print(perf@y.values)
# #