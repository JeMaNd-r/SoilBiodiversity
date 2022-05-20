#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#load(file=paste0(here::here(), "/sdm/SDM_Models_", spID, ".RData"))

#- - - - - - - - - - - - - - - - - - - - -
## Model performance ####
#- - - - - - - - - - - - - - - - - - - - -

mod_eval_template <- data.frame(species = NA, # name of the species (more important for later)
                       roc = NA, # AUC roc
                       prg = NA, # AUC prg
                       cor = NA, # correlation between predicted & true values
                       tss = NA, # True Skills Statisti
                       kappa = NA, # Cohen's kappa
                       thres.maxTSS = NA, # Threshold (presence vs. absence) at max. TSS (max(se+sp))
                       model = NA, # name of the SDM algorithm
                       time_model = NA, # mean computational time for model training
                       time_predict = NA, # mean computational time for prediction
                       bg= NA, # background dataset used
                       no.runs = NA, #how many background data runs to we have
                       n_vali_presence = NA, # number of validation presence records
                       n_vali_absence = NA) # number of validation absence records

# function to calculate TSS
opt.cut <- function(perf, pred){
  cut.ind <- mapply(FUN=function(x, y, p){
    d <- (x - 0)^2 + (y-1)^2
    ind <- which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

sdm_names <- c("gm", "lm1", "lm_subset", "lasso", "ridge", "mars", "maxent", "maxnet", "brt", 
               "brt2", "xgb", "svm", "rf", "rf2", "rf_downsample", "biomod", "ensm")


for(spID in unique(speciesNames[speciesNames$NumCells_2km >= 5,]$SpeciesID)){ try({
print(paste0("Species: ", spID))
mod_eval <- mod_eval_template

for(i in 1:length(sdm_names)){ try({ 
  temp.model <- sdm_names[i]
  number.models <- 1
 
  print("=====================================")
  print(temp.model)

  temp_sdm <- get(load(file=paste0(here::here(), "/results/", Taxon_name, "/temp_files/SDM_",temp.model,"_", spID, ".RData")))
  #if(exists("temp_sdm")==T){  

  #- - - - - - - - - - - - - - - - - - - 
  ## calculate statistics for BIOMOD ####
  if(temp.model=="biomod") {

    myBiomodModelEval <- temp_sdm[["evaluation"]]
    
    # define background dataset (for testing data)
    modelName <- temp_sdm[["bg_data"]]
    
    # identify and load all relevant data files
    temp.files <- list.files(path = paste0(here::here(),"/results/",Taxon_name), 
                             pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
    lapply(temp.files, load, .GlobalEnv)
    
    prediction <- temp_sdm[["validation"]]
    
    # save in output data frame
    mod_eval[i,]$model <- sdm_names[i]
    mod_eval[i,]$species <- spID
     
    try(mod_eval[i,]$time_model <- paste0(temp_sdm[["time_model"]][1], " ", temp_sdm[["time_model"]][2]), silent=T)
    try(mod_eval[i,]$time_predict <- paste0(temp_sdm[["time_predict"]][1], " ", temp_sdm[["time_predict"]][2]), silent=T)
    mod_eval[i,]$bg <- modelName
    mod_eval[i,]$no.runs <- 1
    try(mod_eval[i,]$n_vali_presence <- sum(validation$occ), silent=T)
    try(mod_eval[i,]$n_vali_absence <- (length(validation$occ) - sum(validation$occ)), silent=T)
    
    mod_eval[i,]$tss <- myBiomodModelEval[2,1]
    mod_eval[i,]$roc <- myBiomodModelEval[3,1]
    mod_eval[i,]$kappa <- myBiomodModelEval[1,1]
    mod_eval[i,]$thres.maxTSS <- myBiomodModelEval[2,2] / 1000    

    rm(prg, temp.model, modelName, precrec_obj, validation, temp.tss, temp.kappa, temp.tresh, myBiomodModelEval, prediction, modelName)
    
  }else{ # all except BIOMOD modeling
    
    # define background dataset (for testing data)
    modelName <- temp_sdm[["bg_data"]]
    
    # identify and load all relevant data files
    temp.files <- list.files(path = paste0(here::here(),"/results/",Taxon_name), 
                             pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
    
    if(length(temp.files)>2){
      validation_mean <- data.frame("x"=NA, "y"=NA, "occ"=NA)[0,]
      for(l in 1:(length(temp.files)/2)){
         temp.files.subset <- list.files(path = paste0("./results/",Taxon_name), 
                                         pattern = paste0(modelName, l, "_[[:graph:]]*", spID), full.name = T)
  
         lapply(temp.files.subset, load, .GlobalEnv)
         
         validation$sites <- rownames(validation)
         validation_mean <- rbind(validation_mean, validation[,c("x", "y", "occ", "sites", "SpeciesID")])
      }
       validation <- validation_mean %>% group_by(sites, SpeciesID) %>% 
                        summarise(across(everything(), mean)) %>% as.data.frame()
       rownames(validation) <- validation$sites
       
       
       # load prediction
       prediction_raw <- temp_sdm[["validation"]]
       prediction <- as.numeric(prediction_raw$occ)
       names(prediction) <- prediction_raw$sites
 
     }else{
       lapply(temp.files, load, .GlobalEnv)

       # load prediction
       prediction <- temp_sdm[["validation"]]
     }
    
    
      # define validation
      validation <- validation[,c("x","y","SpeciesID", "occ")] %>% mutate("site" = rownames(validation)) %>%
        filter(site %in% names(prediction)) 
      validation <- validation[match(names(prediction), validation$site),] %>% dplyr::select(occ)
    
      if(sum(validation$occ, na.rm=T)==0) print(paste0("For model ", temp.model, ": There is no presence records in the validation data. Too few occurrence records."))

      # try loop (not stopping if error)
      try({
        ## ROC and PR
        # calculate area under the ROC and PR curves
        precrec_obj <- precrec::evalmod(scores = prediction, labels = validation[,"occ"])      

        ## PRG
        # calculate the PRG curve
        prg_curve <- prg::create_prg_curve(labels = validation[,"occ"], pos_scores = prediction)
      
        #- - - - - - - - - - - - - - - - - - - 
        # scale predictions between 0 and 1 
        prediction <- scales::rescale(prediction, to = c(0,1))
      
        ## TSS and Kappa
        roc.pred <- ROCR::prediction(as.numeric(prediction[names(prediction) %in% rownames(validation)]), validation[rownames(validation) %in% names(prediction),"occ"])
        roc.perf <- ROCR::performance(roc.pred, measure = "tpr", x.measure = "fpr")
      
        # calculate Sensitivity and Specificity
        sen_spe <- opt.cut(roc.perf, roc.pred)
      
        # calculate TSS
        temp.tss <- as.numeric(sen_spe [1,1] +  sen_spe[2,1] - 1 )
      
        # calculate Kappa
        mod.object <- sdm::evaluates(x = validation[,"occ"], p = prediction)
        temp.kappa <- mod.object@threshold_based[mod.object@threshold_based$criteria=="max(se+sp)", "Kappa"] #kappa at max(se+sp)
        temp.tresh <- mod.object@threshold_based[mod.object@threshold_based$criteria=="max(se+sp)", "threshold"] #threshold max(se+sp)

    }, silent=T)
    
    # save in output data frame
    mod_eval[i,]$species <- spID
    try(mod_eval[i,]$roc <- precrec::auc(precrec_obj)[1,4], silent=T)
    try(mod_eval[i,]$prg <- prg::calc_auprg(prg_curve), silent=T)
    try(mod_eval[i,]$cor <- cor(prediction,  validation[,"occ"]), silent=T)
    try(mod_eval[i,]$tss <- temp.tss, silent=T)
    try(mod_eval[i,]$kappa <- temp.kappa, silent=T)
    try(mod_eval[i,]$thres.maxTSS <-  temp.tresh, silent=T) # threshold at max(se+sp)
    mod_eval[i,]$model <- sdm_names[i]
    try(mod_eval[i,]$time_model <- paste0(temp_sdm[["time_model"]][1], " ", temp_sdm[["time_model"]][2]), silent=T)
    try(mod_eval[i,]$time_predict <- paste0(temp_sdm[["time_predict"]][1], " ", temp_sdm[["time_predict"]][2]), silent=T)
    mod_eval[i,]$bg <- modelName
    mod_eval[i,]$no.runs <- temp_sdm[["runs"]]
    try(mod_eval[i,]$n_vali_presence <- sum(validation$occ), silent=T)
    try(mod_eval[i,]$n_vali_absence <- (length(validation$occ) - sum(validation$occ)), silent=T)
  
    rm(temp.kappa, temp.tss, temp.tresh, sen_spe, temp.model, modelName, precrec_obj, validation, temp.files, prediction)
  }

 #} # end if exists loop
 rm(temp_sdm)
})}

mod_eval

#- - - - - - - - - - - - - - - - - - - 
## Save ####

write.csv(mod_eval, file=paste0(here::here(), "/results/", Taxon_name,"/ModelEvaluation_", Taxon_name, "_", spID, ".csv"), row.names = F)
rm(mod_eval)
})}

#- - - - - - - - - - - - - - - - - - - 
## Combine all evaluation tables ####

mod_eval <- mod_eval_template

temp_files <- list.files(paste0(here::here(), "/results/", Taxon_name))
temp_files <- temp_files[stringr::str_detect(temp_files, "ModelEvaluation")]
temp_files

for(file in temp_files){
  temp_eval <- read.csv(file=paste0(here::here(), "/results/", Taxon_name, "/", file))
  mod_eval <- rbind(mod_eval, temp_eval)
}

mod_eval <- mod_eval[!is.na(mod_eval$species),] %>% unique()
str(mod_eval)

# identify best performing model
best_mods <- mod_eval %>% filter(!is.na(tss) & !is.na(roc)) %>% group_by(species) %>% slice_max(n=1, order_by=tss+roc)

# is there any species missing?
print("Is there any species' best algorithm missing?")
setdiff(unique(mod_eval$species), unique(best_mods$species))
print("No? Then please continue. If yes, think about what to add to code to select exisiting prediction.")

# add missing species: take ensm as best model
for(spID in c(setdiff(unique(mod_eval$species), unique(best_mods$species)))){
   best_mods <- rbind(best_mods, mod_eval[mod_eval$species==spID & mod_eval$model=="ensm" & !is.na(mod_eval$model),])
}

print("Which species is still missing?")
setdiff(unique(mod_eval$species), unique(best_mods$species))

# add missing species if ensm failed
for(spID in c(setdiff(unique(mod_eval$species), unique(best_mods$species)))){
   best_model <- sample(mod_eval[mod_eval$species==spID & 
                        !is.na(mod_eval$roc) &
                        mod_eval$roc==max(mod_eval[mod_eval$species==spID,"roc"], na.rm=T) & 
                        !is.na(mod_eval$model),"model"],1)
   best_mods <- rbind(best_mods, mod_eval[mod_eval$species==spID & mod_eval$model==best_model & !is.na(mod_eval$model),])
}

best_mods$best <- 1
best_mods

print("Now everything should be fine. Or?")
setdiff(unique(mod_eval$species), unique(best_mods$species))

# merge best_mods with mod_eval (add identified best model)
mod_eval <- mod_eval %>% full_join(best_mods[,c("species","model", "best")], by=c("model", "species"))

write.csv(mod_eval, file=paste0(here::here(), "/results/ModelEvaluation_", Taxon_name,".csv"), row.names = F)

