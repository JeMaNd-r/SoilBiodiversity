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

mod_eval <- data.frame(species = NA, # name of the species (more important for later)
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

for(i in 1:length(sdm_names)){ 
  temp.model <- sdm_names[i]
  number.models <- 1
 
  print("=====================================")
  print(temp.model)

  try(temp_sdm <- get(load(file=paste0(here::here(), "/results/", Taxon_name, "/temp_files/SDM_",temp.model,"_", spID, ".RData"))), silent=T)
  
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

    # doesn't work because there is no inout validation data... 
    # you may have to run biomod2 twice, one time with full env. space and one time with validation dataset...
    try({
      precrec_obj <- precrec::auc(precrec::evalmod(scores = prediction$layer, 
                                                     labels = validation[,c("x","y","SpeciesID", "occ")]$occ))
      temp.prg <- prg::calc_auprg(prg::create_prg_curve(labels = validation[,c("x","y","SpeciesID", "occ")]$occ, 
                                                               pos_scores = prediction$layer))
      
      prg_curve <- create_prg_curve(labels = validation[,"occ"], pos_scores = prediction$layer)
        
      
      ## TSS and Kappa
      roc.pred <- prediction(prediction[names(prediction) %in% rownames(validation)], validation[rownames(validation) %in% names(prediction),"occ"])
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
    try(mod_eval[i,]$roc <- precrec_obj[precrec_obj$curvetypes=="ROC", "aucs"], silent=T)
    #try(mod_eval[i,]$roc <- myBiomodModelEval["ROC", "Testing.data"])
    try(mod_eval[i,]$prg <- temp.prg, silent=T)
    #try(mod_eval[i,]$prg <- NA)
    try(mod_eval[i,]$cor <- cor(prediction, validation$occ), silent=T)
    #try(mod_eval[i,]$cor <- NA)
    #try(mod_eval[i,]$tss <-  myBiomodModelEval["TSS", "Testing.data"])
    try(mod_eval[i,]$tss <- temp.tss, silent=T)
    #try(mod_eval[i,]$kappa <-  myBiomodModelEval["KAPPA", "Testing.data"])
    try(mod_eval[i,]$kappa <- temp.kappa, silent=T)
    
    # threshold at max TSS (Note: biomod predicion not scaled, their range is 0-1000
    try(mod_eval[i,]$thres.maxTSS <-  myBiomodModelEval[2,2] / 1000, silent=T) 
    try(mod_eval[i,]$thres.maxTSS <-  temp.tresh, silent=T) # threshold at max(se+sp)
     
    mod_eval[i,]$model <- temp.model
    try(mod_eval[i,]$time_model <- paste0(temp_sdm[["time_model"]][1], " ", temp_sdm[["time_model"]][2]), silent=T)
    try(mod_eval[i,]$time_predict <- paste0(temp_sdm[["time_predict"]][1], " ", temp_sdm[["time_predict"]][2]), silent=T)
    mod_eval[i,]$bg <- modelName
    mod_eval[i,]$no.runs <- 1
    try(mod_eval[i,]$n_vali_presence <- sum(validation$occ), silent=T)
    try(mod_eval[i,]$n_vali_absence <- (length(validation$occ) - sum(validation$occ)), silent=T)
    
    if(is.na(mod_eval[i,]$tss)==T) mod_eval[i,]$tss <- myBiomodModelEval[2,1]
    if(is.na(mod_eval[i,]$roc)==T) mod_eval[i,]$roc <- myBiomodModelEval[3,1]
    if(is.na(mod_eval[i,]$kappa)==T) mod_eval[i,]$kappa <- myBiomodModelEval[1,1]
    
    rm(prg, temp.model, modelName, precrec_obj)
    
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
    
      if(sum(validation$occ)==0) print(paste0("For model ", temp.model, ": There is no presence records in the validation data. Too few occurrence records."))

      # try loop (not stopping if error)
      try({
        ## ROC and PR
        # calculate area under the ROC and PR curves
        precrec_obj <- evalmod(scores = prediction, labels = validation[,"occ"])      

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
    mod_eval[i,]$model <- temp.model
    try(mod_eval[i,]$time_model <- paste0(temp_sdm[["time_model"]][1], " ", temp_sdm[["time_model"]][2]), silent=T)
    try(mod_eval[i,]$time_predict <- paste0(temp_sdm[["time_predict"]][1], " ", temp_sdm[["time_predict"]][2]), silent=T)
    mod_eval[i,]$bg <- modelName
    mod_eval[i,]$no.runs <- temp_sdm[["runs"]]
    try(mod_eval[i,]$n_vali_presence <- sum(validation$occ), silent=T)
    try(mod_eval[i,]$n_vali_absence <- (length(validation$occ) - sum(validation$occ)), silent=T)
  }
  
  rm(temp.tss, sen_spe, temp.model, modelName, precrec_obj, temp.kappa, temp.tresh, temp_sdm)
}

mod_eval

#- - - - - - - - - - - - - - - - - - - 
## Save ####

write.csv(mod_eval, file=paste0(here::here(), "/results/", Taxon_name,"/ModelEvaluation_", Taxon_name, "_", spID, ".csv"), row.names = F)

