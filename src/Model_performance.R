#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

load(file=paste0(here::here(), "/sdm/SDM_Models_", spID, ".RData"))

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
                       time = NA, # mean computational time
                       bg= NA, # background dataset used
                       no.runs = NA) #how many background data runs to we have

# function to calculate TSS
opt.cut <- function(perf, pred){
  cut.ind <- mapply(FUN=function(x, y, p){
    d <- (x - 0)^2 + (y-1)^2
    ind <- which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}


for(i in 1:length(SDMs)){ 
  temp.model <- names(SDMs)[[i]]
  number.models <- 1
  
  #- - - - - - - - - - - - - - - - - - - 
  ## calculate statistics for BIOMOD ####
  if(temp.model=="biomod_pred") {
    myBiomodModelEval <- SDMs[["biomod_pred"]][[1]]
    
    # define background dataset (for testing data)
    modelName <- SDMs[[i]][[3]]
    
    # identify and load all relevant data files
    temp.files <- list.files(path = paste0(here::here(),"/results/",Taxon_name), 
                             pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
    lapply(temp.files, load, .GlobalEnv)
    
    precrec_obj <- precrec::auc(precrec::evalmod(scores = SDMs[["biomod_pred"]][[2]], 
                                                     labels = validation_pa$occ))
    temp.prg <- prg::calc_auprg(prg::create_prg_curve(labels = validation_pa$occ, 
                                                             pos_scores = SDMs[["biomod_pred"]][[2]]))
    
    prg_curve <- create_prg_curve(labels = validation_pa[,"occ"], pos_scores = prediction)
      
    
    ## TSS and Kappa
    roc.pred <- prediction(prediction[names(prediction) %in% rownames(validation_pa)], validation_pa[rownames(validation_pa) %in% names(prediction),"occ"])
    roc.perf <- ROCR::performance(roc.pred, measure = "tpr", x.measure = "fpr")
    
    # calculate Sensitivity and Specificity
    sen_spe <- opt.cut(roc.perf, roc.pred)
    
    # calculate TSS
    temp.tss <- as.numeric(sen_spe [1,1] +  sen_spe[2,1] - 1 )
    
    # calculate Kappa
    mod.object <- sdm::evaluates(x = validation_pa[,"occ"], p = prediction)
    temp.kappa <- mod.object@threshold_based[mod.object@threshold_based$criteria=="max(se+sp)", "Kappa"] #kappa at max(se+sp)
    temp.tresh <- mod.object@threshold_based[mod.object@threshold_based$criteria=="max(se+sp)", "threshold"] #threshold max(se+sp)
    
    # save in output data frame
    mod_eval[i,]$species <- spID
    try(mod_eval[i,]$roc <- precrec_obj[precrec_obj$curvetypes=="ROC", "aucs"])
    #try(mod_eval[i,]$roc <- myBiomodModelEval["ROC", "Testing.data"])
    try(mod_eval[i,]$prg <- temp.prg)
    #try(mod_eval[i,]$prg <- NA)
    try(mod_eval[i,]$cor <- cor(SDMs[["biomod_pred"]][[2]], validation_pa$occ))
    #try(mod_eval[i,]$cor <- NA)
    #try(mod_eval[i,]$tss <-  myBiomodModelEval["TSS", "Testing.data"])
    try(mod_eval[i,]$tss <- temp.tss)
    #try(mod_eval[i,]$kappa <-  myBiomodModelEval["KAPPA", "Testing.data"])
    try(mod_eval[i,]$kappa <- temp.kappa)
    
    # threshold at max TSS (Note: biomod predicion not scaled, their range is 0-1000
    try(mod_eval[i,]$thres.maxTSS <-  myBiomodModelEval["TSS", "Cutoff"] / 1000) 
    try(mod_eval[i,]$thres.maxTSS <-  temp.tresh) # threshold at max(se+sp)
     
    mod_eval[i,]$model <- temp.model
    try(mod_eval[i,]$time <- SDMs[[i]][[4]])
    mod_eval[i,]$bg <- modelName
    mod_eval[i,]$no.runs <- number.models
    
    rm(prg, temp.model, modelName, precrec_obj)
    
  }else{  # all except BIOMOD modeling
    
    # if necessary, unlist models
    if(length(SDMs[[i]])>6) {
      number.models <- length(SDMs[[i]])
    }
    
    # define background dataset (for testing data)
    modelName <- SDMs[[i]][[3]]
    
    # identify and load all relevant data files
    temp.files <- list.files(path = paste0(here::here(),"/results/",Taxon_name), 
                             pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
    lapply(temp.files, load, .GlobalEnv)
    
    # load prediction
    prediction <- SDMs[[temp.model]][[2]]
    
    # try loop (not stopping if error)
    try({
      ## ROC and PR
      # calculate area under the ROC and PR curves
      precrec_obj <- evalmod(scores = prediction, labels = validation_pa[,"occ"])
      #print(precrec_obj)
      
      # # plot the ROC and PR curves
      # roc <- autoplot(precrec_obj, curvetype = "ROC")
      # pr <- autoplot(precrec_obj, curvetype = "PR")
      # 
      ## PRG
      # calculate the PRG curve
      prg_curve <- create_prg_curve(labels = validation_pa[,"occ"], pos_scores = prediction)
      
      # # calculate area under the PRG cure (AUC PRG)
      # au_prg <- calc_auprg(prg_curve)
      # #print(au_prg)
      # 
      # # plot the PRG curve
      # prg <- plot_prg(prg_curve)
      
      # # plot all plots into one
      # ggpubr::ggarrange(roc, pr, prg,
      #                   labels=c("", "", "PRG", ""), nrow = 2, ncol=2) # second row
    })
    
    #- - - - - - - - - - - - - - - - - - - 
    # scale predictions between 0 and 1 
    prediction <- scales::rescale(prediction, to = c(0,1))
    
    ## TSS and Kappa
    roc.pred <- prediction(prediction[names(prediction) %in% rownames(validation_pa)], validation_pa[rownames(validation_pa) %in% names(prediction),"occ"])
    roc.perf <- ROCR::performance(roc.pred, measure = "tpr", x.measure = "fpr")
    
    # calculate Sensitivity and Specificity
    sen_spe <- opt.cut(roc.perf, roc.pred)
    
    # calculate TSS
    temp.tss <- as.numeric(sen_spe [1,1] +  sen_spe[2,1] - 1 )
    
    # calculate Kappa
    mod.object <- sdm::evaluates(x = validation_pa[,"occ"], p = prediction)
    temp.kappa <- mod.object@threshold_based[mod.object@threshold_based$criteria=="max(se+sp)", "Kappa"] #kappa at max(se+sp)
    temp.tresh <- mod.object@threshold_based[mod.object@threshold_based$criteria=="max(se+sp)", "threshold"] #threshold max(se+sp)
    
    # save in output data frame
    mod_eval[i,]$species <- spID
    try(mod_eval[i,]$roc <- precrec::auc(precrec_obj)[1,4])
    try(mod_eval[i,]$prg <- prg::calc_auprg(prg_curve))
    try(mod_eval[i,]$cor <- cor(prediction,  validation_pa[,"occ"]))
    try(mod_eval[i,]$tss <- temp.tss)
    try(mod_eval[i,]$kappa <- temp.kappa)
    try(mod_eval[i,]$thres.maxTSS <-  temp.tresh) # threshold at max(se+sp)
    mod_eval[i,]$model <- temp.model
    try(mod_eval[i,]$time <- SDMs[[i]][[4]])
    mod_eval[i,]$bg <- modelName
    mod_eval[i,]$no.runs <- number.models
  }
  
  rm(temp.tss, sen_spe, temp.model, modelName, precrec_obj, temp.kappa, temp.tresh)
}

mod_eval

#- - - - - - - - - - - - - - - - - - - 
## Save ####

write.csv(mod_eval, file=paste0(here::here(), "/results/ModelEvaluation_", Taxon_name, "_", spID, ".csv"), row.names = F)




