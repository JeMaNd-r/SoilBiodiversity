#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

load(file=paste0(here::here(), "/results/", Taxon_name, "/Predicted_SDMs_", spID, ".RData"))

#- - - - - - - - - - - - - - - - - - - - -
## Model performance ####
#- - - - - - - - - - - - - - - - - - - - -

mod_eval <- data.frame(species = NA, 
                       roc = NA,
                       prg = NA,
                       cor = NA,
                       model = NA,
                       time = NA,
                       bg= NA)

for(i in 1:length(names(SDMs))){
  temp.model <- names(SDMs)[[i]]
  
  # # if necessary, unlist models
  # if(length(SDMs[[i]])!=3) {
  #   
  # }
    
  # define background dataset (for testing data)
  modelName <- SDMs[[i]][[3]]
  
  # identify and load all relevant data files
  temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                           pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
  lapply(temp.files, load, .GlobalEnv)
  
  ## ROC and PR
  # load prediction
  prediction <- SDMs[[temp.model]][[2]]
  
  # calculate area under the ROC and PR curves
  precrec_obj <- evalmod(scores = prediction, labels = testing_pa[,"occ"])
  #print(precrec_obj)
  
  # plot the ROC and PR curves
  roc <- autoplot(precrec_obj, curvetype = "ROC")
  pr <- autoplot(precrec_obj, curvetype = "PR")
  
  ## PRG
  # calculate the PRG curve
  prg_curve <- create_prg_curve(labels = testing_pa[,"occ"], pos_scores = prediction)
  
  # calculate area under the PRG cure (AUC PRG)
  au_prg <- calc_auprg(prg_curve)
  #print(au_prg)
  
  # plot the PRG curve
  prg <- plot_prg(prg_curve)
  
  # plot all plots into one
  ggpubr::ggarrange(roc, pr, prg,
                    labels=c("", "", "PRG", ""), nrow = 2, ncol=2) # second row
  
  mod_eval[i,]$species <- spID
  try(mod_eval[i,]$roc <- precrec::auc(precrec_obj)[1,4])
  try(mod_eval[i,]$prg <- prg::calc_auprg(prg_curve))
  try(mod_eval[i,]$cor <- cor(prediction,  testing_pa[,"occ"]))
  mod_eval[i,]$model <- temp.model
  try(mod_eval[i,]$time <- SDMs[[i]][[4]])
  mod_eval[i,]$bg <- modelName
                    
}

write.csv(mod_eval, file=paste0(here::here(), "/results/ModelEvaluation_", Taxon_name, "_", spID, ".csv"), row.names = F)


## calculate additional indices according to equations in Allouche et al. 2006

# Sensitivity


# Specificity 



# TSS


# kappa





