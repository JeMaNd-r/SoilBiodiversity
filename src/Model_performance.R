#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(precrec)
library(ggplot2) # for plotting the curves
# devtools::install_github("meeliskull/prg/R_package/prg")
library(prg)
library(ggpubr)

load(file=paste0(here::here(), "/results/", Taxon_name, "/Predicted_SDMs_", spID, ".RData"))

#- - - - - - - - - - - - - - - - - - - - -
## Model performance ####
#- - - - - - - - - - - - - - - - - - - - -

for(i in length(names(SDMs))){
  temp.model <- names(SDMs)[[i]]
  
  # load relevant files
  modelName <- SDMs[[i]][[2]]
  
  # identify and load all relevant data files
  temp.files <- list.files(path = paste0("./results/",Taxon_name), 
                           pattern = paste0(modelName, "[[:graph:]]*", spID), full.name = T)
  lapply(temp.files, load, .GlobalEnv)
  
  ## ROC and PR
  # load prediction
  prediction <- SDMs[[temp.model]][[1]]
  
  # calculate area under the ROC and PR curves
  precrec_obj <- evalmod(scores = prediction, labels = testing_pa[,"occ"])
  print(precrec_obj)
  
  # plot the ROC and PR curves
  roc <- autoplot(precrec_obj, curvetype = "ROC")
  pr <- autoplot(precrec_obj, curvetype = "PR")
  
  ## PRG
  # calculate the PRG curve for RF down-sampled
  prg_curve <- create_prg_curve(labels = testing_pa[,"occ"], pos_scores = prediction)
  
  # calculate area under the PRG cure
  au_prg <- calc_auprg(prg_curve)
  print(au_prg)
  
  # plot the PRG curve
  prg <- plot_prg(prg_curve)
  
  # plot all plots into one
  # ####TODO: ADD precrec obj to plot as text...?
  ggpubr::ggarrange(roc, pr, prg, plot(label=print(precrec_obj)),
                    labels=c("", "", "PRG", ""), nrow = 2, ncol=2) # second row
                    
}

