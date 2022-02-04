#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#        Colinearity of predictors          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

myExpl <- stack(paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate variable inflation factor (VIF) ####
# VIF is the extent of correlation between one predictor and all others.
# The lower VIF, the better we can tell what predictor contributed (most) to the model

## VIF basen on raw data (explanatory raster stack)
env.vif <- usdm::vif(myExpl)

# which predictors should be excluded?
vifcor(myExpl, th=0.8)  #th = threshold vif for exclusion

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate correlations between predictors ####
# https://github.com/joaofgoncalves/GoncalvesAna_et_al_2021/tree/master/RCODE/PostModelAnalyses
myExpl <- readRDS(file=paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd")) %>% as.data.frame

corMatSpearman <- cor(myExpl, method="spearman") %>% round(2)
corMatPearson <- cor(myExpl, method="pearson") %>% round(2)

#write.csv(corMatSpearman,"./OUT/corMatSpearman.csv")
#write.csv(corMatPearson,"./OUT/corMatPearson.csv")

caret::findCorrelation(corMatPearson, cutoff = 0.8, names=TRUE)
caret::findCorrelation(corMatPearson, cutoff = 0.7, names=TRUE)



