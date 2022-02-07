#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Variable importance              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

## load models ####
load(file=paste0(here::here(), "/results/SDM_Models_", spID, ".RData")) #SDMs

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate variable importance (VI) ####
# create result data frame
var_imp <- data.frame("Predictor"= predictorNames)

# for loop through all models
for(i in 1:length(SDMs)){ try({
  temp.model <- names(SDMs)[[i]]
  number.models <- 1
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## GLM and GAM ####
  if(temp.model == "gm_pred" | temp.model =="lm1_pred" | temp.model =="lm_subset_pred"){
    # calculate VI
    temp.vi <- caret::varImp(SDMs[[temp.model]][[1]], scale=T) #scaled between 0 and 100% 
  
    # check if it looks right
    #print(temp.vi)
    
    # replace Inf values with 0
    temp.vi[is.infinite(temp.vi[,1]),1] <- 0
      
    temp.vi$Predictor <- rownames(temp.vi)
    colnames(temp.vi)[1] <- temp.model
  }
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Lasso and rigde regression ####
  # because of the use of cv.glmnet function, we need to calculate varImp differently
  if(temp.model == "lasso_pred" | temp.model =="ridge_pred"){
    coefs <- coef(SDMs[[temp.model]][[1]]) #!!! simple way...
    
    # remove first element (Intercept)
    temp.vi <- data.frame("Predictor" = coefs@Dimnames[[1]][coefs@i + 1], temp.model = coefs@x)[-1,] 
    colnames(temp.vi)[2] <- temp.model
    
    temp.vi$Predictor <- stringr::str_split_fixed(temp.vi$Predictor, "_", 2)[,1]
    
    temp.vi <- temp.vi %>% group_by(Predictor) %>% summarize_all(mean, na.rm=T) %>% as.data.frame()
    
    # check if it looks right
    #print(temp.vi)
  }
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## MARS ####
  if(temp.model == "mars_pred"){# if necessary, unlist models
  
   temp.vi <- SDMs[[temp.model]][[5]]
   temp.vi <- as.data.frame(temp.vi) %>% dplyr::rename("mars_pred" = "Overall")
   #print(temp.vi)
  }
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ##  ####    
  if(temp.model == ""){   
   # calculate VI
    temp.vi <- caret::varImp(SDMs[[temp.model]][[1]][,j], scale=T,) #scaled between 0 and 100% 
    
    # check if it looks right
    print(temp.vi)
    
    temp.vi$Predictor <- rownames(temp.vi)
    colnames(temp.vi)[1] <- temp.model
    
    # replace NA with 0
    temp.vi[,2] <- temp.vi[,2] %>% tidyr::replace_na(0)
    
    temp.vi.df <- dplyr::full_join(temp.vi.df, temp.vi)
    
    rm(temp.vi) 
  }
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## MaxEnt ####
  if(temp.model == "maxmod_pred"){
    
    temp.results <- SDMs[[temp.model]][[1]]@results
    temp.vi <- as.data.frame(temp.results[str_detect(rownames(temp.results),"permutation.importance"),])
    
    # extract predictor names
    temp.vi$Predictor <- stringr::str_split_fixed(rownames(temp.vi), "[:punct:]", 2)[,1]
    colnames(temp.vi)[1] <- temp.model
    
    # plot variable importance
    #plot(SDMs[[temp.model]][[1]])
    
    # plot response curve
    #response(SDMs[[temp.model]][[1]])
    
    rm(temp.results)
  }    
  
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## MaxNet ####
  if(temp.model == "maxnet_pred"){
      
    temp.results <- SDMs[[temp.model]][[1]]
    
    # use package designed to calculate variable importance of maxnet object
    #fitMaxnet::varImportance(temp.results)
    
    # instead: extract beta coefficients
    temp.vi <- temp.results$betas
    temp.vi <- as.data.frame(temp.vi[str_detect(names(temp.vi),"hinge")]) 
    temp.vi$Predictor <- stringr::str_split_fixed(rownames(temp.vi), "[:punct:]", 3)[,2]
    
    temp.vi <- temp.vi %>% group_by(Predictor) %>% summarize_all(mean, na.rm=T) %>% as.data.frame()
    colnames(temp.vi) <- c("Predictor", temp.model)
    
    # plot response curves
    #plot(temp.results)
    
    rm(temp.results)
  }    
  
  # harmonizes predictor column structure
  temp.vi$Predictor <- as.character(temp.vi$Predictor)
  
  # scale variable importance values between 0 and 10
  temp.vi[,temp.model] <- scales::rescale(abs(temp.vi[,temp.model]), c(0,10))
  
  var_imp <- dplyr::full_join(var_imp, temp.vi)
  
  # replace NA with 0
  var_imp[,ncol(var_imp)][is.na(var_imp[,ncol(var_imp)])] <- 0
  
  rm(temp.vi)
})}

var_imp

#- - - - - - - - - - - - - - - - - - - - - -
## for BIOMOD ####

# https://github.com/joaofgoncalves/GoncalvesAna_et_al_2021/tree/master/RCODE/PostModelAnalyses

# Calculate variable importance across all PA sets, eval runs and algorithms 
varImportance <- biomod2::get_variables_importance(myBiomodModelOut)
#varImportanceByVariableAVG <- apply(varImportance,1,mean)
#varImportanceByVariableSTD <- apply(varImportance,1,sd)

varImportanceByVariableAVG <- apply(varImportance,1,mean, na.rm=TRUE)
varImportanceByVariableSTD <- apply(varImportance,1,sd, na.rm=TRUE)

vimpDF <- data.frame(cnames=names(varImportanceByVariableAVG),
                     vimpAVG = varImportanceByVariableAVG, 
                     varImpSTD=varImportanceByVariableSTD) %>% 
  arrange(desc(vimpAVG))

write.csv(vimpDF, file = paste(getwd(),"/",sp,"/",sp,"_varImportance.csv",sep=""))




#- - - - - - - - - - - - - - - - - - - - - -
## LASSO and RIDGE varImp - long and accurrate way... ####
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