#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Variable importance              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

## load models ####
load(file=paste0(here::here(), "/results/", Taxon_name, "/SDM_Models_", spID, ".RData")) #SDMs

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate variable importance (VI) ####
# create result data frame
var_imp <- data.frame("Predictor"= predictorNames)

# for loop through all models
for(i in 1:length(SDMs)){ 
  temp.model <- names(SDMs)[[i]]
  number.models <- 1
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Lasso and rigde regression ####
  # because of the use of cv.glmnet function, we need to calculate varImp differently
  if(temp.model == "lasso_pred" | temp.model =="ridge_pred"){
    coefs <- coef(SDMs[[temp.model]][[1]]) #!!! simple way...
    
    # remove first element (Intercept)
    temp.vi <- data.frame("Predictor" = coefs@Dimnames[[1]][coefs@i + 1], temp.model = coefs@x)[-1,]
    colnames(temp.vi)[2] <- temp.model
    
    temp.vi$Predictor <- stringr::str_split_fixed(temp.vi$Predictor, "_", 2)[,1]
    
    temp.vi[,2] <- scales::rescale(abs(temp.vi[,2]), c(0,10))
    
    # check if it looks right
    print(temp.vi)
  
  #- - - - - - - - - - - - - - - - - - - - - -
  }else{
    
    #- - - - - - - - - - - - - - - - - - - - - -
    ## MARS ####
    if(temp.model == "mars_pred"){# if necessary, unlist models
    
     temp.vi <- SDMs[[temp.model]][[5]]
     print(temp.vi)
      
    }else if(temp.model == ""){
          
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
      
    }else{
    
    # calculate VI
    temp.vi <- caret::varImp(SDMs[[temp.model]][[1]], scale=T) #scaled between 0 and 100% 
    
    # check if it looks right
    print(temp.vi)
    
    temp.vi$Predictor <- rownames(temp.vi)
    colnames(temp.vi)[1] <- temp.model
  }
  
  # replace NA with 0
  temp.vi[,2] <- temp.vi[,2] %>% tidyr::replace_na(0)
  
  var_imp <- dplyr::full_join(var_imp, temp.vi)
  
  rm(temp.vi)
}

var_imp

#- - - - - - - - - - - - - - - - - - - - - -
## for BIOMOD ####

# https://github.com/joaofgoncalves/GoncalvesAna_et_al_2021/tree/master/RCODE/PostModelAnalyses

# Calculate variable importance across all PA sets, eval runs and algorithms 
varImportance <- get_variables_importance(myBiomodModelOut)
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