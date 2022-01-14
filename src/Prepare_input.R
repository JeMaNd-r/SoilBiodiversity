#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare input data set for SDM       #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load background data (pseudo-absences) for each modeling approach
load(file=paste0(here::here(), "/results/", Taxon_name, "/PA_Env_", Taxon_name, "_", spID, ".RData"))

for(i in 1:(length(bg.list)-1)){ #exclude biomod data
  if(ncol(bg.list[[i]])>(length(covarsNames)+4)){
    # define presence and absence as new colum
    bg.list[[i]]$occ <- 1
    
    # split the datasets into the runs
    bg.temp.list <- vector(mode="list", 
                           length=sum(str_detect(colnames(bg.list[[i]]), "PA")))
    for(j in 1:sum(str_detect(colnames(bg.list[[i]]), "PA"))){
      temp.columns <- c(paste0("PA", j),c("x", "y", "SpeciesID", "occ"), 
                        covarsNames[covarsNames %in% colnames(bg.list[[i]])])
      bg.temp.list[[j]] <- bg.list[[i]][,temp.columns]
      
      bg.temp.list[[j]]$occ <- 0
      bg.temp.list[[j]][bg.temp.list[[j]][,1]==TRUE,]$occ <- 1
    }
    bg.list[[i]] <- bg.temp.list
  }else{
    # rename column with presence-absence information
    bg.list[[i]] <- rename(bg.list[[i]], "PA1"=colnames(bg.list[[i]])[str_detect(colnames(bg.list[[i]]), "data.species")])
    
    # define presence and absence as new column
    bg.list[[i]]$occ <- 1
    bg.list[[i]][is.na(bg.list[[i]]$PA1),]$occ <- 0
  }
}  

## split all datasets into training, testing and validation data
# define function to do so
split.data <- function(x){
  training <- x
  
  # split data into training (60%), testing (20%) and validation data (20%)
  set.seed(1)
  random.rows <- sample(1:nrow(training), nrow(training))
  
  testing_env <- training[random.rows[1:round(0.2*nrow(training))], 
                          covarsNames[covarsNames %in% colnames(x)]]
  testing_pa <- training[random.rows[1:round(0.2*nrow(training))], 
                         c("x","y", "SpeciesID", "occ")]
  
  validation.data <- training[random.rows[round(0.2*nrow(training)):round(0.4*nrow(training))],]
  
  training <- training[random.rows[round(0.4*nrow(training)):nrow(training)],]
  
  # subset uncorrelated covariates
  training <- training[, c("occ", covarsNames[covarsNames %in% colnames(x)])]
  testing_env <- testing_env[, covarsNames[covarsNames %in% colnames(x)]]
  
  # # convert the categoricals to factor
  # training$vegsys <- as.factor(training$vegsys)
  # testing_env$vegsys <- as.factor(testing_env$vegsys)
  
  # normalize the covariates (exept categorical)
  # *notice: not all the models are fitted on normalized data in
  # the main analysis! Please check Valavi et al. 2021.
  for(v in covarsNames[covarsNames %in% colnames(x)]){
    meanv <- mean(training[,v], na.rm=T)
    sdv <- sd(training[,v], na.rm=T)
    training[,v] <- (training[,v] - meanv) / sdv
    testing_env[,v] <- (testing_env[,v] - meanv) / sdv
  }
  
  # # print the first few rows and columns
  # print(paste0("Head of the prepared dataset for ", modelName, "."))
  # print(training[1:5, 1:5])
  
  # save all datasets
  save(training, file=paste0(here::here(), "/results/", Taxon_name, "/TrainingData_", modelName, runningNumber, "_", Taxon_name,"_", spID, ".RData"))
  save(testing_pa, file=paste0(here::here(), "/results/", Taxon_name, "/TestingData_pa_", modelName,runningNumber, "_", Taxon_name,"_", spID, ".RData"))
  save(testing_env, file=paste0(here::here(), "/results/", Taxon_name, "/TestingData_env_", modelName,runningNumber, "_", Taxon_name,"_", spID, ".RData"))
  save(validation.data, file=paste0(here::here(), "/results/", Taxon_name, "/ValidationData_", modelName,runningNumber, "_", Taxon_name,"_", spID, ".RData"))
}

# now we run the function for all data
for(k in 1:(length(bg.list)-1)){ #exclude biomod data
  modelName <- names(bg.list)[k]
  
  # if we have only a dataframe in the current list element, its easy
  if(class(bg.list[[k]])!="list"){
    runningNumber <- 1
    split.data(bg.list[[k]])
    
    # if we have a list of dataframes, we have to specify the right list layer
  }else{
    for(l in 1:length(bg.list[[k]])){
      runningNumber <- l
      split.data(bg.list[[k]][[l]])
    }
  }
}

## For BIOMOD ####
modelName <- "bg.biomod"
runningNumber <- 1

training <- bg.list$bg.biomod

# take testing anjd validation data from bg.glm
load(file=paste0(here::here(), "/results/", Taxon_name, "/TestingData_pa_", "bg.glm",runningNumber, "_", Taxon_name,"_", spID, ".RData"))
load(file=paste0(here::here(), "/results/", Taxon_name, "/TestingData_env_", "bg.glm",runningNumber, "_", Taxon_name,"_", spID, ".RData"))
load(file=paste0(here::here(), "/results/", Taxon_name, "/ValidationData_", "bg.glm",runningNumber, "_", Taxon_name,"_", spID, ".RData"))

# save datasets for BIOMOD
save(training, file=paste0(here::here(), "/results/", Taxon_name, "/TrainingData_", modelName, runningNumber, "_", Taxon_name,"_", spID, ".RData"))
save(testing_pa, file=paste0(here::here(), "/results/", Taxon_name, "/TestingData_pa_", modelName,runningNumber, "_", Taxon_name,"_", spID, ".RData"))
save(testing_env, file=paste0(here::here(), "/results/", Taxon_name, "/TestingData_env_", modelName,runningNumber, "_", Taxon_name,"_", spID, ".RData"))
save(validation.data, file=paste0(here::here(), "/results/", Taxon_name, "/ValidationData_", modelName,runningNumber, "_", Taxon_name,"_", spID, ".RData"))



print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
print(paste0("The model input data is now saved for ", Taxon_name, ": ", spID, "."))
print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
