#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare input data set for SDM       #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

## needed for biomod again
# environmental (explanatory) variables as raster file
Env_norm <- raster::stack(paste0(here::here(), "/results/EnvPredictor_2km_normalized.grd"))
#Env_norm <- stack(Env_norm)

#- - - - - - - - - - - - - - - - - - - - - -
## define function to split data ####
split.data <- function(x){
  training <- x
  
  # split data into training (80%, including testing data), and validation data (20%)
  set.seed(1)
  random.rows <- sample(1:nrow(training), nrow(training))
  
  validation <- training[random.rows[1:round(0.2*nrow(training))], 
                             c("x","y", "SpeciesID", "occ", covarsNames[covarsNames %in% colnames(x)])]
  
  training <- training[random.rows[round(0.2*nrow(training)):nrow(training)],]
  
  # subset uncorrelated covariates
  training <- training[, c("occ", covarsNames[covarsNames %in% colnames(x)])]
  
  # # convert the categoricals to factor
  # training$vegsys <- as.factor(training$vegsys)
  # validation_env$vegsys <- as.factor(validation_env$vegsys)

  # # print the first few rows and columns
  # print(paste0("Head of the prepared dataset for ", modelName, "."))
  # print(training[1:5, 1:5])
  
  # save all datasets
  save(training, file=paste0(here::here(), "/results/", Taxon_name, "/TrainingData_", modelName, runningNumber, "_", Taxon_name,"_", temp_species, ".RData"))
  save(validation, file=paste0(here::here(), "/results/", Taxon_name, "/ValidationData_", modelName, runningNumber, "_", Taxon_name,"_", temp_species, ".RData"))
  }

#- - - - - - - - - - - - - - - - - - - - - -
## parallelize ####
# Calculate the number of cores
no.cores <- detectCores()/2; no.cores

# Initiate cluster used in foreach function
registerDoParallel(no.cores)

# check if right data are loaded
if(length(speciesNames$NumCells_2km) == 0){
  print("Error: Please make sure to load the Species List containing the column NumCells (from results folder)!")
}

# for loop
foreach(temp_species=speciesNames[speciesNames$NumCells_2km >=5,]$SpeciesID, 
        .export = c("Taxon_name", "covarsNames"), 
        .packages = c("tidyverse")) %dopar% {
  
  try({
    
    # load background data (pseudo-absences) for each modeling approach
    load(file=paste0(here::here(), "/results/", Taxon_name, "/PA_Env_", Taxon_name, "_", temp_species, ".RData")) #bg.list
    
    for(i in 1:(length(bg.list)-1)){ #exclude biomod data
      
      # for models with multiple background data runs
      if(ncol(bg.list[[i]])>(length(names(Env_norm))+4)){
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
          bg.temp.list[[j]][bg.temp.list[[j]][,1]==TRUE | is.na(bg.temp.list[[j]][,1]),]$occ <- 1
        }
        bg.list[[i]] <- bg.temp.list
        
      # for models with single background runs  
      }else{
        # rename column with presence-absence information
        bg.list[[i]] <- rename(bg.list[[i]], "PA1"=colnames(bg.list[[i]])[stringr::str_detect(colnames(bg.list[[i]]), "data.species")])
        
        # define presence and absence as new column
        bg.list[[i]]$occ <- 1
        bg.list[[i]][is.na(bg.list[[i]]$PA1),]$occ <- 0
      }
    }  
    
    ## split all datasets into training, testing and validation data
    # but first, copy bg.glm to choose to normalize/not normalize 
    bg.list[["bg.maxent"]] <- bg.list[["bg.glm"]]
    bg.list[["bg.lasso"]] <- bg.list[["bg.glm"]]
    
    # now we run the function for all data
    for(k in 1:(length(bg.list))){
      modelName <- names(bg.list)[k]
      runningNumber <- 1
      
      # GLM/ GAM
      if(modelName == "bg.glm"){ 
        split.data(bg.list[[k]])
      }
      
      # MaxEnt or Lasso/ Ridge     
      if(modelName == "bg.maxent" | modelName == "bg.lasso"){
        split.data(bg.list[[k]])
      }
      
      # MARS/ RF/ BRT/ SVM
      if(modelName == "bg.mars" | modelName == "bg.rf" | modelName == "bg.mda"){
          # if we have only a dataframe in the current list element, its easy
          if(class(bg.list[[k]])!="list"){
            runningNumber <- 1
            split.data(bg.list[[k]])
        
          # if we have a list of dataframes, we have to specify the right list layer
          }else{ for(l in 1:length(bg.list[[k]])){
              runningNumber <- l
              split.data(bg.list[[k]][[l]])
          }}
      }
      
    ## For BIOMOD ####
    if(modelName == "bg.biomod"){
      
      runningNumber <- 1
      
      ## re-create background data with bg.glm background data as evaluation data
      # load testing background data
      load(file=paste0(here::here(), "/results/", Taxon_name, "/ValidationData_bg.glm1_", Taxon_name,"_", temp_species, ".RData"))  #validation
      
      # tmp <- Sys.time()
      # bg.biomod <- biomod2::BIOMOD_FormatingData(resp.var = bg.list[["bg.biomod"]][["myResp"]],
      #                                            expl.var = myExpl,
      #                                            resp.xy = bg.list[["bg.biomod"]][["myRespCoord"]],
      #                                            resp.name = bg.list[["bg.biomod"]][["myRespName"]],
      #                                            PA.nb.rep = bg.list[["bg.biomod"]][["temp.runs"]],
      #                                            PA.nb.absences = bg.list[["bg.biomod"]][["temp.number"]],
      #                                            PA.strategy = bg.list[["bg.biomod"]][["temp.strategy"]],
      #                                            #eval.expl.var = validation_env,
      #                                            #eval.resp.var = validation_pa[,"occ"], 
      #                                            #eval.resp.xy = validation_pa[,c("x", "y")],
      #                                            )
      # temp.time <- Sys.time() - tmp
      # 
      
      training <- bg.list[["bg.biomod"]]
      
      ## take validation data from bg.glm 
      # NOTE: validation only needed if no evaluation data provided above (eval.expl etc. ...)!
      
      # save datasets for BIOMOD
      save(training, file=paste0(here::here(), "/results/", Taxon_name, "/TrainingData_", modelName, runningNumber, "_", Taxon_name,"_", temp_species, ".RData"))
      save(validation, file=paste0(here::here(), "/results/", Taxon_name, "/ValidationData_", modelName,runningNumber, "_", Taxon_name,"_", temp_species, ".RData"))
     }   
    }
  }) # end of try loop
  
  print(paste0("The model input data is now saved for ", Taxon_name, ": ", temp_species, "."))
  
}

# at the end
stopImplicitCluster()
