#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare input data set for SDM       #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# load background data (pseudo-absences) for each modeling approach
bg.glm <- read.csv(file=paste0(here::here(), "/results/BackgroundData_GLM_", Taxon_name, ".csv"))
bg.mars <- read.csv(file=paste0(here::here(), "/results/BackgroundData_MARS_", Taxon_name, ".csv"))
bg.mda <- read.csv(file=paste0(here::here(), "/results/BackgroundData_MDA_", Taxon_name, ".csv"))
bg.rf <- read.csv(file=paste0(here::here(), "/results/BackgroundData_RF_", Taxon_name, ".csv"))

# define presence and absence as new column
bg.glm$occ <- 1
bg.glm[is.na(bg.glm$bg.glm.data.species),]$occ <- 0

bg.rf$occ <- 1
bg.rf[is.na(bg.rf$bg.rf.data.species),]$occ <- 0

# do the same for the other data sets, but this time, split the datasets into the runs
bg.mars$occ <- 0
bg.mars.list <- vector(mode="list", length=sum(str_detect(colnames(bg.mars), "PA")))
for(i in 1:sum(str_detect(colnames(bg.mars), "PA"))){
  temp.columns <- c(paste0("PA", i),c("x", "y", "SpeciesID", "occ"), covarsNames[covarsNames %in% colnames(bg.mars)])
  bg.mars.list[[i]] <- bg.mars[,temp.columns]
  
  bg.mars.list[[i]][bg.mars[,1]==TRUE,]$occ <- 1
}
  
bg.mda$occ <- 0
bg.mda.list <- vector(mode="list", length=sum(str_detect(colnames(bg.mda), "PA")))
for(i in 1:sum(str_detect(colnames(bg.mda), "PA"))){
  temp.columns <- c(paste0("PA", i),c("x", "y", "SpeciesID", "occ"), covarsNames[covarsNames %in% colnames(bg.mda)])
  bg.mda.list[[i]] <- bg.mda[,temp.columns]
  
  bg.mda.list[[i]][bg.mda[,1]==TRUE,]$occ <- 1
}

# combine the presence and background points
training.glm <- rbind(pr, bg.glm)
training.rf <- rbind(pr, bg.rf)
training.mars1 <- rbind(pr, bg.mars)
training.mars2 <- rbind(pr, bg.glm)
training.mars3 <- rbind(pr, bg.glm)
training.mars1 <- rbind(pr, bg.glm)
training.mars1 <- rbind(pr, bg.glm)
training.mars1 <- rbind(pr, bg.glm)
training.mars1 <- rbind(pr, bg.glm)
training.mars1 <- rbind(pr, bg.glm)

## do for all background datasets:
for(j in c("bg.glm", "bg.mda", "bg.mars", "bg.rf")){

  training <- Values(paste0("training.", j))
  
  # split data into training (60%), testing (20%) and validation data (20%)
  set.seed(1)
  random.rows <- sample(1:nrow(training), nrow(training))
  
  testing_env <- training[random.rows[1:round(0.2*nrow(training))], predictorNames]
  testing_pa <- training[random.rows[1:round(0.2*nrow(training))], c("x","y", "SpeciesID", "occ")]
  
  validation.data <- training[random.rows[round(0.2*nrow(training)):round(0.4*nrow(training))],]
  
  training <- training[random.rows[round(0.4*nrow(training)):nrow(training)],]
  
  # subset uncorrelated covariates
  training <- training[, c("occ", covarsNames)]
  testing_env <- testing_env[, covarsNames]
  
  # # convert the categoricals to factor
  # training$vegsys <- as.factor(training$vegsys)
  # testing_env$vegsys <- as.factor(testing_env$vegsys)
  
  # normalize the covariates (exept categorical)
  # *notice: not all the models are fitted on normalized data in
  # the main analysis! Please check Valavi et al. 2021.
  for(v in covarsNames){
    meanv <- mean(training[,v], na.rm=T)
    sdv <- sd(training[,v], na.rm=T)
    training[,v] <- (training[,v] - meanv) / sdv
    testing_env[,v] <- (testing_env[,v] - meanv) / sdv
  }
  
  # print the first few rows and columns
  print("Head of the prepared dataset.")
  training[1:5, 1:5]
  
  # save all datasets
  write.csv(training, file=paste0(here::here(), "/results/TrainingData_", j, "_", Taxon_name, ".csv"), row.names = F)
  write.csv(testing_pa, file=paste0(here::here(), "/results/TestingData_pa_", j, "_", Taxon_name, ".csv"), row.names = F)
  write.csv(testing_env, file=paste0(here::here(), "/results/TestingData_env_", j, "_", Taxon_name, ".csv"), row.names = F)
  write.csv(validation.data, file=paste0(here::here(), "/results/ValidationData_", j, "_", Taxon_name, ".csv"), row.names = F)
}


