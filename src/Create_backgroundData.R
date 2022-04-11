#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Background data/ Pseudo-absences     #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# according to recommendations of Barbet-Massin et al. (2012): large number 
# (e.g. 10,000) of pseudo-absences

# environmental (explanatory) variables as raster file
myExpl <- stack(paste0(here::here(), "/results/EnvPredictor_2km.grd"))

# response variable (i.e., species occurrences) in wide format
#mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_1km_", Taxon_name, ".csv"))
mySpeciesOcc <- read.csv(file=paste0(here::here(), "/results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))

## parallelize
# Calculate the number of cores
no.cores <- detectCores()/2; no.cores

# Initiate cluster used in foreach function
registerDoParallel(no.cores)

#- - - - - - - - - - - - - - - - - - - - - - - 
## For loop through all species ####
foreach(myRespName = speciesNames[speciesNames$NumCells >= 5,]$SpeciesID, .export = c("mySpeciesOcc"), 
        .packages = c("biomod2", "tidyverse")) %dopar% {
          
  # define response variable index
  myResp <- as.numeric(mySpeciesOcc[,myRespName])
 
  # get NAs id
  na.id <- which(is.na(myResp))
  
  # remove NAs to enforce PA sampling to be done on explanatory rasters
  myResp <- myResp[-na.id]
  myRespCoord <- mySpeciesOcc[-na.id,c('x','y')]
  
  # create summary table for model settings
  model.settings <- data.frame(SpeciesID=myRespName, model="X", strategy="test", 
                               No.runs=1, No.points=1, min.distance=1, run.time=0)[0,]
  
  ## GLM, GAM, Biomod ####
  # randomly performs consistently well, excepted when presences are climatically 
  #   biased for which ‘2°far’ is the best method
  # 10,000 PA or a minimum of 10 runs with 1,000 PA with an equal weight for 
  #   presences and absences
  temp.runs <- 1
  temp.number <- 10000
  temp.strategy <- "random"
  
  tmp <- Sys.time()
  bg.glm <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                 expl.var = myExpl,
                                 resp.xy = myRespCoord,
                                 resp.name = myRespName,
                                 PA.nb.rep = temp.runs,
                                 PA.nb.absences = temp.number,
                                 PA.strategy = temp.strategy)
  temp.time <- Sys.time() - tmp
  
  if(temp.runs==1){
    bg.glm <- cbind(bg.glm@data.species, bg.glm@coord, bg.glm@data.env.var)
  }else{
    bg.glm <- cbind(bg.glm@PA, bg.glm@coord, bg.glm@data.env.var)
  }
  bg.glm$SpeciesID <- myRespName
  
  str(bg.glm)
  #print("Note: TRUE in PA means presence.")
  
  # add model settings into summary table
  temp.dat <- data.frame(SpeciesID=myRespName, model="GLM.GAM", strategy=temp.strategy,
                         No.runs=temp.runs, No.points=temp.number, 
                         min.distance=NA, run.time=temp.time)
  model.settings <- rbind(model.settings, temp.dat)
  
  rm(temp.strategy, temp.runs, temp.number, temp.min.dist, temp.time)
  
  ## MARS ####	
  # random performs consistently well, except when presences are climatically 
  #    biased for which ‘2°far’ is the best method 
  # minimum of 10 runs with 100 PA
  temp.runs <- 10
  temp.number <- 100
  temp.strategy <- "random"
  
  tmp <- Sys.time()
  bg.mars <- BIOMOD_FormatingData(resp.var = myResp,
                                  expl.var = myExpl,
                                  resp.xy = myRespCoord,
                                  resp.name = myRespName,
                                  PA.nb.rep = temp.runs,
                                  PA.nb.absences = temp.number,
                                  PA.strategy = temp.strategy)
  temp.time <- Sys.time() - tmp
  
  bg.mars <- cbind(bg.mars@PA, bg.mars@coord, bg.mars@data.env.var)
  bg.mars$SpeciesID <- myRespName
  
  
  # add model settings into summary table
  temp.dat <- data.frame(SpeciesID=myRespName, model="MARS",strategy=temp.strategy, 
                         No.runs=temp.runs, No.points=temp.number, 
                         min.distance=NA, run.time=temp.time)
  model.settings <- rbind(model.settings, temp.dat)
  
  rm(temp.strategy, temp.runs, temp.number, temp.min.dist, temp.time)
  
  # ## MDA ####	
  # # ‘2°far’ performs consistently better with few presences, ‘SRE’ performs 
  # #   better with a large number of presences; 
  # #   randomly performs consistently well with spatially biased presences 	
  # # minimum of 10 runs with 100 PA with an equal weight for presences and absences
  # if(length(myResp) < 500){
  #   temp.runs <- 10
  #   temp.number <- 500
  #   temp.strategy <- "disk"
  #   temp.min.dist <- 222222
  # 
  #   }else{
  #     temp.runs <- 10
  #     temp.number <- 500
  #     temp.strategy <- "sre"
  #     temp.min.dist <- NULL
  # }
  # 
  # tmp <- Sys.time()
  # bg.mda <- BIOMOD_FormatingData(resp.var = myResp,
  #                                 expl.var = myExpl,
  #                                 resp.xy = myRespCoord,
  #                                 resp.name = myRespName,
  #                                 PA.nb.rep = temp.runs,
  #                                 PA.nb.absences = temp.number,
  #                                 PA.strategy = temp.strategy, 
  #                                PA.dist.min = temp.min.dist)
  # temp.time <- Sys.time() - tmp
  # 
  # bg.mda <- cbind(bg.mda@PA, bg.mda@coord, bg.mda@data.env.var)
  # bg.mda$SpeciesID <- myRespName
  # 
  # # add model settings into summary table
  # temp.dat <- data.frame(SpeciesID=myRespName, model="MDA", strategy=temp.strategy, 
  #                        No.runs=temp.runs, No.points=temp.number, 
  #                        min.distance=NA, run.time=temp.time)
  # model.settings <- rbind(model.settings, temp.dat)
  # 
  # rm(temp.strategy, temp.runs, temp.number, temp.min.dist, temp.time)
  
  ## CTA, BRT, RF ####	
  # ‘2°far’ performs consistently better with few presences, ‘SRE’ performs 
  #   better with a large number of presences 	
  # same as number of presences, 10 runs when less than 1000 PA with an equal 
  # weight for presences and absences
  if(length(myResp) < 500){
    temp.runs <- 10
    temp.number <- 100
    temp.strategy <- "disk"
    temp.min.dist <- 2000 #in meters
      
    }else{
      if(length(myResp) < 1000){
        temp.runs <- 10
        temp.number <- length(myResp)
        temp.strategy <- "sre"
        temp.min.dist <- NULL
      }else{
        temp.runs <- 1
        temp.number <- length(myResp)
        temp.strategy <- "sre"
        temp.min.dist <- NULL
    }
  }
  
  tmp <- Sys.time()
  bg.rf <- BIOMOD_FormatingData(resp.var = myResp,
                                  expl.var = myExpl,
                                  resp.xy = myRespCoord,
                                  resp.name = myRespName,
                                  PA.nb.rep = temp.runs,
                                  PA.nb.absences = temp.number,
                                  PA.strategy = temp.strategy, 
                                PA.dist.min = temp.min.dist)
  temp.time <- Sys.time() - tmp
  
  if(temp.runs==1){
    bg.rf <- cbind(bg.rf@data.species, bg.rf@coord, bg.rf@data.env.var)
    bg.rf$SpeciesID <- myRespName
  }else{
    bg.rf <- cbind(bg.rf@PA, bg.rf@coord, bg.rf@data.env.var)
    bg.rf$SpeciesID <- myRespName
  }
  
  # add model settings into summary table
  temp.dat <- data.frame(SpeciesID=myRespName, model="RF.BRT.CTA",strategy=temp.strategy, 
                         No.runs=temp.runs, No.points=temp.number, 
                         min.distance=NA, run.time=temp.time)
  model.settings <- rbind(model.settings, temp.dat)
  
  rm(temp.strategy, temp.runs, temp.number, temp.min.dist, temp.time)
  
  #- - - - - - - - - - - - - - - - - - - 
  ## BIOMOD ####
  temp.runs <- 1
  temp.number <- 10000 * 0.8 # take only 80% of the data as the other ones will
  # be splitted in "Prepare_input.r" in 80% training and 20% testing...
  temp.strategy <- "random"
  
  tmp <- Sys.time()
  bg.biomod <- biomod2::BIOMOD_FormatingData(resp.var = myResp,
                                             expl.var = myExpl,
                                             resp.xy = myRespCoord,
                                             resp.name = myRespName,
                                             PA.nb.rep = temp.runs,
                                             PA.nb.absences = temp.number,
                                             PA.strategy = temp.strategy)
  temp.time <- Sys.time() - tmp
  
  # NOTE: testing data and validation data will be used from bg.glm
  # therefore, we need to save some parameters for later
  bg.biomod <- list("bg.biomod"=bg.biomod, "myResp"=myResp, "myRespCoord"=myRespCoord, "myRespName"=myRespName,
                    "temp.runs"=temp.runs, "temp.number"=temp.number, "temp.strategy"=temp.strategy)
  
  # add model settings into summary table
  temp.dat <- data.frame(SpeciesID=myRespName, model="BIOMOD", strategy=temp.strategy,
                         No.runs=temp.runs, No.points=temp.number, 
                         min.distance=NA, run.time=temp.time)
  model.settings <- rbind(model.settings, temp.dat)
  
  rm(temp.strategy, temp.runs, temp.number, temp.min.dist, temp.time)
  
  print(model.settings)
  print("The strategy used for GLM/GAM can be used for other models than listed.")
  
  
  # save model setting summary for later
  save(model.settings, file=paste0(here::here(), "/results/", Taxon_name, "/BackgroundData_modelSettings_", Taxon_name, "_", myRespName, ".RData"))
  
  # save background data
  bg.list <- list(bg.glm, bg.mars, bg.rf, bg.biomod) #bg.mda,
  names(bg.list) <- c("bg.glm", "bg.mars", "bg.rf", "bg.biomod") #"bg.mda",
  save(bg.list, file=paste0(here::here(), "/results/", Taxon_name, "/PA_Env_", Taxon_name, "_", myRespName, ".RData"))
  # write.csv(bg.glm, file=paste0(here::here(), "/results/", Taxon_name, "/PA_Env_GLM_", Taxon_name, "_", myRespName, ".csv"), row.names = F)
  # write.csv(bg.mars, file=paste0(here::here(), "/results/", Taxon_name, "/PA_Env_MARS_", Taxon_name, "_", myRespName, ".csv"), row.names = F)
  # #write.csv(bg.mda, file=paste0(here::here(), "/results/", Taxon_name, "/PA_Env_MDA_", Taxon_name, "_", myRespName, ".csv"), row.names = F)
  # write.csv(bg.rf, file=paste0(here::here(), "/results/", Taxon_name, "/PA_Env_RF_", Taxon_name, "_", myRespName, ".csv"), row.names = F)
}
          
# at the end
stopImplicitCluster()
          
          
