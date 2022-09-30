function (data, NbRunEval = NULL, DataSplit, DataSplitTable = NULL, 
    Prevalence = 0.5, weighting.score, models, tune = FALSE, 
    modeling.id = as.character(format(Sys.time(), "%s")), models.options = NULL, 
    which.biva = NULL, parallel = FALSE, cleanup = FALSE, Yweights = NULL) 
{
    if (!weighting.score %in% c("AUC", "TSS", "Boyce", "Kappa", 
        "SomersD")) {
        stop("weighting score not supported! Choose one of the following: AUC, TSS, Boyce, Kappa or SomersD")
    }
    if (data@has.data.eval) {
        stop("Evaluation with independant data is not supported yet!")
    }
    if (weighting.score == "AUC" | weighting.score == "SomersD") {
        models.eval.meth <- "ROC"
    }
    if (weighting.score == "Kappa") {
        models.eval.meth <- "KAPPA"
    }
    if (weighting.score == "TSS") {
        models.eval.meth <- "TSS"
    }
    if (weighting.score == "Boyce") {
        models.eval.meth <- "ROC"
    }
    if (is.null(models.options)) {
        models.options <- BIOMOD_ModelingOptions()
        models.options@GBM$n.trees <- 1000
        models.options@GBM$interaction.depth <- 4
        models.options@GBM$shrinkage <- 0.005
        models.options@GAM$select <- TRUE
        models.options@CTA$control$cp <- 0
        models.options@ANN$size <- 8
        models.options@ANN$decay <- 0.001
        models.options@MARS$interaction.level <- 0
        models.options@MARS$nprune <- 2
        models.options@MAXENT.Phillips$product <- FALSE
        models.options@MAXENT.Phillips$threshold <- FALSE
        models.options@MAXENT.Phillips$betamultiplier <- 0.5
        models.options@GLM$test <- "none"
    }
    if ("MAXENT.Phillips" %in% models) {
        if (!file.exists(paste(models.options@MAXENT.Phillips$path_to_maxent.jar, 
            "maxent.jar", sep = "/"))) {
            stop("maxent.jar file not found!")
        }
    }
    models <- sort(models)
    iniwd <- getwd()
    on.exit(setwd(iniwd))
    dir.create(paste("./ESM.BIOMOD.output", data@sp.name, sep = "_"))
    newwd <- paste(getwd(), "/ESM.BIOMOD.output_", data@sp.name, 
        sep = "")
    setwd(newwd)
    combinations <- combn(colnames(data@data.env.var), 2)
    if (is.null(which.biva)) {
        which.biva <- 1:ncol(combinations)
    }
    mydata <- data
    if (is.null(DataSplitTable)) {
        mod.prep.dat <- .Models.prepare.data(mydata, NbRunEval, 
            DataSplit, Yweights = NULL, Prevalence = Prevalence, 
            do.full.models = TRUE)
        if (length(dim(mod.prep.dat[[1]]$calibLines)) == 3) {
            calib.lines <- mod.prep.dat[[1]]$calibLines[, , 1]
        }
        if (length(dim(mod.prep.dat[[1]]$calibLines)) == 2) {
            calib.lines <- mod.prep.dat[[1]]$calibLines
        }
        rm(mod.prep.dat)
    }
    else {
        calib.lines <- DataSplitTable
    }
    if (is.null(NbRunEval)) {
        if (ncol(calib.lines) > 1) {
            if (sum(!as.data.frame(calib.lines)[, ncol(calib.lines)]) == 
                0) {
                NbRunEval <- ncol(calib.lines) - 1
            }
            else {
                NbRunEval <- ncol(calib.lines)
            }
        }
        else {
            if (sum(!calib.lines) == 0) {
                NbRunEval <- ncol(calib.lines) - 1
            }
            else {
                NbRunEval <- ncol(calib.lines)
            }
        }
    }
    mymodels <- list()
    if (parallel == FALSE) {
        for (k in which.biva) {
            mydata@data.env.var <- data@data.env.var[, colnames(data@data.env.var) %in% 
                combinations[, k]]
            mydata@sp.name <- paste("ESM.BIOMOD", k, sep = ".")
            if (tune == TRUE) {
                models.options <- BIOMOD_tuning(data = mydata, 
                  models = models[models != "RF"], models.options = models.options, 
                  Yweights = Yweights)$models.options
            }
            mymodels[[k]] <- "failed"
            try(mymodels[[k]] <- BIOMOD_Modeling(data = mydata, 
                models = models, models.options = models.options, 
                models.eval.meth = models.eval.meth, DataSplitTable = as.matrix(calib.lines), 
                Prevalence = Prevalence, rescal.all.models = FALSE, 
                do.full.models = TRUE, VarImport = 0, modeling.id = modeling.id, 
                Yweights = NULL))
            if (cleanup != FALSE) {
                removeTmpFiles(h = cleanup)
            }
        }
    }
    if (parallel == TRUE) {
        mymodels <- foreach(k = which.biva, .packages = c("biomod2", 
            "raster")) %dopar% {
            setwd(newwd)
            mydata@data.env.var <- data@data.env.var[, colnames(data@data.env.var) %in% 
                combinations[, k]]
            mydata@sp.name <- paste("ESM.BIOMOD", k, sep = ".")
            if (cleanup != FALSE) {
                removeTmpFiles(h = cleanup)
            }
            if (tune == TRUE) {
                models.options <- BIOMOD_tuning(data = mydata, 
                  models = models[models != "RF"], models.options = models.options)$models.options
            }
            BIOMOD_Modeling(data = mydata, models = models, models.options = models.options, 
                models.eval.meth = models.eval.meth, DataSplitTable = as.matrix(calib.lines), 
                Prevalence = Prevalence, rescal.all.models = TRUE, 
                do.full.models = TRUE, VarImport = 0, modeling.id = modeling.id, 
                Yweights = Yweights)
        }
    }
    output <- list(modeling.id = modeling.id, models. = grep(modeling.id, 
        gtools::mixedsort(list.files(getwd(), "models.out", recursive = TRUE, 
            full.names = TRUE)), value = TRUE), models = models, 
        calib.lines = calib.lines, NbRunEval = NbRunEval, data = data, 
        wd = getwd(), which.biva = which.biva, mymodels = mymodels)
    save(output, file = paste("ESM_Modeling..models", modeling.id, 
        "out", sep = "."))
    setwd(iniwd)
    return(output)
}
