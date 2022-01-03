# see Tutorial https://cran.r-project.org/web/packages/SSDM/vignettes/SSDM.html#ensemble_species_distribution_models_(esdms)

gui()

## Environmental (explanatory) variables as raster file
Env <- stack(paste0(here::here(), "/results/EnvPredictor_", Taxon_name, ".grd"))

# crop to Europe
Env <- crop(Env, extent_Europe)
Env

# Load occurrence data, thinning of data to resolution of environmental data
Occ <- read.csv(file=paste0(here::here(), "/results/Occurrences_", Taxon_name, ".csv"))

Occ <- data.frame(SPECIES=Occ$SpeciesID, LONGITUDE=Occ$decimalLongitude, LATITUDE=Occ$decimalLatitude)
Occ <- Occ[complete.cases(Occ),]

## Define species names
# check first if there are enough occurrence data for each species
sp.occ <- as.data.frame(table(Occ$SPECIES))

# define species 
sp.occ <- sp.occ[sp.occ$Freq>=100,"Var1"]

# remove species with less than XX occurrences
Occ <- Occ[Occ$SPECIES %in% sp.occ,]
head(Occ)

# define which algorithms should be used
#setAlgorithms <- c("CTA", "MARS", "MAXENT")
setAlgorithms <- "all"

# run stacked SDMs
SSDM <- stack_modelling(setAlgorithms, Occ, Env, rep = 1, ensemble.thresh = 0, metric="TSS", axes.metric = "AUC",
                        Xcol = "LONGITUDE", Ycol = "LATITUDE", ensemble.metric = "AUC",
                        PA=list(nb=10000, strat="disk"), #create pseudo-absences
                        cv="k-fold", cv.param=c(5,2), #cross-validation
                        Spcol = "SPECIES", method = "pSSDM", verbose = FALSE,  endemism = c("CWEI", "NbOcc"),
                        # for the individual algorithms:
                        test="AIC", epsilon=10e-08, # GLM & GAM 
                        maxit=500, # GLM & GAM & ANN
                        degree=2, # MARS
                        trees=250, final.leave=2, # GBM & RF
                        algocv=5, tresh.shrink=1e-03) # GBM & CTA

save(SSDM, file=paste0("SSDM_", Taxon_name, ".RData"))

plot(SSDM@diversity.map, main = paste0("SSDM for ", Taxon_name, "\nwith ", paste(setAlgorithms, collapse=" & "), " algorithms"))

# Model accuracy
knitr::kable(SSDM@evaluation)

# Six evaluation metrics are computed: 
# (1) the species richness error, i.e. the difference between the predicted and observed species richness; 
# (2) the assemblage prediction success, i.e. the proportion of correct predictions; 
# (3) the assemblage Cohen’s kappa, i.e. the proportion of specific agreement; 
# (4) the assemblage specificity, i.e. the proportion of true negatives (species that are both predicted and observed as being absent); 
# (5) the assemblage sensitivity, i.e. the proportion of true positives (species that are both predicted and observed as present); and 
# (6) the Jaccard index, a widely used metric of community similarity.

# Variable importance
knitr::kable(SSDM@variable.importance)

# endemism map
plot(SSDM@endemism.map, main = paste0("Endemism map for ", Taxon_name, "\nwith ", paste(setAlgorithms, collapse=" & "), " algorithms"))

# summary plots
plot(SSDM)
