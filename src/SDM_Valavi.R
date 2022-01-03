#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(mgcv) # for GAM
library(gam)  # for GLM (!)


# load all datasets
trainig <- read.csv(file=paste0(here::here(), "/results/TrainingData", Taxon_name, ".csv"))
testing_pa <- read.csv(file=paste0(here::here(), "/results/TestingData_pa_", Taxon_name, ".csv"))
testing_env <- read.csv(file=paste0(here::here(), "/results/TestingData_env_", Taxon_name, ".csv"))
validation.data <- read.csv(file=paste0(here::here(), "/results/ValidationData", Taxon_name, ".csv"))

# calculate Infinitely Weighted Logistic Regression (IWLR) for GLM and GAM models suggested by Fithian and Hastie (2013)
iwp <- (10^6)^(1-training$occ)


#- - - - - - - - - - - - - - - - - - - - -
## GAM ####
#- - - - - - - - - - - - - - - - - - - - -

# calculating the case weights (equal weights)
# the order of weights should be the same as presences and backgrounds in the training data
prNum <- as.numeric(table(training$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)

# define formula
form <- paste0("occ ~ ", paste0(paste0("s(", covars, ")"), collapse=" + "))

# general settings
tmp <- Sys.time()
set.seed(32639)

# run model
gm <- mgcv::gam(formula = as.formula(form),
                data = training,
                family = binomial(link = "logit"),
                weights = wt,
                method = "REML")

# get running time
Sys.time() - tmp

# check the appropriateness of Ks
# the model parameter k should not be higher than the number of unique values in each covariate
gam.check(gm)

plot(gm, pages = 1, rug = TRUE, shade = TRUE)


#- - - - - - - - - - - - - - - - - - - - -
## GLM ####
#- - - - - - - - - - - - - - - - - - - - -

# calculating the weights
# the order of weights should be the same as presences and backgrounds in the training data
prNum <- as.numeric(table(training$occ)["1"]) # number of presences
bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
wt <- ifelse(training$occ == 1, 1, prNum / bgNum)

# the base glm model with linear terms
lm1 <- glm(occ ~., data = training, weights = wt, family = binomial(link = "logit"))
summary(lm1)

# model scope for subset selection
mod_scope <- list("cti" = ~1 + cti + poly(cti, 2),
                  "disturb" = ~1 + disturb + poly(disturb, 2),
                  "mi" = ~1 + mi + poly(mi, 2),
                  "rainann" = ~1 + rainann + poly(rainann, 2),
                  "raindq" = ~1 + raindq + poly(raindq, 2),
                  "rugged" = ~1 + rugged + poly(rugged, 2),
                  "soildepth" = ~1 + soildepth + poly(soildepth, 2),
                  "soilfert" = ~1 + soilfert + poly(soilfert, 2),
                  "solrad" = ~1 + solrad + poly(solrad, 2),
                  "tempann" = ~1 + tempann + poly(tempann, 2),
                  "topo" = ~1 + topo + poly(topo, 2),
                  "vegsys" = ~1 + vegsys)
tmp <- Sys.time()
set.seed(32639)
lm_subset <- step.Gam(object = lm1,
                      scope = mod_scope,
                      direction = "both",
                      data = training, # this is optional - see details below
                      trace = FALSE)
Sys.time() - tmp

summary(lm_subset)


#- - - - - - - - - - - - - - - - - - - - -
## Regularized regressions: lasso and ridge regression ####
##- - - - - - - - - - - - - - - - - - - - -





