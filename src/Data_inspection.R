#- - - - - - - - - - - - - - - - -#
#  Data inspection (correlation)  #
#                                 #
#     author: Romy Zeiss          #
#       date: 18.08.2021          #
#- - - - - - - - - - - - - - - - -#

library(here)  #instead of setwd()
library(tidyverse)
library(corrplot)

I:\eie\==PERSONAL\Macroecology\Students\Jessica\ANgela

# load data
data <- readr::read_csv("I:/eie/==PERSONAL/Macroecology/Students/Jessica/ANgela/2k.txt")

str(data)
colnames(data)

#- - - - - - - - - - - - - - - - - - - - -
## Check the correlations between variables ####

cor.species <- cor(data[,3:47])
cor.env <- cor(data[complete.cases(data$NDVI),48:72])
# note: some NA in NDVI -> have to be removed

# species
pdf(file=here::here("results","Corrplot_earthworm_species.pdf"))
corrplot.mixed(cor.species, number.cex=0.5, cl.cex=0.4, tl.cex=0.5)
dev.off()

# environmental variables
pdf(file=here::here("results","Corrplot_earthworm_env.pdf"))
corrplot.mixed(cor.env, number.cex=0.5, cl.cex=0.4, tl.cex=0.5)
dev.off()

