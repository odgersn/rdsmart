library(terra)
library(tidyverse)
library(mlr3verse)

load("G:/R/dsmart/data/dalrymple_covariates.rda")
load("G:/R/dsmart/data/dalrymple_composition.rda")
load("G:/R/dsmart/data/dalrymple_lookup.rda")
load("G:/R/dsmart/data/dalrymple_observations.rda")
load("G:/R/dsmart/data/dalrymple_polygons.rda")
dalrymple_covariates <- rast(dalrymple_covariates)
dalrymple_polygons <- vect(dalrymple_polygons)

covariates <- dalrymple_covariates
polygons <- dalrymple_polygons
composition <- dalrymple_composition
rate <- 15
reals <- 10
observations <- dalrymple_observations
method.sample = "by_polygon"
method.allocate = "weighted"
method.model = NULL
args.model = NULL
strata = NULL
outputdir = getwd()
stub = NULL
type = "response"
predict = TRUE
factors = NULL
nprob = 3

for(i in list.files(file.path(getwd(), "R"), full.names = TRUE)) source(i)
for(i in list.files(file.path(getwd(), "src"), full.names = TRUE)) Rcpp::sourceCpp(i)

dsm <- dsmart(covariates, polygons, composition, 15, 10, observations, type = type, stub = type)
