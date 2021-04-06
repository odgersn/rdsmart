# library(terra)
# library(tidyverse)
# library(mlr3verse)

data("dalrymple_covariates")
data("dalrymple_composition")
data("dalrymple_lookup")
data("dalrymple_observations")
data("dalrymple_polygons")

covariates <- terra::rast(dalrymple_covariates)
polygons <- terra::vect(dalrymple_polygons)
composition <- dalrymple_composition
rate <- 15
reals <- 3
observations <- dalrymple_observations
method.sample = "by_polygon"
method.allocate = "weighted"
method.model = NULL
args.model = NULL
strata = NULL
outputdir = getwd()
stub = NULL
type = "prob"
predict = TRUE
factors = NULL
nprob = 3

# for(i in list.files(file.path(getwd(), "R"), full.names = TRUE)) source(i)
# for(i in list.files(file.path(getwd(), "src"), full.names = TRUE)) Rcpp::sourceCpp(i)

dsm <- dsmart(covariates, polygons, composition, 15, 3, observations, type = type, stub = type)
