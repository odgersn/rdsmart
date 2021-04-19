# Load datasets
data("dalrymple_covariates")
data("dalrymple_composition")
data("dalrymple_lookup")
data("dalrymple_observations")
data("dalrymple_polygons")

# Generate variables for function call
covariates <- terra::rast(dalrymple_covariates)
polygons <- terra::vect(dalrymple_polygons)
composition <- dalrymple_composition
rate <- 15
reals <- 3
observations <- dalrymple_observations
method.sample <- "by_polygon"
method.allocate <- "weighted"
strata <- NULL
outputdir <- getwd()
stub <- NULL
type <- "prob"
predict <- TRUE
factors <- NULL
nprob <- 3

# Generate model parameters. method.model = NULL defaults to C5.0 tree models,
# and args.model is a named list of model parameters
method.model <- "classif.C50"
args.model <- list(CF = 0.5)

# To view model specific parameters:
mlr3::lrn(method.model)$param_set

# Run DSMART algorithm
dsm <- dsmart(covariates, polygons, composition, 15, 3, observations, type = type, stub = type)
