library(raster)

setwd("/home/brendo/myWork/dsmart/rPackage/dsmart/pkg/data")
list.files()

#dsmart outputs
load("dsmartOutMaps.rda")
dsmartOutMaps

#dsmart covariates
load("dsT_covariates.rda")
dsT_covariates

#dsmart composition
load("dsT_composition.rda")
dsT_composition

#dsmart polygons
load("dsT_polygons.rda")
dsT_polygons

#dsmart polygons
load("dsT_lookup.rda")
dsT_lookup


#Run Functions
dsmartR(rLocs= dsT_covariates, nprob = 2, sepP=TRUE, lookup=dsT_lookup , cpus=1)


