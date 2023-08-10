library(terra)

# Load covariate rasters to a SpatRaster
dalrymple_covariates <-
  rast(list.files(here::here("data-raw", "data", "covariates"), pattern = "tif$", full.names = TRUE)) %>% 
  wrap()

# Write data to rda file
usethis::use_data(dalrymple_covariates, overwrite = TRUE)
