library(terra)

# Load polygons to a SpatVector
dalrymple_polygons <-
  vect(here::here("data-raw", "data", "polygons", "dlr_polys_alb_7km.shp")) %>% 
  wrap()

# Write data to rda file
usethis::use_data(dalrymple_polygons, overwrite = TRUE)
