# rdsmart 2.1.3

* Fix bug in selection of mlr3 models using `disaggregate()` argument `method.model`
* Correct `disaggregate()` and `sort_stack_values()` documentation
* Remove dataset `dalrymple_realisations` to remove dependency on `raster` package
* Add package `mlr3proba` to Imports

# rdsmart 2.1.2

## Bug fixes

* Spatial disaggregation using `method.sample = "by_area"` works properly again. Replaced use of deprecated function `terra::area()` with `terra::expanse()`. This *should* mean that rdsmart should be compatible with spatial data in a longitude/latitude CRS (but it's not tested).

# rdsmart 2.1.1

* Repackaged demonstration data `dalrymple_covariates` and `dalrymple_polygons` as (wrapped) SpatRaster and SpatVector objects respectively. Use `terra::unwrap()` to unwrap before use. Updated example code accordingly.
* Added vignette by Matthew Coghill to explain major underlying changes to package from v2.1.0.

# rdsmart 2.1.0

* Replaced `raster` package functions with equivalent `terra` package functions. This requires at least `terra` version 1.1-17 (CRAN).
* Wrote C++ sorting and ordering algorithms
* Commented out the `.onAttach()` message since `terra` takes care of previous issues
* Updated `NAMESPACE` and `DESCRIPTION` to handle package updates. Bumped minimum R version to 3.5.0 to mirror `terra` minimum version. Removed `C50` import since `disaggregate` will install and load external model packages.
* Cleaned up man folder

# rdsmart 2.0.28

* Created `NEWS.md`.
* Updated doi links in documentation for `disaggregate()`, `dsmart()` and `summarise()` to https.
* Implemented `shannon_entropy()` to compute Shannon's entropy (rather than compute it inside `summarise()`) and fixed a bug in its previous implementation.
* Improved documentation for `confusion_index()`.
* The function `summarise()` now uses `confusion_index()` and `shannon_entropy()` to compute the confusion index and Shannon's entropy rather than computing them itself.
* Implemented `.onAttach()` to display a message when the package is loaded to explain that writing rasters may fail under GDAL >= 3 and PROJ >= 6.