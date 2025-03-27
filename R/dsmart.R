#' Disaggregating and harmonising soil map units through resampled 
#' classification trees
#' 
#' `dsmart` performs the spatial disaggregation of a soil choropleth map. 
#' The underlying functions are `disaggregate`, which performs the spatial 
#' downscaling, and `summarise`, which computes the soil class 
#' probabilities and n-most-probable soil classes.
#' 
#' @param covariates A SpatRaster of *scorpan* environmental 
#'   covariates to calibrate the `C50` classification trees against. See 
#'   *Details* for more information.
#' @param polygons A `SpatVector` containing the soil map unit polygons
#'   that will be disaggregated. The first field of the data frame must be an 
#'   integer that identifies each polygon.
#' @param composition A `data.frame` that contains information on the 
#'   soil-class composition of each polygon in `polygons`. Each row 
#'   contains information about one soil class component of one polygon, which 
#'   belongs to one soil map unit. First field contains the integer that 
#'   identifies the polygon. Second field contains a code that identifies the 
#'   soil map unit that the polygon belongs to.
#'   
#'   If `strata = NULL` (the default), third column contains a code that 
#'   identifies the soil class and fourth column contains a number in the range 
#'   `(0, 100)` that identifies the proportion of the **map unit** 
#'   that the soil class corresponds to. See the example data 
#'   `data(dalrymple_composition)`.
#'   
#'   If `strata` is a SpatRaster, third column contains an integer 
#'   that identifies the stratum in `strata`, fourth column contains a code
#'   that identifies the soil class and fifth column contains a number in the
#'   range `(0, 100)` that identifies the proportion of the
#'   **stratum** that the soil class corresponds to.
#' @param rate An integer that identifies the number of virtual samples to draw 
#'   from each polygon in each realisation. If `method.sample = 
#'   "by_polygon"`, the number of samples to draw from each polygon in 
#'   `polygons`. If `method.sample = "by_area"`, the sampling density 
#'   in number of samples per square kilometer.
#' @param reals An integer that identifies the number of realisations of the 
#'   soil class distribution that DSMART should compute.
#' @param observations *optional* A `data.frame` that contains actual 
#'   observations of the soil class at locations across the soil map area. These
#'   data augment the virtual samples and are used in each realisation. Each row
#'   contains information about one soil class observation. First and second 
#'   fields contain the *x-* and *y-*components of the observation's 
#'   spatial location. Third field is the soil class code. See *Details*.
#' @param method.sample Identifies the sampling method. Valid values are 
#'   `"by_polygon"` (the default), in which case the same number of samples
#'   are taken from each polygon; or `"by_area"`, in which case the number 
#'   of samples per polygon depends on the area of the polygon.
#' @param method.allocate Method of allocation of virtual samples to soil 
#'   classes. Valid values are `"weighted"`, for weighted-random allocation
#'   to a soil class from within the virtual sample's map unit; 
#'   `"random_mapunit"`, for completely random allocation to a soil class 
#'   from within the virtual sample's map unit; and `"random_all"`, for 
#'   completely random allocation to a soil class from within the entire map 
#'   area.
#' @param method.model Method to be used for the classification model. If no
#'   value is passed, a C5.0 decision tree is built. Otherwise, the value must
#'   match a valid 'learner' argument in the mlr3::lrn() function.
#' @param method.args A named list of arguments to be passed to the mlr3 learner
#'   object. The list will modify the learner's 'param_set' which controls the 
#'   behavior of the model. Named arguments are passed directly to the train 
#'   function and predictive model. To view a model's given parameter set, use
#'   `mlr3::lrn(method.model)$param_set`
#' @param strata *optional* An integer-valued SpatRaster that will 
#'   be used to stratify the allocation of virtual samples to soil classes. 
#'   Integer values could represent classes of slope position (e.g. crest, 
#'   backslope, footslope, etc.) or land use (e.g. cropland, native vegetation, 
#'   etc.) or some other variable deemed to be an important discriminator of the
#'   occurrence of soil classes within a map unit.
#' @param nprob At any location, disaggregated soil class predictions can be 
#'   ranked according to their probabilities of occurrence. `rdsmart` can 
#'   map the class predictions, and their probabilities, at any rank. 
#'   `nprob` is an integer that identifies the number of probability ranks 
#'   to map. For example, if `nprob = 3`, DSMART will map the first-, 
#'   second- and third-most-probable soil classes and their probabilities of 
#'   occurrence.
#' @param outputdir A character string that identifies the location of the main 
#'   output directory. The folder `output` and its subfolders will be 
#'   placed here. Default is the current working directory, `getwd()`.
#' @param stub *optional* A character string that identifies a short name 
#'   that will be prepended to all output.
#' @param factors A character vector with the names of the covariates that should
#'   be treated as factors.
#' @param type A character vector to specify the type of the predictions. By
#'   default, "response" class predictions are used. Each realisation will 
#'   produce a map of soil classes, and class probabilities will be calculated 
#'   from the frequency of each class across the realisations. If set to "prob", 
#'   each realisation will produce a SpatRaster with the probabilities of each 
#'   class. The final class probabilities will be calculated by averaging the 
#'   class probabilities across the realisation.
#'   
#' @return A list that aggregates metadata about the current run of
#'   `disaggregate` and `summarise`.
#'   
#' @references McBratney, A.B., Mendonca Santos, M. de L., Minasny, B., 2003. On
#'   digital soil mapping. Geoderma 117, 3--52. doi: 
#'   [10.1016/S0016-7061(03)00223-4](https://doi.org/10.1016/S0016-7061(03)00223-4)
#'   
#'   Odgers, N.P., McBratney, A.B., Minasny, B., Sun, W., Clifford, D., 2014. 
#'   DSMART: An algorithm to spatially disaggregate soil map units, *in:* 
#'   Arrouays, D., McKenzie, N.J., Hempel, J.W., Richer de Forges, A., 
#'   McBratney, A.B. (Eds.), GlobalSoilMap: Basis of the Global Spatial Soil 
#'   Information System. Taylor & Francis, London, pp. 261--266.
#'   
#'   Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., 2014. 
#'   Disaggregating and harmonising soil map units through resampled 
#'   classification trees. Geoderma 214, 91--100. doi: 
#'   [10.1016/j.geoderma.2013.09.024](https://doi.org/10.1016/j.geoderma.2013.09.024)
#'   
#' @examples 
#' # Load datasets
#' data(dalrymple_composition)
#' data(dalrymple_covariates)
#' data(dalrymple_observations)
#' data(dalrymple_polygons)
#' 
#' # Run dsmart without adding observations
#' dsmart(unwrap(dalrymple_covariates), unwrap(dalrymple_polygons), dalrymple_composition,
#'  rate = 15, reals = 10)
#' 
#' # Run dsmart with extra observations
#' dsmart(unwrap(dalrymple_covariates), unwrap(dalrymple_polygons), dalrymple_composition,
#'  observations = dalrymple_observations, rate = 15, reals = 10)
#' 
#' @export
#' 
dsmart <- function(covariates, polygons, composition, rate = 15, reals = 100, 
                   observations = NULL, method.sample = "by_polygon", 
                   method.allocate = "weighted",
                   method.model = NULL, args.model = NULL,
                   strata = NULL, nprob = 3,
                   outputdir = getwd(), stub = NULL,
                   factors = NULL, type = "response")
{
  # Create list to store output
  output <- base::list()
  
  # Save start time
  output$timing <- base::list(start = Sys.time())
  
  # Set stub to "" if NULL
  if(is.null(stub))
  {
    stub <- ""
    
  } else if (stub == "") {
  
    stub <- ""
    
  } else if(!(substr(stub, nchar(stub), nchar(stub)) == "_")) {
    
    stub <- paste0(stub, "_")
  }
  
  # Create output directory
  outputdir <- file.path(outputdir)
  dir.create(file.path(outputdir, "output"), showWarnings = FALSE)
  
  # Write function call to text file as a means of preserving the parameters
  # that were submitted to the dsmart function for the current run.
  # We should think of a less clunky way to do it---usually the call is returned
  # from the called function as a list element together with other output.
  base::match.call() %>%
    base::deparse() %>%
    base::write(file = file.path(outputdir, "output", "dsmart_function_call.txt"))
  
  # Carry out spatial disaggregation
  output$disaggregate <- disaggregate(covariates, polygons, composition,
                                      rate = rate, reals = reals, 
                                      observations = observations,
                                      method.sample = method.sample,
                                      method.allocate = method.allocate,
                                      method.model = method.model,
                                      args.model = args.model,
                                      strata = strata, outputdir = outputdir,
                                      stub = stub, factors = factors, type = type)
  
  
  real_files <- list.files(path = file.path(outputdir, "output", "realisations"),
                           pattern = paste0(stub, ".*.tif$"), full.names = TRUE)
  
  # If raw class predictions are used, load realisations to SpatRaster
  if(type != "prob")
  {
    realisations <- terra::rast(real_files)
  } else {
    # If probabilistic class predictions are used, load realisations to a list
    # of SpatRaster objects.
    realisations <- lapply(real_files, terra::rast)
  }
  
  # Load lookup table
  lookup <- read.table(file.path(outputdir, "output", paste0(stub,"lookup.txt")),
                       header = TRUE, sep = ",")
  
  # Summarise the results of the spatial disaggregation
  output$summarise <- summarise(realisations, lookup, n.realisations = reals,
                                nprob = nprob, 
                                outputdir = outputdir, stub = stub, type = type)
  
  # Save finish time
  output$timing$finish <- Sys.time()
  
  # Cleanup temp files
  suppressWarnings(terra::tmpFiles(old = TRUE, remove = TRUE))
  
  # Return output
  return(output)
}
