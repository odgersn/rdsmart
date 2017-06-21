#' Disaggregating and harmonising soil map units through resampled 
#' classification trees
#' 
#' This function, together with companion function \code{dsmartR} implements the
#' DSMART (Disaggregating and harmonising soil map units through resampled 
#' classification trees) algorithm as described in Odgers et al. (2014). This is
#' the workhorse function that involves multiple resampling, C5 decision tree 
#' model fitting, and subsequent mapping in order to realize potential candidate
#' soil classes within aggregated soil mapping units. There is also added 
#' facility to incorporate point observed data into the algorithm too.
#' 
#' There is no prescription for the layers that compose \code{covariates} other 
#' than that they should represent the \emph{scorpan} factors (McBratney 
#' \emph{et al.}, 2003) as effectively as possible. The order of the layers in 
#' the RasterStack is not important. See \code{data(dsT_covariates)} for an 
#' example.
#' 
#' DSMART assumes, but does not currently check, that \code{covariates}, 
#' \code{polygons} and any \code{observations} are projected to the same 
#' coordinate reference system.
#' 
#' DSMART assumes that the soil classes in \code{composition} and those in 
#' \code{observations} belong to the same taxonomic level of the same soil 
#' classification system. For example if the soil classes in \code{composition} 
#' are all soil series, those in \code{observations} should also be series 
#' rather than, for example, orders or great groups.
#' 
#' \code{dsmart} produces a number of outputs which are saved to subdirectories 
#' in the current working directory. The base folder for the outputs is 
#' \code{output}. Inside this folder, rasters of the realisations are saved into
#' the \code{realisations} subfolder and the classification trees are saved into
#' the \code{trees} subfolder.
#' 
#' @param covariates A \code{RasterStack} of \emph{scorpan} environmental 
#'   covariates to calibrate the \code{C50} classification trees against. See 
#'   \emph{Details} for more information.
#' @param polygons A \code{SpatialPolygonsDataFrame} containing the soil map 
#'   unit polygons that will be disaggregated. The first field of the data frame
#'   must be an integer that identifies each polygon.
#' @param composition A \code{data.frame} that contains information on the 
#'   soil-class composition of each polygon in \code{polygons}. Each row 
#'   contains information about one soil class component of one polygon, which 
#'   belongs to one soil map unit. First field contains the integer that 
#'   identifies the polygon. Second field contains a code that identifies the 
#'   soil map unit that the polygon belongs to.
#'   
#'   If \code{strata = NULL} (the default), third column contains a code that 
#'   identifies the soil class and fourth column contains a number in the range 
#'   \code{(0, 100)} that identifies the proportion of the \strong{map unit} 
#'   that the soil class corresponds to. See the example data 
#'   \code{data(dalrymple_composition)}.
#'   
#'   If \code{strata} is a \code{RasterLayer}, third column contains an integer 
#'   that identifies the stratum in \code{strata}, fourth column contains a code
#'   that identifies the soil class and fifth column contains a number in the 
#'   range \code{(0, 100)} that identifies the proportion of the 
#'   \strong{stratum} that the soil class corresponds to.
#' @param rate An integer that identifies the number of virtual samples to draw 
#'   from each polygon in each realisation. If \code{method.sample = 
#'   "by_polygon"}, the number of samples to draw from each polygon in 
#'   \code{polygons}. If \code{method.sample = "by_area"}, the sampling density 
#'   in number of samples per square kilometre.
#' @param reals An integer that identifies the number of realisations of the 
#'   soil class distribution that DSMART should compute.
#' @param observations \emph{Optional} A \code{data.frame} that contains actual 
#'   observations of the soil class at locations across the soil map area. These
#'   data augment the virtual samples and are used in each realisation. Each row
#'   contains information about one soil class observation. First and second 
#'   fields contain the \emph{x-} and \emph{y-}components of the observation's 
#'   spatial location. Third field is the soil class code. See \emph{Details}.
#' @param method.sample Identifies the sampling method. Valid values are 
#'   \code{"by_polygon"} (the default), in which case the same number of samples
#'   are taken from each polygon; or \code{"by_area"}, in which case the number 
#'   of samples per polygon depends on the area of the polygon.
#' @param method.allocate Method of allocation of virtual samples to soil 
#'   classes. Valid values are \code{"weighted"}, for weighted-random allocation
#'   to a soil class from within the virtual sample's map unit; 
#'   \code{"random_mapunit"}, for completely random allocation to a soil class 
#'   from within the virtual sample's map unit; and \code{"random_all"}, for 
#'   completely random allocation to a soil class from within the entire map 
#'   area.
#' @param strata \emph{optional} An integer-valued \code{RasterLayer} that will 
#'   be used to stratify the allocation of virtual samples to soil classes. 
#'   Integer values could represent classes of slope position (e.g. crest, 
#'   backslope, footslope, etc.) or land use (e.g. cropland, native vegetation, 
#'   etc.) or some other variable deemed to be an important discriminator of the
#'   occurrence of soil classes within a map unit.
#' @param outputdir A character string that identifies the location of the main 
#'   output directory. The folder \code{output} and its subfolders will be 
#'   placed here. Default is the current working directory, \code{getwd()}.
#' @param stub \emph{optional} A character string that identifies a short name 
#'   that will be prepended to all output.
#' @param cpus An integer that identifies the number of CPU processors to use 
#'   for parallel processing.
#'   
#' @examples
#' # Load datasets
#' data(dalrymple_composition)
#' data(dalrymple_covariates)
#' data(dalrymple_observations)
#' data(dalrymple_polygons)
#' 
#' # Run disaggregate without adding observations
#' disaggregate(dalrymple_covariates, dalrymple_polygons, dalrymple_composition,
#'  rate = 15, reals = 10, cpus = 6)
#' 
#' # Run disaggregate with extra observations
#' disaggregate(dalrymple_covariates, dalrymple_polygons, dalrymple_composition,
#'  observations = dalrymple_observations, rate = 15, reals = 10)
#' 
#' @references McBratney, A.B., Mendonca Santos, M. de L., Minasny, B., 2003. On
#'   digital soil mapping. Geoderma 117, 3--52. doi: 
#'   \href{http://dx.doi.org/10.1016/S0016-7061(03)00223-4}{10.1016/S0016-7061(03)00223-4}
#'   
#'   
#'   
#'   
#'   
#'   
#'   Odgers, N.P., McBratney, A.B., Minasny, B., Sun, W., Clifford, D., 2014. 
#'   DSMART: An algorithm to spatially disaggregate soil map units, \emph{in:} 
#'   Arrouays, D., McKenzie, N.J., Hempel, J.W., Richer de Forges, A., 
#'   McBratney, A.B. (Eds.), GlobalSoilMap: Basis of the Global Spatial Soil 
#'   Information System. Taylor & Francis, London, pp. 261--266.
#'   
#'   Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., 2014. 
#'   Disaggregating and harmonising soil map units through resampled 
#'   classification trees. Geoderma 214, 91--100. doi: 
#'   \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{10.1016/j.geoderma.2013.09.024}
#'   
#'   
#'   
#'   
#'   
#'   
#' @export

disaggregate <- function(covariates, polygons, composition, rate = 15,
                         reals = 100, observations = NULL,
                         method.sample = "by_polygon", 
                         method.allocate = "weighted", strata = NULL,
                         outputdir = getwd(), stub = NULL, cpus = 1)
{
  # Check arguments before proceeding
  messages <- c("Attention is required with the following arguments:\n")
  if(!(class(covariates) == "RasterStack"))
  {
    messages <- append(messages, "'covariates': Not a valid RasterStack.\n")
  }
  if(!(class(polygons) == "SpatialPolygonsDataFrame"))
  {
    messages <- append(messages, 
                       "'polygons': Not a valid SpatialPolygonsDataFrame.\n")
  }
  if(!(class(composition) == "data.frame"))
  {
    messages <- append(messages, "'composition': Not a valid data.frame.\n")
  }
  if(rate <= 0)
  {
    messages <- append(messages, "'n': Value must be greater than 0.\n")
  }
  if(reals <= 0)
  {
    messages <- append(messages, "'reals': Value be greater than 0.\n")
  }
  if(cpus <= 0)
  {
    messages <- append(messages, "cpus must be greater than 0.\n")
  }
  if(!(is.null(observations)))
  {
    if(!(class(observations) == "data.frame"))
    {
      messages <- append(messages, "'observations': Not a valid data.frame.\n")
    }
  }
  if(!(file.exists(outputdir)))
  {
    messages <- append(messages, "'outputdir': Output directory does not exist.\n")
  }
  if(!(is.null(strata)))
  {
    if(!(class(strata) == "RasterLayer"))
    {
      messages <- append(messages, "'strata': Not a valid RasterLayer.\n")
    }
  }
  if(length(messages) > 1)
  {
    stop(messages)
  }
  
  # Strip trailing / of outputdir, if it exists
  if(substr(outputdir, nchar(outputdir), nchar(outputdir) + 1) == "/")
  {
    outputdir <- substr(outputdir, 1, nchar(outputdir) - 1)
  }
  
  # Set stub to "" if NULL
  if(is.null(stub))
  {
    stub <- ""
    
  } else if (stub == "") {
    
    stub <- ""
    
  } else if(!(substr(stub, nchar(stub), nchar(stub)) == "_")) {
    
    stub <- paste0(stub, "_")
  }
  
  # Create subdirectories to store results in
  dir.create(paste0(outputdir, "/output/"), showWarnings = FALSE)
  dir.create(paste0(outputdir, "/output/realisations"), showWarnings = FALSE)
  dir.create(paste0(outputdir, "/output/trees"), showWarnings = FALSE)
  
  # Generate lookup table
  if(!(is.null(strata))) {
    names(composition) <- c("poly", "mapunit", "stratum", "soil_class", "proportion")
  } else {
    names(composition) <- c("poly", "mapunit", "soil_class", "proportion")
  }
  
  # Write covariate names to file
  write.table(names(covariates), paste0(outputdir, "/output/", stub,
                                        "covariate_names.txt"),
              quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
  
  # Get samples for all realisations
  message(paste0("Generating samples for ", reals, " realisations"))
  samples <- data.frame()
  
  if(is.null(strata))
  {
    # Get samples without allocation stratification
    samples <- .getVirtualSamples(covariates, polygons, composition,
                                  n.realisations = reals, rate = rate,
                                  method.sample = method.sample,
                                  method.allocate = method.allocate)
  } else {
    samples <- .getStratifiedVirtualSamples(covariates, polygons, composition,
                                            strata, n.realisations = reals,
                                            rate = rate, 
                                            method.sample = method.sample, 
                                            method.allocate = method.allocate)
  }
  
  # If there are observations, get their covariates
  if(!(is.null(observations)))
  { 
    names(observations) <- c("x", "y", "class")
    observations <- .observations(observations, covariates)
    write.table(observations, paste0(outputdir, "/output/", stub,
                                     "observations_with_covariates.txt"),
                sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  # Write samples to text file
  write.table(samples, paste0(outputdir, "/output/", stub,
                              "virtual_samples.txt"),
              sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  # We submit the target classes to C5.0 as a factor. To do that, we need to 
  # make sure that the factor has the full set of levels, since due to the 
  # randomness of the allocation procedure we can't assume that they are all
  # represented in the samples drawn for a particular realisation. If the levels
  # are not specified correctly and they are not all represented in all sets of
  # samples, it is possible that the integer values that the class predictions
  # are coded to by raster::predict will not refer to the same soil class from
  # realisation to realisation.
  levs <- as.character(unique(composition$soil_class))
  
  if(!(is.null(observations)))
  {
    # Union map unit levels with observation levels
    levs <- base::union(levs, as.character(unique(observations$soil_class)))
  }
  
  # Sort levels
  levs <- sort(levs)
  
  # Generate lookup table
  lookup = as.data.frame(levs)
  lookup$code = seq(from=1, to=nrow(lookup), by=1)
  colnames(lookup) = c("name", "code")
  
  # Write lookup table to file
  write.table(lookup, paste0(outputdir, "/output/", stub,"lookup.txt"),
              sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  # Process realisations
  for (j in 1:reals)
  {
    message(paste0("Realisation ", j))
    
    # Column number of samples' first covariate
    startcol <- numeric(0)
    if(is.null(strata))
    {
      startcol = 8
    } else {
      startcol = 9
    }
    
    # Extract sample covariates for the current realisation
    s <- samples[which(samples$realisation == j), startcol:ncol(samples)]
    soil_class <- as.character(samples$soil_class[which(samples$realisation == j)])
    
    # Concatenate observations if they exist
    if(!(is.null(observations)))
    {
      s <- rbind(s, observations[, 8:ncol(observations)])
      soil_class <- append(soil_class, as.character(observations$soil_class))
    }
    
    # Sort levels and convert soil_class back to factor
    soil_class <- base::factor(soil_class, levels = levs)

    # Fit classification tree
    tree = C50::C5.0(s, y = soil_class)
    
    # Save tree to text file
    out <- utils::capture.output(summary(tree))
    cat(out, file = paste0(outputdir, "/output/trees/", stub, "tree_",
                           formatC(j, width = nchar(reals), format = "d",
                                   flag = "0"), ".txt"),
        sep = "\n", append = TRUE)
    
    # Save tree to rdata file
    save(tree, file = paste0(outputdir, "/output/trees/", stub, "tree_",
                             formatC(j, width = nchar(reals), format = "d",
                                     flag = "0"), ".RData"))
    
    # Predict realisation and save it to raster
    raster::beginCluster(cpus)
    r1 <- raster::clusterR(covariates, predict, args = list(tree),
                           filename = paste0(outputdir, "/output/realisations/",
                                             stub, "realisation_", formatC(j, width = nchar(reals), format = "d", flag = "0"), ".tif"),
                           format = "GTiff", overwrite = TRUE, datatype = "INT2S")
    raster::endCluster()
  }
}

#END