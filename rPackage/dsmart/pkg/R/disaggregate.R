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
#'   soil map unit that the polygon belongs to. Third column contains a code 
#'   that identifies the soil class. Fourth column contains a number in the 
#'   range \code{(0, 100)} that identifies the proportion of the map unit that 
#'   the soil class corresponds to.
#' @param n An integer that identifies the number of virtual samples to draw 
#'   from each polygon in each realisation.
#' @param reals An integer that identifies the number of realisations of the 
#'   soil class distribution that DSMART should compute.
#' @param observations \emph{Optional} A \code{data.frame} that contains actual 
#'   observations of the soil class at locations across the soil map area. These
#'   data augment the virtual samples and are used in each realisation. Each row
#'   contains information about one soil class observation. First and second 
#'   fields contain the \emph{x-} and \emph{y-}components of the observation's 
#'   spatial location. Third field is the soil class code. See \emph{Details}.
#' @param allocate Method of allocation of virtual samples to soil classes.
#'   Valid values are \code{"weighted"}, for weighted-random allocation to a
#'   soil class from within the virtual sample's map unit; 
#'   \code{"random-mapunit"}, for completely random allocation to a soil class
#'   from within the virtual sample's map unit; and \code{"random-all"}, for
#'   completely random allocation to a soil class from within the entire map 
#'   area.
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
#'  n = 15, reals = 10, cpus = 6)
#' 
#' # Run disaggregate with extra observations
#' disaggregate(dalrymple_covariates, dalrymple_polygons, dalrymple_composition,
#'  observations = dalrymple_observations, n = 15, reals = 10)
#' 
#' @references McBratney, A.B., Mendonca Santos, M. de L., Minasny, B., 2003. On
#'   digital soil mapping. Geoderma 117, 3--52. doi: 
#'   \href{http://dx.doi.org/10.1016/S0016-7061(03)00223-4}{10.1016/S0016-7061(03)00223-4}
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
#' @export

disaggregate <- function(covariates, polygons, composition, n = 15,
                         reals = 100, observations = NULL,
                         allocate = "weighted", outputdir = getwd(),
                         stub = NULL, cpus = 1)
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
  if(n <= 0)
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
    messages <- append(messages, "'outputdir': Output directory does not exist.")
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
  }
  else
  {
    stub <- paste0(stub, "_")
  }
  
  # Create subdirectories to store results in
  dir.create(paste0(outputdir, "/output/"), showWarnings = FALSE)
  dir.create(paste0(outputdir, "/output/realisations"), showWarnings = FALSE)
  dir.create(paste0(outputdir, "/output/trees"), showWarnings = FALSE)
  
  # Generate lookup table
  names(composition) <- c("poly", "mapunit", "soil_class", "proportion")
  lookup = as.data.frame(sort(unique(composition$soil_class)))
  lookup$code = seq(from=1, to=nrow(lookup), by=1)
  colnames(lookup) = c("name", "code")
  
  write.table(lookup, paste0(outputdir, "/output/", stub,"lookup.txt"),
              sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  # Get samples for all realisations
  message(paste0("Generating samples for ", reals, " realisations"))
  samples <- .getVirtualSamples(covariates, polygons, composition, n.realisations = reals,
                                n.samples = n, method = allocate)
  
  # If there are observations, get their covariates
  if(!(is.null(observations)))
  { 
    names(observations) <- c("x", "y", "class")
    observations <- .observations(observations, covariates)
  }
  
  # Write samples to text file
  write.table(rbind(samples, observations), paste0(outputdir, "/output/", stub, "samples.txt"),
              sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  # Process realisations
  for (j in 1:reals)
  {
    message(paste0("Realisation ", j))
    
    # Extract samples for the current realisation
    s <- samples[which(samples$r == j), ]
    
    # Concatenate observations if they exist
    if(!(is.null(observations)))
    {
      s <- rbind(s, observations)
    }
    
    # Fit classification tree
    tree = C50::C5.0(s[, 7:ncol(s)], y = s$soil_class)
    
    # Save tree to text file
    out <- utils::capture.output(summary(tree))
    cat(out, file = paste0(outputdir, "/output/trees/", stub, "tree_", formatC(j, width = nchar(reals), format = "d", flag = "0"), ".txt"),
        sep = "\n", append = TRUE)
    
    # Save tree to rdata file
    save(tree, file = paste0(outputdir, "/output/trees/", stub, "tree_", formatC(j, width = nchar(reals), format = "d", flag = "0"), ".RData"))
    
    # Predict realisation and save it to raster
    raster::beginCluster(cpus)
    r1 <- raster::clusterR(covariates, predict, args = list(tree),
                           filename = paste0(outputdir, "/output/realisations/", stub, "r_",
                                             formatC(j, width = nchar(reals), format = "d", flag = "0"), ".tif"),
                           format = "GTiff", overwrite = TRUE, datatype = "INT2S")
    raster::endCluster()
  }
}

#END