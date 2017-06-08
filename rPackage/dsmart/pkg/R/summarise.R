#' Summarise the results of a DSMART spatial disaggregation.
#' 
#' \code{summarise} summarises the results of the spatial disaggregation of a 
#' polygon soil map in several ways. First, it computes the probabilities of 
#' occurrence of the soil classes that occur across all the map's map units. 
#' Second, it ranks the soil class predictions according to their probabilities 
#' of occurrence and maps the \emph{n}-most-probable soil classes at each grid 
#' cell, and their probabilities. Finally, it computes the degree of confusion 
#' between the most probable and second-most-probable soil classes.
#' 
#' @param realisations A \code{RasterStack} where each layer contains one 
#'   realisation of the soil class distribution across the soil map area, as 
#'   produced by \code{\link{disaggregate}}.
#' @param lookup A two-column \code{data.frame} containing a mapping between the
#'   integer soil class codes in the layers of \code{realisations}, and the soil
#'   class codes defined by the map unit composition \code{data.frame} used as
#'   an argument to \code{disaggregate} and \code{dsmart}. \code{lookup} is the
#'   same lookup table that is produced by \code{dsmart}. First column is the
#'   soil class code of the map unit composition; second column is the integer
#'   soil class code.
#' @param n.realisations An integer that identifies the number of realisations
#'   of the soil class distribution that were computed by \code{disaggregate}.
#'   Default value is \code{raster::nlayers(realisations)}.
#' @param nprob At any location, disaggregated soil class predictions can be 
#'   ranked according to their probabilities of occurence. \code{rdsmart} can 
#'   map the class predictions, and their probabilities, at any rank. 
#'   \code{nprob} is an integer that identifies the number of probability ranks 
#'   to map. For example, if \code{n = 3}, DSMART will map the first-, second- 
#'   and third-most-probable soil classes and their probabilities of occurrence.
#' @param cpus An integer that identifies the number of CPU processors to use 
#'   for parallel processing.
#' @param outputdir A character string that identifies the location of the main 
#'   output directory. The folder \code{output} and its subfolders will be 
#'   placed here. Default is the current working directory, \code{getwd()}.
#' @param stub \emph{optional} A character string that identifies a short name
#'   that will be prepended to all output.
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
#' @examples
#' # Load datasets
#' data(dalrymple_lookup)
#' data(dalrymple_realisations)
#' 
#' # Summarise
#' summarise(dalrymple_realisations, dalrymple_lookup, nprob = 5, cpus = 6)
#' 
#' @export

summarise <- function(realisations, lookup, n.realisations = raster::nlayers(realisations),
                      nprob = 3, cpus = 1, outputdir = getwd(), stub = NULL)
{
  # Check arguments before proceeding
  messages <- c("Attention is required with the following arguments:\n")
  if(!(class(realisations) == "RasterStack"))
  {
    messages <- append(messages, "'realisations': Not a valid RasterStack.\n")
  }
  if(!(class(lookup) == "data.frame"))
  {
    messages <- append(messages, "'lookup': Not a valid data.frame.\n")
  }
  if(n.realisations <= 0)
  {
    messages <- append(messages, "'n.realisations': Value must be greater than 0.\n")
  }
  if(nprob <= 0)
  {
    messages <- append(messages, "'nprob': Value must be greater than 0.\n")
  }
  if(cpus <= 0)
  {
    messages <- append(messages, "'cpus': Value must be greater than 0.\n")
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
    
  } else if (stub == "") {
    
    stub <- ""
    
  } else if(!(substr(stub, nchar(stub), nchar(stub)) == "_")) {
    
    stub <- paste0(stub, "_")
  }
  
  # Set up output directories
  dir.create(paste0(outputdir, "/output/probabilities/"), showWarnings = FALSE)
  dir.create(paste0(outputdir, "/output/mostprobable/"), showWarnings = FALSE)
  
  # Parameter to pass to counts function as a global variable
  param <- nrow(lookup)
  assign("param", param, envir = .GlobalEnv)
  
  # Compute counts
  raster::beginCluster(cpus)
  counts <- raster::clusterR(realisations, calc,
                             args = list(fun = function(x) {tabulate(x, nbins = param)}),
                             export = "param")
  raster::endCluster()
  
  # Parameter to pass to probabilities function as a global variable
  assign("n.realisations", n.realisations, envir = .GlobalEnv)
  
  # Compute probabilities
  # probs = counts / n.realisations is faster on small datasets.
  raster::beginCluster(cpus)
  probs <- raster::clusterR(counts, calc,
                            args = list(fun = function(x) {x / n.realisations}),
                            export = "n.realisations")
  raster::endCluster()
  
  # Write probabilities to raster files
  for(i in 1:raster::nlayers(probs))
  {
    raster::writeRaster((probs[[i]]),
                        filename = paste0(outputdir, "/output/probabilities/", stub, "prob_", lookup$name[which(lookup$code == i)], ".tif"),
                        format = "GTiff", overwrite = TRUE)
  }
  
  # Compute the class indices of the n-most-probable soil classes
  raster::beginCluster(cpus)
  ordered.indices = raster::clusterR(counts, calc, 
                                     args = list(fun = function(x) order(x, decreasing = TRUE, na.last = TRUE)))
  raster::endCluster()
  
  # Compute the class probabilities of the n-most-probable soil classes
  raster::beginCluster(cpus)
  ordered.probs = raster::clusterR(probs, calc, 
                                   args = list(fun = function(x) sort(x, decreasing = TRUE, na.last = TRUE)))
  raster::endCluster()
  
  for (i in 1:nprob)
  {
    # Write ith-most-probable soil class raster to file
    raster::writeRaster(ordered.indices[[i]],
                        filename = paste0(outputdir, "/output/mostprobable/", stub, "mostprob_",
                                          formatC(i, width = nchar(nrow(lookup)), format = "d", flag = "0"), "_class.tif"),
                        format = "GTiff", overwrite = TRUE)
    
    # Write ith-most-probable soil class probability raster to file
    raster::writeRaster(ordered.probs[[i]],
                        filename = paste0(outputdir, "/output/mostprobable/", stub, "mostprob_",
                                          formatC(i, width = nchar(nrow(lookup)), format = "d", flag = "0"), "_probs.tif"),
                        format = "GTiff", overwrite = TRUE)
  }

  # Compute the confusion index between the most probable and second-most-probable soil classes
  # and write it to file
  raster::beginCluster(cpus)
  confusion <- raster::clusterR(ordered.probs, fun = function(x) (1 - (x[[1]] - x[[2]])),
                                filename = paste0(outputdir, "/output/mostprobable/", stub, "confusion.tif"), format = "GTiff", overwrite = TRUE)
  raster::endCluster()
}

#END