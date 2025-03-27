#' Summarise the results of a DSMART spatial disaggregation.
#' 
#' `summarise` summarises the results of the spatial disaggregation of a 
#' polygon soil map in several ways. First, it computes the probabilities of 
#' occurrence of the soil classes that occur across all the map's map units. 
#' Second, it ranks the soil class predictions according to their probabilities 
#' of occurrence and maps the *n*-most-probable soil classes at each grid 
#' cell, and their probabilities. Finally, it computes Shannon's entropy on the
#' class probabilities, and the degree of confusion between the most probable 
#' and second-most-probable soil classes.
#' 
#' @param realisations A SpatRaster where each layer contains one 
#'   realisation of the soil class distribution across the soil map area, as 
#'   produced by [disaggregate()]. If probabilistic predictions are
#'   used (`type = "prob"`), a list of RasterBrick objects with predicted
#'   class probabilities must be passed.
#' @param lookup A two-column `data.frame` containing a mapping between the
#'   integer soil class codes in the layers of `realisations`, and the soil
#'   class codes defined by the map unit composition `data.frame` used as
#'   an argument to `disaggregate` and `dsmart`. `lookup` is the
#'   same lookup table that is produced by `dsmart`. First column is the
#'   soil class code of the map unit composition; second column is the integer
#'   soil class code.
#' @param n.realisations An integer that identifies the number of realisations
#'   of the soil class distribution that were computed by `disaggregate`.
#'   Default value is `terra::nlyr(realisations)`.
#' @param nprob At any location, disaggregated soil class predictions can be 
#'   ranked according to their probabilities of occurence. `rdsmart` can 
#'   map the class predictions, and their probabilities, at any rank. 
#'   `nprob` is an integer that identifies the number of probability ranks 
#'   to map. For example, if `n = 3`, DSMART will map the first-, second- 
#'   and third-most-probable soil classes and their probabilities of occurrence.
#' @param outputdir A character string that identifies the location of the main 
#'   output directory. The folder `output` and its subfolders will be 
#'   placed here. Default is the current working directory, `getwd()`.
#' @param stub *optional* A character string that identifies a short name
#'   that will be prepended to all output.
#' @param type A character vector to specify the type of the predictions to be 
#'   summarised. By default, "response" class predictions are used. If set to 
#'   "prob", probabilistic predictions are used.
#'
#' @return A list that contains metadata about the current run of
#'   `summarise`.
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
#' data(dalrymple_lookup)
#' data(dalrymple_realisations)
#' 
#' # Summarise
#' summarise(dalrymple_realisations, dalrymple_lookup, nprob = 5)
#' 
#' @export

summarise <- function(realisations, lookup, 
                      n.realisations = ifelse(is.list(realisations), length(realisations), 
                                              terra::nlyr(realisations)),
                      nprob = 3, outputdir = getwd(), stub = NULL, type = "response")
{
  
  # Create list to store output
  output <- base::list()
  
  # Save start time
  output$timing <- base::list(start = Sys.time())
  
  # Check arguments before proceeding
  messages <- c("Attention is required with the following arguments:\n")
  if(type != "prob")
  {
    if(!(class(realisations) == "SpatRaster"))
    {
      messages <- append(messages, "'realisations': Not a valid SpatRaster\n")
    }
  }else{
    if(is.list(realisations) == FALSE) {
      messages <- append(messages, "'realisations' must be a list of SpatRaster objects when probabilistic predictions are used.'.\n")
    }else{
      if(sum(unlist(lapply(realisations, function(x) class(x) != "SpatRaster"))) > 0) {
        messages <- append(messages, "'realisations' must be a list of SpatRaster objects when probabilistic predictions are used.'.\n")
      }
    }
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
  if(!(file.exists(outputdir)))
  {
    messages <- append(messages, "'outputdir': Output directory does not exist.")
  }
  if(length(messages) > 1)
  {
    stop(messages)
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
  
  # Save function call
  output$call <- base::match.call()
  
  # Save parameters
  output$parameters <- base::list(n.realisations = n.realisations,
                                  nprob = nprob, stub = stub, type = type)
  
  # Set up output directories
  outputdir <- file.path(outputdir)
  dir.create(file.path(outputdir, "output"), showWarnings = FALSE)
  dir.create(file.path(outputdir, "output", "probabilities"), showWarnings = FALSE)
  dir.create(file.path(outputdir, "output", "mostprobable"), showWarnings = FALSE)
  
  # Save output locations
  output$locations <- base::list(root = file.path(outputdir, "output"),
                                 probabilities = file.path(outputdir, "output", "probabilities"),
                                 mostprobable = file.path(outputdir, "output", "mostprobable"))
  
  # Make sure lookup table column names are correct
  names(lookup) <- c("name", "code")
  
  # If raw predictions are used, calculate class probabilities by counting.
  if(type != "prob"){
    # Compute counts
    counts <- terra::app(realisations, function(x) {
      if(is.na(sum(x))) {
        rep(NA, nrow(lookup))
      } else {
        tabulate(x, nbins = nrow(lookup))
      }
    }, wopt = list(names = lookup$name))
    
    # Compute probabilities and write probabilities to raster files
    probs <- (counts / n.realisations) %>% 
      terra::writeRaster(filename = file.path(outputdir, "output", "probabilities",
                                              paste0(stub, "prob_", names(.), ".tif")), 
        overwrite = TRUE)
  }else{
    # If probabilistic predictions are used, calculate class probabilities by averaging
    # the predicted probabilities across the realisations.
    # If only one realisation is used, no averaging is needed.
    if(length(realisations) == 1 | n.realisations == 1){
      probs <- terra::writeRaster(realisations[[1]], 
                                  filename = file.path(outputdir, "output", "probabilities", 
                                                       paste0(stub, "prob_", names(realisations[[1]]), ".tif")), 
                                  overwrite = TRUE)
    }else{
      probs <- rast(lapply(1:nrow(lookup), function(i) {
        rlist <- mean(rast(lapply(realisations, "[[", i)), na.rm = TRUE) %>% 
          terra::writeRaster(file.path(
            outputdir, "output", "probabilities", 
            paste0(stub, "prob_", lookup$name[which(lookup$code == i)], ".tif")),
            overwrite = TRUE, wopt = list(names = lookup$name[which(lookup$code == i)]))
      }))
    }
  }
  
  # Compute the class indices of the n-most-probable soil classes
  # assign("nprob", nprob, envir = .GlobalEnv)
  if(type != "prob")
  {
    # If raw class predictions are used, use "counts" for indicing.
    ordered.indices <- order_stack_values(counts, n = nprob)
    
  }else{
    # If probabilistic predictions are used, use "probs" for indicing.
    ordered.indices <- order_stack_values(probs, n = nprob)
  }
  
  # Compute the class probabilities of the n-most-probable soil classes
  ordered.probs <- sort_stack_values(probs, n = max(2, nprob))
  
  # Write nprob soil class rasters to files as factors
  levels(ordered.indices) <- lapply(1:nprob, function(x) 
    lookup %>% dplyr::select(code, name))
  
  ordered.ind.names <- paste0(
    stub, "mostprob_", formatC(1:nprob, width = nchar(nrow(lookup)), 
                               format = "d", flag = "0"), "_class")
  
  ordered.indices <- terra::writeRaster(
    ordered.indices, filename = file.path(
      outputdir, "output", "mostprobable",
      paste0(ordered.ind.names, ".tif")),
    overwrite = TRUE, wopt = list(datatype = "INT1U", names = ordered.ind.names))
  
  # Write most probable layers
  ordered.prob.names <- paste0(
    stub, "mostprob_", formatC(1:nprob, width = nchar(nrow(lookup)), 
                               format = "d", flag = "0"), "_probs")
  
  ordered.probs.out <- subset(ordered.probs, 1:nprob) %>% 
    writeRaster(filename = file.path(
      outputdir, "output", "mostprobable",
      paste0(ordered.prob.names, ".tif")),
      overwrite = TRUE, wopt = list(names = ordered.prob.names))

  # Compute the confusion index
  confusion <- confusion_index(ordered.probs) %>% 
    writeRaster(filename = file.path(outputdir, "output", "mostprobable",
                                     paste0(stub, "confusion.tif")),
                overwrite = TRUE)
  
  # Compute Shannon's entropy on the class probabilities
  shannon <- shannon_entropy(ordered.probs) %>% 
    writeRaster(filename = file.path(outputdir, "output", "mostprobable",
                                     paste0(stub, "shannon.tif")),
                overwrite = TRUE)
  
  # Save finish time
  output$timing$finish <- Sys.time()
  
  # Return output
  return(output)
}

#END
