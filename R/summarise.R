#' Summarise the results of a DSMART spatial disaggregation.
#' 
#' \code{summarise} summarises the results of the spatial disaggregation of a 
#' polygon soil map in several ways. First, it computes the probabilities of 
#' occurrence of the soil classes that occur across all the map's map units. 
#' Second, it ranks the soil class predictions according to their probabilities 
#' of occurrence and maps the \emph{n}-most-probable soil classes at each grid 
#' cell, and their probabilities. Finally, it computes Shannon's entropy on the
#' class probabilities, and the degree of confusion between the most probable 
#' and second-most-probable soil classes.
#' 
#' @param realisations A \code{RasterStack} where each layer contains one 
#'   realisation of the soil class distribution across the soil map area, as 
#'   produced by \code{\link{disaggregate}}. If probabilistic predictions are
#'   used (\code{type = "prob"}), a list of RasterBrick objects with predicted
#'   class probabilities must be passed.
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
#' @param outputdir A character string that identifies the location of the main 
#'   output directory. The folder \code{output} and its subfolders will be 
#'   placed here. Default is the current working directory, \code{getwd()}.
#' @param stub \emph{optional} A character string that identifies a short name
#'   that will be prepended to all output.
#' @param prob A character vector to specify the type of the predictions to be 
#'   summarised. By default, raw class predictions are used. If set to "prob",
#'   probabilistic predictions are used.
#'
#' @return A list that contains metadata about the current run of
#'   \code{summarise}.
#'   
#' @references McBratney, A.B., Mendonca Santos, M. de L., Minasny, B., 2003. On
#'   digital soil mapping. Geoderma 117, 3--52. doi: 
#'   \href{https://doi.org/10.1016/S0016-7061(03)00223-4}{10.1016/S0016-7061(03)00223-4}
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
#'   \href{https://doi.org/10.1016/j.geoderma.2013.09.024}{10.1016/j.geoderma.2013.09.024}
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

summarise <- function(realisations, lookup, 
                      n.realisations = ifelse(is.list(realisations), length(realisations), 
                                              terra::nlyr(realisations)),
                      nprob = 3, outputdir = getwd(), stub = NULL, type = "response")
{
  # TO DELETE ONCE PACKAGE IS WRITTEN:
  # Rcpp::sourceCpp("./src/order.cpp")
  # Rcpp::sourceCpp("./src/sort.cpp")
  
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
  ordered.indices.factors <- rast(lapply(1:nprob, function(x) {
    
    # Subset from stack
    ordered.indices.factor <- as.factor(subset(ordered.indices, x))
    
    # Get factor labels and codes
    lookup_labels <- data.frame(levels = levels(ordered.indices.factor)[[1]]$levels,
                                code = levels(ordered.indices.factor)[[1]]$labels)
    
    # Merge with lookup table and rearrange
    label_lookup <- lookup %>% 
      merge(lookup_labels) %>% 
      dplyr::select(levels, name, code)
    
    # Write full index to file
    write.csv(label_lookup, file.path(
      outputdir, "output", "mostprobable", 
      paste0(stub, "mostprob_",
             formatC(x, width = nchar(nrow(lookup)), format = "d", flag = "0"),
             "_class_lut.csv")),
      row.names = FALSE)
    
    # Apply factor levels and names to raster
    levels(ordered.indices.factor) <- label_lookup[, c(1, 2)]
    
    # Write out factored raster
    ordered.indices.factor <- terra::writeRaster(
      ordered.indices.factor, filename = file.path(
        outputdir, "output", "mostprobable", 
        paste0(stub, "mostprob_",
               formatC(x, width = nchar(nrow(lookup)), format = "d", flag = "0"),
               "_class.tif")), 
      overwrite = TRUE, wopt = list(datatype = "INT1U"))
    
    # Write associated probability file
    ordered.probs.subset <- subset(ordered.probs, x) %>% 
      writeRaster(filename = file.path(outputdir, "output", "mostprobable",
                                       paste0(stub, "mostprob_",
                                              formatC(x, width = nchar(nrow(lookup)), format = "d", flag = "0"),
                                              "_probs.tif")),
                  overwrite = TRUE)
    
    return(ordered.indices.factor)
  }))

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