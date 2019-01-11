#' Sample the polygons of a SpatialPolygonsDataFrame.
#' 
#' \code{.getVirtualSamples} samples the polygons of a SpatialPolygonsDataFrame.
#' It accomplishes three tasks: (i) it draws samples from within each polygon, 
#' (ii) extracts the values of the covariates at the sample locations, and (iii)
#' through \code{.allocate}, allocates each sample to a soil class. See 
#' \emph{Details} for more information.
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
#' @param n.realisations An integer that identifies the number of realisations 
#'   of the soil class distribution that DSMART should compute.
#' @param rate An integer that identifies the number of virtual samples to draw 
#'   from each polygon in each realisation. If \code{method.sample =
#'   "by_polygon"}, the number of samples to draw from each polygon in
#'   \code{polygons}. If \code{method.sample = "by_area"}, the sampling density
#'   in number of samples per square kilometre.
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
#' @param cpus An integer that identifies the number of CPU processors to use 
#'   for parallel processing.
#' 
.getVirtualSamples <- function(covariates, polygons, composition,
                               n.realisations = 100, rate = 15,
                               method.sample = "by_polygon", 
                               method.allocate = "weighted", cpus = 1)
{
  # Initialise cluster
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)
  
  samples <- foreach::foreach(poly.id = polygons@data[, 1],
                               .packages = c('raster', 'sp'),
                              .export = ".allocate") %dopar% #It was necessary to export the .allocate function when testing on my own PC, but it may have to be deleted later.
  {
    # Subset a polygon
    poly <- base::subset(polygons, polygons@data[, 1] == poly.id)
    
    # If sample = "area", determine the correct number of samples to take
    n.samples <- 0
    if(method.sample == "by_area") {
    
      # Compute area of polygon in square kilometres
      # FOR NOW, assumes that CRS of polygons is projected and with units of m.
      # The area function of the raster package also gives the area in square
      # metres for longitude/latitude coordinate systems.
      area <- raster::area(poly) / 1e6
      
      if(area < 1.0)
      {
        area <- 1.0
      }
      
      # Compute number of samples to take
      n.samples <- rate * base::trunc(area)
    } else if(method.sample == "by_polygon") {
      n.samples <- rate
    } else stop("Sampling method \"", method.sample, "\" is unknown")
    
    # Extract covariates of all grid cells that intersect with the polygon
    # Retain only those cells that do not have NA in their covariates
    poly.cells <- as.data.frame(raster::extract(covariates, poly, 
                                                cellnumbers = TRUE))
    poly.cells <- poly.cells[base::which(stats::complete.cases(poly.cells)), ]
    
    # Sample grid cells with replacement for ALL realisations in one go
    poly.samples <- poly.cells[base::sample(base::nrow(poly.cells),
                                            replace = TRUE,
                                            size = n.samples * n.realisations), ]
    
    # Allocate all samples to a soil class
    soil_class <- character()
    if(method.allocate == "weighted") {
      # Weighted random allocation
      poly.classes <- base::as.character(composition[base::which(composition[, 1] == poly.id), 3])
      poly.weights <- composition[base::which(composition[, 1] == poly.id), 4]
      
      if(length(poly.classes) == 0) {
        stop(paste0("No map unit composition for polygon ", poly.id))
      } else {
        soil_class <- .allocate(poly.classes, n = n.samples * n.realisations, 
                                method = "weighted", weights = poly.weights)
      }
    } else if(method.allocate == "random-mapunit") {
      # Random-within_map unit allocation
      poly.classes <- as.character(composition[which(composition[, 1] == poly.id), 3])
      
      if(length(poly.classes) == 0) {
        stop(paste0("No map unit composition for polygon ", poly.id))
      } else {
        soil_class <- .allocate(poly.classes, n = n.samples * n.realisations, 
                                method = "random", weights = NULL)
      }
    } else if(method.allocate == "random-all") {
      # Random-within_map allocation
      poly.classes <- as.character(unique(composition[, 3]))
      
      if(length(poly.classes) == 0) {
        stop("No soil classes available to allocate to")
      } else {
        soil_class <- .allocate(poly.classes, n = n.samples * n.realisations, 
                                method = "random", weights = NULL)
      }
    } else stop("Allocation method is unknown")
    
    # Add realisation id, spatial coordinates, soil class to sampled grid cells
    xy <- as.data.frame(raster::xyFromCell(covariates, poly.samples$cell))
    meta <- list(realisation = base::rep(1:n.realisations, times = n.samples),
                 type = base::rep("virtual", nrow(xy)),
                 sampling = base::rep(method.sample, base::nrow(xy)),
                 allocation = base::rep(method.allocate, base::nrow(xy)))
    poly.samples <- cbind(as.data.frame(meta), xy, soil_class, 
                          poly.samples[, 2:base::ncol(poly.samples)])
    
    # Add polygon samples to list
    #samples <- append(samples, list(poly.samples))
    return(poly.samples)
  }
  
  parallel::stopCluster(cl)
  
  # Merge polygon sample data frames
  samples <- data.table::rbindlist(samples)
  
  return(samples)
}

#'Sample the polygons of a SpatialPolygonsDataFrame.
#'
#'
.getStratifiedVirtualSamples <- function(covariates, polygons, composition, strata,
                                     n.realisations = 100, rate = 15,
                                     method.sample = "by_polygon", 
                                     method.allocate = "weighted")
{
  # Make sure composition column names are formatted properly
  if(ncol(composition) == 4) {
    names(composition) <- c("poly", "mapunit", "soil_class", "proportion")
  } else if (ncol(composition) == 5) {
    names(composition) <- c("poly", "mapunit", "stratum", "soil_class", "proportion")
  } else stop("Map unit composition in unknown format.")
  
  # Empty list to hold samples
  samples <- list()
  
  # Process each polygon in polygons
  for(poly.id in polygons@data[, 1])
  {
    # Subset a polygon
    poly <- subset(polygons, polygons@data[, 1] == poly.id)
    
    # If sample = "area", determine the correct number of samples to take
    n.samples <- 0
    if(method.sample == "by_area") {
      # Compute area of polygon in square kilometres
      # FOR NOW, assumes that CRS of polygons is projected and with units of m.
      # The area function of the raster package also gives the area in square
      # metres for longitude/latitude coordinate systems.
      area <- raster::area(poly) / 1e6
      
      if(area < 1.0)
      {
        area <- 1.0
      }
      
      # Compute number of samples to take
      n.samples <- rate * base::trunc(area)
    } else if(method.sample == "by_polygon") {
      n.samples <- rate
    } else stop("Sampling method \"", method.sample, "\" is unknown")
    
    # Extract covariates of all grid cells that intersect with the polygon
    # Retain only those cells that do not have NA in their covariates
    poly.cells <- as.data.frame(raster::extract(covariates, poly, cellnumbers = TRUE))
    poly.cells <- poly.cells[which(complete.cases(poly.cells)), ]
    
    # Sample grid cells with replacement for ALL realisations in one go
    poly.samples <- poly.cells[sample(nrow(poly.cells), replace = TRUE, size = n.samples * n.realisations), ]
    
    # Get coordinates of sampled grid cells
    xy <- as.data.frame(raster::xyFromCell(covariates, poly.samples$cell))
    
    # Determine what strata they belong to
    poly.samples.strata <- extract(strata, xy)
    
    # Allocate soils within strata
    poly.samples <- cbind(poly.samples.strata, poly.samples[, 2:ncol(poly.samples)])
    names(poly.samples)[names(poly.samples) == "poly.samples.strata"] <- "stratum"
    poly.samples <- poly.samples[order(poly.samples$stratum), ]
    
    soil_class <- character()
    for(stratum in unique(poly.samples[, 1]))
    {
      # Work out the number of samples to allocate
      stratum.n <- length(which(poly.samples[, 1] == stratum))
      
      # Candidate classes
      stratum.classes <- as.character(composition[which((composition$poly == poly.id) & (composition$stratum == stratum)), 4])
      
      if(length(stratum.classes) == 0) {
        
        stop(paste0("No classes available in stratum ", stratum, " of polygon ", poly.id))
        
      } else if(method.allocate == "weighted") {
        
        # Weighted-random allocation within stratum
        if(!(ncol(composition) == 5)) {
          
        stop("Weighted-random allocation specified but no stratum weights available.")
          
        } else {
        # Class weights (proportions)
        stratum.weights <- composition[which((composition$poly == poly.id) &
                                               (composition$stratum == stratum)), 5]
          
        # Perform allocation
        alloc <- .allocate(stratum.classes, n = stratum.n, method = "weighted",
                           weights = stratum.weights)
          
        soil_class <- append(soil_class, alloc)
          
        } 
      } else if(method.allocate == "random") {
          
        # Completely random allocation within stratum
          
        # Perform allocation
        alloc <- .allocate(stratum.classes, n = stratum.n, method = "random")
        soil_class <- append(soil_class, alloc)
      }
    }
    
    # Attach soil classes to samples
    poly.samples <- cbind(xy, soil_class, poly.samples[, 1:ncol(poly.samples)])
    
    # Shuffle samples row-wise so that they're not in order of stratum
    poly.samples <- poly.samples[base::sample(nrow(poly.samples)), ]
    
    # Add realisation id, spatial coordinates, soil class to sampled grid cells
    meta <- list(realisation = rep(1:n.realisations, times = n.samples),
                 type = base::rep("virtual", nrow(xy)),
                 sampling = base::rep(method.sample, nrow(xy)),
                 allocation = base::rep(method.allocate, nrow(xy)))
    poly.samples <- cbind(as.data.frame(meta), poly.samples)
    
    # Add polygon samples to list
    samples <- base::append(samples, list(poly.samples))
    
  }
  
  # Merge polygon sample data frames
  samples <- data.table::rbindlist(samples)
  
  return(samples)
}

#' Allocate samples to soil classes
#' 
#' \code{.allocate} randomly or weighted-randomly draws a sample of size
#' \code{n} from \code{classes}.
#' 
#' @param classes The classes to allocate to.
#' @param n The number of samples to allocate.
#' @param method The method of allocation. Valid values are \code{"weighted"},
#'   for weighted-random allocation using the weights in \code{weights}, and
#'   \code{"random"} for random allocation (the default).
#'
.allocate <- function(classes, n = 15, method = "random", weights = NULL)
{
  # Check parameter values before proceeding
  if((base::length(classes) == 0) | (base::is.null(classes))) {
    
    stop("Classes are not specified.")
    
  }
  
  if(n < 1) {
    
    stop("n must be greater than 0.")
    
  }
  
  if (method == "weighted") {
    
    if(base::is.null(weights)) {
      
      stop("Weighted random allocation specified but no weights supplied.")
      
    } else if(!(base::length(classes) == base::length(weights))) {
      
      stop("Number of classes is not the same as number of weights.")
      
    }
  }

  ### This is useful functionality but is too slow  
  # if(base::all(base::is.na(classes)) == TRUE) {
  #   
  #   # Check whether all values in classes are NA, in the event that a map
  #   # unit's composition is present in the composition file but undefined.
  #   # Return a vector of NAs of length n (deal with its implications
  #   # appropriately in .getVirtualSamples or .getStratifiedVirtualSamples)
  #   return(base::rep(NA, times = n))
  # }
  
  # Check whether all values in classes or weights are NA, in the event that a
  # map unit's composition is present but undefined
  if(base::all(base::is.na(classes)) == TRUE) {
    
    # Return a vector of NAs of length n (deal with it appropriately in
    # .getVirtualSamples or .getStratifiedVirtualSamples)
    return(base::rep(NA, times = n))
  }
  
  # Perform allocation
  allocation <- character()
  if(method == "weighted"){
    
    # Weighted-random allocation
    # Draw from Dirichlet distribution
    s <- gtools::rdirichlet(1, weights)
    
    # Make prediction
    allocation = base::sample(classes, size = n, replace = TRUE, prob = s[1, ])
  }
  else if(method == "random")
  {
    # Completely random allocation
    # Make prediction
    allocation = base::sample(classes, size = n, replace = TRUE)
  }
  else stop("Allocation method \"", method, "\" is unknown")
  
  return(allocation)
}


.observations <- function(observations, covariates)
{
  # Rename columns for consistency with virtual samples
  base::colnames(observations) <- c("x", "y", "soil_class")
  
  # Identify coordinate fields of observations data frame
  sp::coordinates(observations) <- c(1, 2)
  
  # Extract covariates at observation locations
  o.covariates <- raster::extract(covariates, observations)
  
  # Join covariates back to observations
  meta <- list(realisation = numeric(length = nrow(observations)),
               type = base::rep("actual", nrow(observations)),
               sampling = base::rep("observed", nrow(observations)),
               allocation = base::rep("observed", nrow(observations)))
  obs <- cbind(as.data.frame(meta), as.data.frame(observations), o.covariates)
  
  # Return only those observations with no NA in the covariates
  obs <- obs[which(complete.cases(obs)), ]
  return(obs)
}