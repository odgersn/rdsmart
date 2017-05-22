#' Get samples
#' 
#' 
#'
.getVirtualSamples <- function(covariates, polygons, composition, n.realisations = 100, n.samples = 15, method = "weighted")
{
  # Empty data frame to hold samples
  samples <- data.frame()
  
  # Process each polygon in polygons
  for(poly.id in polygons@data[, 1])
  {
    # Subset a polygon
    poly <- subset(polygons, polygons@data[, 1] == poly.id)
    
    # Extract covariates of all grid cells that intersect with the polygon
    # Retain only those cells that do not have NA in their covariates
    poly.cells <- as.data.frame(raster::extract(covariates, poly, cellnumbers = TRUE))
    poly.cells <- poly.cells[which(complete.cases(poly.cells)), ]
    
    # Sample grid cells with replacement for ALL realisations in one go
    poly.samples <- poly.cells[sample(nrow(poly.cells), replace = TRUE, size = n.samples * n.realisations), ]
    
    # Allocate all samples to a soil class
    soil_class <- .allocate(composition, poly.id, n = n.samples * n.realisations, method = method)
    
    # Add realisation id, spatial coordinates, soil class to sampled grid cells
    xy <- as.data.frame(raster::xyFromCell(covariates, poly.samples$cell))
    meta <- list(r = rep(1:n.realisations, times = n.samples),
                 type = base::rep("virtual", nrow(xy)),
                 method = base::rep(method, nrow(xy)))
    poly.samples <- cbind(as.data.frame(meta), xy, soil_class, poly.samples[, 2:ncol(poly.samples)])
    
    # Add polygon samples to master data frame
    samples = rbind(samples, poly.samples)
  }
  
  # Return master data frame of samples
  return(samples)
}


#' Allocate samples to soil classes
#' 
#' \code{allocate} generates a character vector of soil class codes of length \code{n}. The soil classes belong to the soil map polygon
#' with id \code{polygon} as described in \code{composition}. The soil classes can be sampled by several different methods, as described in
#' \emph{Details}.
#' 
#' The default allocation method is \code{method = "weighted"}, which weighted-randomly samples the soil classes of a given polygon
#' using the soil classes' proportions of occurrence as weights. Thus if three soil classes A, B and C are in a soil map unit in
#' proportions of 70\%, 20\% and 10\% respectively, at each draw, class A would have a 70\% chance of being selected at random.
#' 
#' With \code{method = "random-mapunit"}, classes are drawn at complete random from among those that exist in the given polygon. This
#' effectively disregards the proportion of the map unit that each class is assumed to occur in. \code{method = "random-all"} works
#' the same way except the classes are drawn from among all those in the whole soil map, regardless of whether they occur in the given polygon
#' or not.
#' 
#' @param composition A \code{data.frame} containing the soil class composition of all polygons in a choropleth soil map.
#' @param polygon The id of the polygon from which to draw the soil classes.
#' @param n The number of samples to draw.
#' @param method The method of allocation. Permissible values are \code{"weighted"} (the default), \code{"random-mapunit"} and
#' \code{"random-all"}. See \emph{Details} for information on how the modes operate.
#'
.allocate <- function(composition, polygon, n = 15, method = "weighted")
{
  # Subset composition data frame
  polygon.comp <- subset(composition, composition$poly == polygon)
  
  if(nrow(polygon.comp) > 0)
  {
    if(method == "weighted")
    {
      # Draw from Dirichlet distribution
      s <- gtools::rdirichlet(1, polygon.comp$proportion)
      
      # Make prediction
      classes = base::sample(polygon.comp$soil_class, size = n, replace = TRUE, prob = s[1, ])
      
      return(classes)
    }
    else if (method == "random-mapunit")
    {
      # Make prediction
      classes = sample(polygon.comp$soil_class, size = n, replace = TRUE)
      
      return(classes)
    }
    else if (method == "random-all")
    {
      classes = sample(unique(composition$soil_class), size = n, replace = TRUE)
      
      return(classes)
    }
    else stop("mode \"", weighted, "\" is unknown")
  }
  else
  {
    stop(paste0("No polygon composition available for polygon with id ", polygon, "."))
  }
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
  r <- numeric(length = nrow(observations))
  meta <- list(r = numeric(length = nrow(observations)),
               type = base::rep("actual", nrow(observations)),
               method = base::rep("observed", nrow(observations)))
  obs <- cbind(as.data.frame(meta), as.data.frame(observations), o.covariates)
  
  # Return only those observations with no NA in the covariates
  obs <- obs[which(complete.cases(obs)), ]
  return(obs)
}