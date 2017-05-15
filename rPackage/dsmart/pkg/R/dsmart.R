


#' Disaggregating and harmonising soil map units through resampled classification trees
#' 
#' This function, together with companion function \code{dsmartR} implements
#' the DSMART (Disaggregating and harmonising soil map units through resampled
#' classification trees) algorithm as described in Odgers et al. (2014). This
#' is the workhorse function that involves multiple resampling, C5 decision
#' tree model fitting, and subsequent mapping in order to realize potential
#' candidate soil classes within aggregated soil mapping units. There is also 
#' added facility to incorporate point observed data into the algorithm too.
#'
#' @param covariates A \code{RasterStack} of \emph{scorpan} environmental
#' covariates to both guide the C5 model fitting and provide the
#' necessary environmental data for the spatial mapping of predicted soil
#' classes.
#' @param polygons See \code{data(dsT_polygons)} for an example of 
#' requirements. Each polygon needs to have the following necessary attribute
#' data: A numeric unique map unit identifier, a character map unit name or 
#' code. And coupled columns of map unit classes and their respective 
#' proportions.
#' @param composition See \code{data(dsT_composition)} for an example of
#' requirements. First column indicates a numeric unique map unit identifier.
#' Second column indicates character map unit name or code. Third column 
#' indicates map unit compositions. And fourth column indicates the percentage
#' proportion of map unit compositions. All map unit codes and their subsequent
#' compositions are row stacked together.
#' @param n numeric; number of samples to take from each soil mapping polygon 
#' for C5 model fitting.
#' @param reals numeric; number of C5 modeling fitting and mapping realisations
#' to implement.
#' @param cpus numeric; number of compute nodes to use. Default is 1.
#' @param obsdat {\code{data.frame}; Should be just 3 columns. Columns 1 and 2 
#' are the spatial coordinates, while column 3 is the target variable class. It
#' is assumed that the coordinates are on the same projection as the map to be
#' disaggregated and associated covariates. It is also assumed that the 
#' observed target variable classes correspond to the same classification level
#' as the map unit compositions. For example MapUnit1 contains soilclassA, 
#' soilClassB, and soilClassC; the observed data should therefore correspond to
#' either soilclassA, soilClassB, and soilClassC etc.
#'
#' @return
#' @export
#'
#' @examples
dsmart <- function(covariates, polygons, composition, n = 15, reals = 100, cpus = 1, obsdat = NULL)
{
  # Check arguments
  if(n < 1)
  {
    stop("n must be greater than 0")
  }
  if(reals < 1)
  {
    stop("reals must be greater than 0")
  }
  if(cpus < 1)
  {
    stop("cpus must be greater than 0")
  }
  
  # Generate lookup table
  names(composition) <- c("poly", "mapunit", "soil_class", "proportion")
  lookup = as.data.frame(sort(unique(composition$soil_class)))
  lookup$code = seq(from=1, to=nrow(lookup), by=1)
  colnames(lookup) = c("name", "code")
  
  if(!is.null(obsdat)){names(obsdat)<- c("x", "y", "class")}
  
  # Create subdirectories to store results in
  #model_lists <- vector("list", reals) #empty list 
  dir.create("output/", showWarnings = FALSE)
  dir.create("output/realisations", showWarnings = FALSE)
  dir.create("output/trees", showWarnings = FALSE)
  #dir.create("output/summaries", showWarnings = FALSE)
  
  #strg<- paste(getwd(), "/output/rasters/", sep = "")
  #strm<- paste(getwd(), "/output/models/", sep = "")
  #strs<- paste(getwd(), "/output/summaries/", sep = "")
  write.table(lookup, paste("classLookupTable.txt", sep = ""), sep = ",", col.names = TRUE, row.names = FALSE) 
  
  #pb <- txtProgressBar(min = 0, max = reals, style = 3)
  
  # Get samples for all realisations
  message(paste0("Generating samples for ", reals, " realisations"))
  samples <- getSamples(covariates, shp, composition, n.realisations = reals,
                        n.samples = n, method = "weighted")
  
  # Write samples to text file
  
  for (j in 1:reals)
  {
    message(paste0("Realisation ", j))
    
    # Extract samples for the current realisation
    s <- samples[which(samples$r == j), ]
    
    # Fit classification tree
    tree = C50::C5.0(s[, 6:ncol(s)], y = s$soil_class)
    #model_lists[[j]] <- tree
    
    # Save tree to text file
    out <- utils::capture.output(summary(tree))
    cat(out, file = paste0("output/trees/tree_", formatC(j, width = nchar(reals), format = "d", flag = "0"), ".txt"),
        sep = "\n", append = TRUE)
    
    # Save tree to rdata file
    save(tree, file = paste0("output/trees/tree_", formatC(j, width = nchar(reals), format = "d", flag = "0"), ".RData"))
    
    # Predict realisation and save it to raster
    raster::beginCluster(cpus)
    r1 <- raster::clusterR(covariates, predict, args = list(tree),
                           filename = paste0("output/realisations/r_", formatC(j, width = nchar(reals), format = "d", flag = "0"), ".tif"),
                           format = "GTiff", overwrite = TRUE, datatype = "INT2S")
    raster::endCluster()
    #plot(r1)
    #setTxtProgressBar(pb, j)
  }
  
  #Save models to file
  #save(model_lists, file = paste(paste(getwd(),"/dsmartOuts/",sep=""),"dsmartModels.RData", sep="") )
  
  #close(pb)
  #message(paste(paste("DSMART outputs can be located at:",getwd(), sep=" "), "/output/",sep="") )
}

#END