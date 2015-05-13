\name{dsmart}
\docType{methods}
\alias{dsmart}
\title{Disaggregating and harmonising soil map units through resampled classification trees}
\description{This function, together with companion function \code{dsmartR} implements the the dsmart (Disaggregating and harmonising soil map units through resampled classification trees) algorithm as described in Odgers et al. (2014). This is the workhorse function that involves multiple resampling, C5 decision tree model fitting, and subsequent mapping in order to realize potential candidate soil classes within aggregated soil mapping units.}

\usage{dsmart(covariates = NULL, polygons = NULL, composition = NULL, n=NULL, reals = NULL, cpus=1)}
\arguments{
  \item{covariates}{\code{RasterStack}; A stack of grids or rasters that are surrogates for environmental variables to both guide the C5 model fitting and provide the necessary environmental data for the spatial mapping of predicted soil classes.} 
  \item{polygons}{\code{SpatialPolygonsDataFrame}; See \code{data(dsT_polygons)} for an example of requirements. Each polygon needs to have the following necessary attribute data: A numeric unique map unit identifier, a character map unit name or code. And couple columns of map unit classes and their respective proportions.}     
  \item{composition}{\code{data.frame}; See \code{data(dsT_composition)} for an example of requirements. First column indicates a numeric unique map unit identifier. Second column indicates character map unit name or code. Third column incites map unit compositions. And fourth column indicates the percentage proportion of map unit compositions. All map unit codes and their subsequent compositions are row stacked together.}
  \item{n}{numeric; number of samples to take from each soil mapping polygon for C5 model fitting.}
  \item{reals}{numeric; number of C5 modeling fitting and mapping realisations to implement. }
  \item{cpus}{numeric; number of compute nodes to use. Default is 1.}
}

\value{This function returns a number of outputs, all of which are saved to the working directory in a folder called \code{dsmartOuts}. 
\describe{
   \item{\code{dsmartModels.Rdata}}{This \code{R} object contains in a list the model object from each realisation of the algorithm.}
   \item{\code{models}}{folder; This folder contains the summary output of each C5 model fitted from each realisation of the algorithm. }
   \item{\code{rasters}}{folder; This folder contains the the rasters each soil classes created from each model realisation and are labelled sequentially. This folder also contains the file \code{classLookupTable.txt} which is a lookup table to link unique numeric labels to soil class character names.}  
}}

\note{The \code{dsmart} function is a companion function to \code{dsmartR}. Although this function can be run in parallel mode, increasing the number of realisations can drastically result in longer computation times.}

\author{Nathan Odgers and Brendan Malone}
\references{
\itemize{
\item Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., (2014) \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{Disaggregating and harmonising soil map units through resampled classification trees}. Geoderma, 214-215: 91-100.
\item Odgers, N.P., Sun, McBratney, A.B., Minasny, B., (2015) \href{http://dx.doi.org/10.1016/j.geoderma.2014.09.009}{Digital soil property mapping and uncertainty estimation using soil class probability rasters}. Geoderma, 237:238 190-198.
}


\examples{
library(dsmart)



#load the relevent data

#Covariates
data(dsT_covariates)

#Polygons
data(dsT_polygons)

#Map unit compositions
data(dsT_composition)

#Run dsmart (with 15 samples per polygon, 5 C5 model realisations, using 1 compute node)
dsmart(covariates = dsT_covariates, polygons = dsT_polygons, composition = dsT_composition, n=15, reals = 5, cpus=1)

}
\keyword{methods}
