\name{dsmartR}
\docType{methods}
\alias{dsmartR}
\title{Disaggregating and harmonising soil map units through resampled classification trees}
\description{This function, together with companion function \code{dsmart} implements the the dsmart (Disaggregating and harmonising soil map units through resampled classification trees) algorithm as described in Odgers et al. (2014). Essentially this function interogates C5 model predicted soil class outputs (rasters) from \code{dsmart} to create maps of the n-most probable soil classes, together with probability raster of each candidate soil class.}

\usage{dsmartR(rLocs= NULL, nprob = NULL, sepP=FALSE, lookup= NULL, cpus=1,param=NULL,param2=NULL)}
\arguments{
  \item{rLocs}{\code{RasterStack}; a stack of the individual soil class rasters that are created from the companion \code{dsmart} function.} 
  \item{nprob}{numeric; the number of n-most probable soil class maps to create. 1 indicates the most probable soil class most only. 2 would indicate the most probable and second most probable soil class maps are created. A \code{RasterBrick} of all probable rasters (the number equals the number of candidate soil classes) is automatically created, yet requesting the n-most probable will create individually the requested number of n-probable maps.}     
  \item{sepP}{logical; If \code{TRUE} probability rasters of each candidate soil class are created.}
  \item{lookup}{\code{data.frame}; An output from the companion \code{dsmart} function. This data links a unique numerical value to each candidate soil class code such as a character name.}
  \item{cpus}{numerical; numeric; number of compute nodes to use. Default is 1.}
}

\value{This function returns 3 folders to the working directory: nProbable, probabilities, and counts. nProbable contains the n-most probable soil class maps. probabilities contains the probability rasters for each candidte soil class if requested. counts contains a \code{RasterBrick} of the count of predicted candidate soil classes from the output realisations of \code{dsmart}. A \code{list} of the n-most probable rasters (\code{RasterStack}) and probabilitiy rasters (if requested) are saved to the working memory.}

\note{The \code{dsmart} function is a companion function to \code{dsmartR}. Together both functions are computationally demanding. The reprieve for this is that very little memory is stored in the global environment, with rasters being written directly to file. Similarly for greater efficency is is also better to keep input rasters (\code{rLocs}) in file, rather than load to memory.}

\author{Nathan Odgers and Brendan Malone}

\references{
\itemize{
\item Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., (2014) \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{Disaggregating and harmonising soil map units through resampled classification trees}. Geoderma, 214-215: 91-100.
\item Odgers, N.P., Sun, McBratney, A.B., Minasny, B., (2015) \href{http://dx.doi.org/10.1016/j.geoderma.2014.09.009}{Digital soil property mapping and uncertainty estimation using soil class probability rasters}. Geoderma, 237:238 190-198.
}


\examples{
library(dsmart)
library(raster)


#set parameters
data(dsT_lookup) #lookup table
data(dsmartOutMaps) #Results from 5 realiations of dsmart


#run function getting 3 most probable and creating probability raster using 1 compute core.
dsmartR(rLocs= dsmartOutMaps, nprob = 3, sepP=TRUE, lookup= dsT_lookup, cpus=1)

}
\keyword{methods}
