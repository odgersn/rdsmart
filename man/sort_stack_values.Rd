% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.r
\name{sort_stack_values}
\alias{sort_stack_values}
\title{Sort raster stack values}
\usage{
sort_stack_values(r, n = nlyr(r), decreasing = TRUE)
}
\arguments{
\item{r}{A SpatRaster whose layers contain the values to be sorted.}

\item{n}{An integer that identifies the number of layers in \code{r} that
should be sorted. Default is to sort all layers in \code{r}.}

\item{decreasing}{A boolean that indicates whether values should be sorted in
decreasing order (\code{TRUE}) or not (\code{FALSE}).}
}
\value{
A SpatRaster containing the sorted data.
}
\description{
#' This function sorts the values in a raster stack. Sorting is performed on
the vector of values at each grid cell. See \code{\link[base:sort]{base::sort()}} for more
information about how ordering works. A custom sorting algorithm written in
C++ is used to decrease processing time
}
