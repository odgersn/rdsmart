\name{Edgeroi Land Class Points}
\docType{data}
\alias{edgeLandClass}
\title{Point data of estimated land classes from the Edgeroi District, NSW, Australia.}
\description{A \code{dataframe} of 500 point locations in the Edgeroi District, NSW, Australia (30.11S 149.66E) where there is an estimated Land Classification. The Land classification is based upon an unsupervised classification of unknown date Landsat 7 ETM+ spectral bands. The points represent a random sample of grid locations of the land classification map. The are 6 land classes and have a numeric identifier. The land classes roughly correspond to: 1) dense forest; 2) open forest; 3) water bodies; 4) woody vegetation and native grassland; 5) Irrigated cropping; and 6) dry land cropping. The CRS of the points is WGS84 UTM Zone 55.}

\usage{data(edgeLandClass)}
\format{
\code{edgeLandClass} is a 500 row \code{dataframe} with locational (2 columns) and land class information (1 column)
}
\details{This data frame is typical point data of environmental phenomena. }

\references{
\itemize{
\item This data is copyright of the Soil Security Lab, The University of Sydney. If any part of this data is to be used in any publication or report, please provide a citation:
Soil Security Laboratory, 2015. Use R for Digital Soil Mapping Manual. The University of Sydney, Sydney, Australia.
}
}

\examples{
library(ithir)

data(edgeLandClass)

#data summary
summary(edgeLandClass$LandClass)


}
\keyword{datasets}
