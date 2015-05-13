\name{Goodness of fit measures}
\docType{methods}
\alias{goof}
\title{Goodness of fit measures}
\description{This function performs a number of model goodness of fit measures which include the root mean square error (RMSE), mean square error (MSE),  prediction bias, coefficient of determination (R squared), concordance correlation coefficient (CCC), ratio of performance to deviation (RPD), ratio of performance to interquartile distance (RPIQ), and residual variance estimates. These goodness of fit measures are used in both digital soil mapping and soil chemometric studies to test the relative performance of competing models.}

\usage{goof(observed,predicted, plot.it = FALSE, type="DSM")}
\arguments{
  \item{observed}{numeric; a vector or matrix of values that are actual observations of some phenomenon.} 
  \item{predicted}{numeric; a vector or matrix of values of predictions of the phenomenon that was observed .}     
  \item{plot.it}{logical; If TRUE an xy-plot of the observations and predictions will be generated. }
  \item{type}{character; Selection from either ``DSM'' or ``spec'' to generate the goodness of fit statistics of greatest relevance to either digital soil mapping and soil spectroscopy respectively.}
}

\value{Returns a dataframe with the goodness of fit statistics. The column headers of the \code{dataframe} are: \code{R2} (coefficient of determination), \code{concordance} (concordance correlation coefficient), \code{MSE} (mean square error), \code{RMSE} (root mean square error), \code{bias} (prediction bias), \code{MSEc} (mean residual variance), \code{RMSEc} (root mean residual variance), \code{RPD} (ratio of performance to deviation), \code{RPIQ} (ratio of performance to interquartile distance). An xy-plot will also be generated if requested. }

\note{These goodness of fit measures are not exclusive to digital soil mapping or soil spectroscopy.}

\author{Brendan Malone}
\references{
\itemize{
\item Bellon-Maurel, V., Fernandez-Ahumada, E., Palagos, B., Roger, J., McBratney, A., (2010) \href{http://dx.doi:10.1016/j.trac.2010.05.006}{Critical review of chemometric indicators commonly used for assessing the quality of the prediction of soil attributes by NIR spectroscopy}. Trends in Analytical Chemistry, 29(9): 1073-1081.
\item Hastie, T., Tibshirani, R., Friedman, J., (2009)  The Elements of Statisitcal Learning. Springer Series in Statistics.}
}


\examples{
library(ithir)
library(MASS)


#some data
data(USYD_soil1)

#fit a linear model
mod.1 <- lm(CEC ~ clay, data = USYD_soil1 , y = TRUE, x = TRUE)

#Goodness of fit
goof(observed = mod.1$y, predicted = mod.1$fitted.values, plot.it = TRUE)

}
\keyword{methods}
