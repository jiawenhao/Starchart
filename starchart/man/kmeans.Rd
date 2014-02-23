\name{kmeans}
\alias{kmeans}
\docType{data}
\title{Rodinia Kmeans swap Kernel Design Samples}
\description{Performance and power measurements of 500 randomized kernel
designs running Rodinia's Kmeans swap kernel on AMD Radeon HD 7970.}
\usage{data(kmeans)}
\format{
  A data frame with 500 observations on the following 5 variables.
  \describe{
    \item{\code{ppb}}{a numeric vector}
    \item{\code{tpp}}{a numeric vector}
    \item{\code{consec}}{a numeric vector}
    \item{\code{runtime}}{a numeric vector}
    \item{\code{wattage}}{a numeric vector}
  }
}
\source{http://www.princeton.edu/~wjia/starchart.pdf}
\references{http://www.princeton.edu/~wjia/starchart/}
\examples{
data(kmeans)
print(kmeans)
}
\keyword{datasets}
