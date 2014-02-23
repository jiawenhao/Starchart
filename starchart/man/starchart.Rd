\name{starchart}
\alias{starchart}
\alias{starchart.default}
\alias{print.starchart}
\alias{summary.starchart}
\alias{print.summary.starchart}
\alias{predict.starchart}
\title{Recursive Partitioning Using Regression Trees}
\description{This package implements Starchart, a regression tree-based tool
  that recursively partitions a (GPU) program's software and hardware design
  space using design parameter values as predictors and measured performance, 
  power, or another metric as responses.}
\usage{
\method{starchart}{default}(x, y, splittable = 10, fitted = FALSE, ...)
\method{print}{starchart}(x, ...)
\method{summary}{starchart}(object, ...)
\method{predict}{starchart}(object, data = NULL, ...)
}
\arguments{
  \item{x}{ a data frame containing design parameter values of sampled designs. }
  \item{y}{ a numeric vector of responses such as performance or power.}
  \item{splittable}{ stop splitting leaf nodes with less than these many samples. }
  \item{fitted}{ calculate fitted values and residuals for x, may take some time. }
  \item{object}{ a Starchart tree object. }
  \item{data}{ unknown design samples whose responses need to be predicted. }
  \item{\dots}{ additional options. }
}
\value{
  An object of class \code{starchart}, a list including elements
  \item{nodes}{ the name of the parameter each node splits on }
  \item{thresholds}{ the threshold value each node splits on }
  \item{lefts}{ the row numbers of left child nodes }
  \item{rights}{ the row numbers of right child nodes }
}
\references{http://www.princeton.edu/~wjia/starchart/}
\author{Wenhao Jia}
\seealso{\code{\link{plot.starchart}}}
\examples{
data("kmeans")
tree = starchart(kmeans[1:200, 1:3], kmeans[1:200, 4])
summary(tree)
print(tree)
plot(tree)
}
\keyword{tree}
