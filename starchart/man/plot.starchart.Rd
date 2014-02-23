\name{plot.starchart}
\alias{plot.starchart}
\title{Plot a Starchart Tree}
\description{This packages plots a Starchart tree using R's built-in hclust 
and dendrogram support.}
\usage{
\method{plot}{starchart}(x, edge.labels = FALSE, height.sse = FALSE, top.splits = 0, ...)
}
\arguments{
  \item{x}{ a Starchart tree. }
  \item{edge.labels}{ overlay split conditions on internal edges. }
  \item{height.sse}{ draw nodes at heights proportional to SSEs. }
  \item{top.splits}{ only draw the top N splits, 0 to draw the whole tree. }
  \item{\dots}{ additional options. }
}
\references{http://www.princeton.edu/~wjia/starchart/}
\author{Wenhao Jia}
\seealso{\code{\link{starchart}}}
\examples{
data("kmeans")
tree = starchart(kmeans[1:200, 1:3], kmeans[1:200, 4])
plot(tree, edge.labels = TRUE, top.splits = 5)}
\keyword{tree}
