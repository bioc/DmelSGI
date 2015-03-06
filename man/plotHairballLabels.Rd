\name{plotHairballLabels}
\alias{plotHairballLabels}
\title{Adds the cluster labels to a graph}
\description{Adds the cluster labels to a graph}
\usage{
plotHairballLabels(g, co, Labels, Col)
}
\arguments{
  \item{g}{An igraph object for the graph.}
  \item{co}{A n x 2 matrix of layout coordinates as returned by the igraph layout algorithms.}
  \item{Labels}{A list with a vector of gene names per cluster as they appear as verteces in the igraph objects. The list element names are printed to the graph plot.}
  \item{Col}{The colors of the cluster names.}
}
\value{Nothing is returned, but as a site-effect, the labels for the hairball are added to a plot.}
\author{Bernd Fischer}
\seealso{\link{DmelSGI-package}}
\examples{
print(plotHairballLabels)
}
\keyword{ datasets }
