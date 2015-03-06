\name{toRaster}
\alias{toRaster}
\title{Converts a real valued matrix in a matrix of color codes printable by grid.raster}
\description{A matrix of real values in a matrix are converted in a matrix of RBG values that can be printed by \code{\link{grid.raster}}.}
\usage{
toRaster(x, cuts, col)
}
\arguments{
  \item{x}{A real valued matrix.}
  \item{cuts}{Break points for the color values. Length of \code{cuts} has to be the length of \code{col} plus one.}
  \item{col}{A vector of color values, e.g. as produced by \code{\link{colorRampPalette}}.}
}
\value{Returns a matrix of RBG color values that can be printed by \code{\link{grid.raster}}.}
\author{Bernd Fischer}
\seealso{\code{\link{grid.raster}},\code{\link{colorRampPalette}},\link{DmelSGI-package}}
\examples{
print(toRaster)
}
\keyword{ datasets }
