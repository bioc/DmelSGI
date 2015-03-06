\name{myHeatmap}
\alias{myHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Draws a heatmap for a three dimensional array}
\description{Draws a heatmap for a three dimensional array as e.g. the three dimensional genetic interaction cube.}
\usage{
myHeatmap(x, cuts, col, fontsize = 18, colnames = TRUE, rownames = FALSE)
}
\arguments{
  \item{x}{A three dimensional array.}
  \item{cuts}{break points for mapping the values in \code{x} to \code{col}. The length of \code{cuts} is one larger than the length of \code{col}.}
  \item{col}{A color bar as returned by \code{\link{colorRampPalette}}.}
  \item{fontsize}{The size of the text labels.}
  \item{colnames}{Logical. If TRUE, the column names are printed.}
  \item{rownames}{Logical. If TRUE, the row names are printed.}
}
\details{This function is used to draw the heatmap of the three dimensional genetic interaction cube.}
\value{Nothing is returned, but the function plots a heatmap as a site-effect.}
\author{Bernd Fischer and Wolfgang Huber}
\seealso{\link{DmelSGI-package}}
\examples{
print(myHeatmap)
}
\keyword{ datasets }
