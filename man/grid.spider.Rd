\name{grid.spider}
\alias{grid.spider}
\alias{grid.spider.legend}
\title{Spider plot}
\description{A grid function to draw a spider plot.}
\usage{
grid.spider(v, col, col.arms = "black", dlim = NULL)
grid.spider.legend(vn, col.arms = "black", dlim = NULL)
}
\arguments{
  \item{v}{A vector of numbers to be presented in the spider plot.}
  \item{vn}{A vector of dimension names that are represented by the spider arms. Has the same length as v.}
  \item{col}{The color of the polygon area.}
  \item{col.arms}{The color of the background of the spider arms.}
  \item{dlim}{A vector with two values. Limits of the spider arm axis.}
}
\details{These function draw a grid spider plot or a legend for the spider arms.}
\value{Both functions return an invisible NULL, but they have an site-effect that draws a spider plot using grid.}
\author{Bernd Fischer}
\seealso{ \code{\link{orderSpiderAxis}},\link{DmelSGI-package} }
\examples{
  print(grid.spider)
  print(grid.spider.legend)
}
\keyword{ datasets }
