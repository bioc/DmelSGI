\name{orderSpiderAxis}
\alias{orderSpiderAxis}
\title{Orders the axis of a spider plot.}
\description{Solves a traveling saleperson problem to optimally order the arms of a spider plot.}
\usage{
orderSpiderAxis(X)
}
\arguments{
  \item{X}{A d x n matrix of values for n instances that are plotted in the spider plots with d arms.}
}
\details{The arms are ordered such that two neighboring spider arms are similar to each other.}
\value{An integer vector with an optimal order of the spider arms.}
\author{Bernd Fischer}
\seealso{\code{\link{grid.spider}},\link{DmelSGI-package}}
\examples{
print(orderSpiderAxis)
}
\keyword{ datasets }
