\name{plot2Phenotypes}
\alias{plot2Phenotypes}
\alias{plotPIdata}
\title{Plot directed epistatic interactions.}
\description{Plots the data to estimate a directed epistatic interactions.}
\usage{
plot2Phenotypes(X, gt, gq, f1, f2, length = 1, ...)
plotPIdata(X, gt, gq, show = "summary", ...)
}
\arguments{
  \item{X}{A 6-dimensional array (phenotype x [xt, xq, pi] x target genes x targetDesigns x query genes x queryDesigns).}
  \item{gt}{The target gene name.}
  \item{gq}{The query gene name.}
  \item{f1}{The first phenotypic feature.}
  \item{f2}{The second phenotypic feature.}
  \item{length}{Length of arrow head.}
  \item{show}{Either show='summary' to show the mean over all dsRNA designs or it is a vector of length 2 that specifies the two dsRNAs to show.}
  \item{\dots}{Other arguments passed to plot.}
}
\details{\code{plot2Phenotypes} shows a plot showing the two phenotypes on the axis. 2 x 2 arrows for the single gene effects of the two dsRNA designs for the two genes are shown in green and purple, the expected double knock-down effects under the non-interacting model for the 4 combinations of dsRNA designs in gray, and the four measured double knock-down effects in black. The black arrows are the genetic interactions.

\code{plotPIdata} shows two scatter plots for the fit of the vector of pairwise interaction scores across all phenotypes as a function of the single gene effects. Each dot represents one phenotype.}
\value{Nothing is returned, but the function draws a plot as a site-effect.}
\author{Bernd Fischer}
\seealso{\link{DmelSGI-package}}
\examples{
print(plot2Phenotypes)
print(plotPIdata)
}
\keyword{ datasets }
