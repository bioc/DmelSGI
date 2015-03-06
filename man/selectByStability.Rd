\name{selectByStability}
\alias{selectByStability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select features by stability.
}
\description{
Features are selected in a greedy manner. A linear model is fitted to estimate the remove the contribution of each features that can already be explaind by previously selected features. The next features is selected such, that the residuals are maximal correlated between replicates.
}
\usage{
selectByStability(subsample,
                  preselect = c("4x.count",
                                "4x.ratioMitotic",
                                "10x.meanNonmitotic.cell.0.s.area"),
                  Rdim = 40,
                  verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subsample}{A list with subsampled data as produced by \code{\link{subSampleForStabilitySelectionFct}}. \code{data(subSampleForStabilitySelection, package='DmelSGI')} provides the dataset used in the paper.}
  \item{preselect}{The names of features that are preselected, e.g. the features '4x.count' and '4x.ratioMitotic' should be selected, because they have a special interpretation.}
  \item{Rdim}{The maximum number of selected features.}
  \item{verbose}{If TRUE, more output is provided.}
}
\value{
\item{selected }{The names of the selected features. \code{Rdim} features will be reported.}
\item{correlation }{The correlation of the residual features.}
\item{ratioPositive }{The fraction of positively correlated residual features in each step of the selection process. The features with ratioPositive > 0.5 should be selected.}
\item{correlationAll }{The correlation of the all residual features in each step of the selection process.}
}
\author{
Bernd Fischer
}
\seealso{\link{DmelSGI-package}}
\examples{
  print(selectByStability)
}
\keyword{ datasets }
