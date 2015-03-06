\name{subSampleForStabilitySelection}
\alias{subSampleForStabilitySelection}
\docType{data}
\title{
A subsampled dataset for use stability selection.
}
\description{
This dataset contains a subsample of the interaction screen for use in the function \code{\link{selectByStability}}. It contains the data matrix $D$ with 3000 experiments x 2 replicates (dsRNA designs) x 162 features.
}
\usage{data(subSampleForStabilitySelection)}
\value{The dataset of samples used for feature selection.}
\format{
  The format is:\cr
List of 3\cr
 $ D        : num [1:3000, 1:2, 1:162] -1.447 0.351 0.44 ...\cr
 $ Sample   : int [1:3000] 77465 94252 95176 ...\cr
 $ phenotype: chr [1:162(1d)] "4x.count" "4x.countpH3" ...\cr
}
\seealso{
\link{selectByStability},\link{subSampleForStabilitySelectionFct},\link{DmelSGI-package}
}
\examples{
data(subSampleForStabilitySelection)
str(subSampleForStabilitySelection)
}
\keyword{datasets}
