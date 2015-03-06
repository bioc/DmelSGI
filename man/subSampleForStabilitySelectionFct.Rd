\name{subSampleForStabilitySelectionFct}
\alias{subSampleForStabilitySelectionFct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subsampling the data for stability selection.
}
\description{
The data is subsampled to reduce computation time and memory demand for stability selection. 10000 experiments are selected by chance to estimate the most stable directions.
}
\usage{
subSampleForStabilitySelectionFct(fileMatrixData, N = 10000, random.seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileMatrixData}{Filename of the HDF5 file containing the data array}
  \item{N}{Number of experiments used for stabilit selection}
  \item{random.seed}{If not NULL, the random.seed is set before sampling, to generate a reproducible analysis script}
}
\details{For each dsRNA design, query gene, and feature, the median value is substracted and the data are devided by the median deviation.}
\value{
  \item{D }{An array of dimension N x 2 x F, where N is the given number of experiments used, 2 represents the two query dsRNA designs, and F is the number of features.}
  \item{Sample }{The index of the sampled elements.}
  \item{phenotype }{The name of the features.}
}
\author{
Bernd Fischer
}
\seealso{\link{selectByStability},\link{subSampleForStabilitySelection},\link{DmelSGI-package}}
\examples{
  print(subSampleForStabilitySelectionFct)
}
\keyword{ datasets }
