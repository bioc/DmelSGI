\name{applyDimensionReduction}
\alias{applyDimensionReduction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subsets the features in a genetic interaction dataset in HDF5 format. }
\description{Subsets the features in a genetic interaction dataset in HDF5 format. The features are selected by \code{\link{selectByStability}} beforehand.
}
\usage{
  applyDimensionReduction(fileMatrixData, fileNew, selected, 
                          verbose = TRUE, overwrite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileMatrixData}{(Input) data matrix in HDF5 format.}
  \item{fileNew}{(Output) data matrix in HDF5 format with a subset of features.}
  \item{selected}{The names of the selected features that will be subsetted.}
  \item{verbose}{Prints more output on screen.}
  \item{overwrite}{If TRUE, overwrite an existing output file (fileNew), otherwise stops.}
}
\value{NULL is returned. As a side effect the HDF5 file 'fileNew' is created and data matrices with subsetted features will be written to it.}
\author{
Bernd Fischer
}
\seealso{
\code{\link{selectByStability}},\link{stabilitySelection},\link{DmelSGI-package}
}
\examples{
  print(applyDimensionReduction)
}
\keyword{ datasets }
