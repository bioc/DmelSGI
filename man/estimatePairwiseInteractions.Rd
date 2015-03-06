\name{estimatePairwiseInteractions}
\alias{estimatePairwiseInteractions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimates pairwise interaction scores.}
\description{Estimates pairwise interaction scores for a large, multi-dimensional combinatorial screen.}
\usage{
estimatePairwiseInteractions(fileMatrixData,
                             filePI,
                             verbose = TRUE,
                             overwrite = FALSE,
                             useSKD = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileMatrixData}{Filename of HDF5 file containing the combinatoril screening data (input).}
  \item{filePI}{Filename of HDF5 file to save the pairwise interaction matrix (output).}
  \item{verbose}{Prints more output on screen.}
  \item{overwrite}{If TRUE, overwrite an existing output file (filePI), otherwise stops.}
  \item{useSKD}{If TRUE, the negative controls are used to estimate the overall effect.}
}
\details{Estimates the pairwise interaction scores for each feature and each batch by calling the function \code{\link{mymedpolish}}.}
\value{Returns TRUE. As a sideeffect, the array with interaction scores is stored in the HDF5 file 'filePI'.}
\references{
Horn T, Sandmann T, Fischer B, Axelsson E, Huber W, Boutros M (2011). \emph{Mapping of signaling networks through synthetic genetic interaction analysis by RNAi}. Nature Methods 8: 341-346.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{mymedpolish}},\link{DmelSGI-package}
}
\examples{
  print(estimatePairwiseInteractions)
}
\keyword{ datasets }
