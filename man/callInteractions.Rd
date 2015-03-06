\name{callInteractions}
\alias{callInteractions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A statistical test to call pairwise interactions from interaction scores.}
\description{Using the four replicates per gene pair (from the two-by-two dsRNA designs) the null hypothesis the the interaction score is zero is tested by a moderated t-test (R-package limma). The p-values are adjusted by the method of Benjamini-Hochberg. The adjusted p-values and the pairwise interaction scores are stored in an HDF5 file.}
\usage{
callInteractions(filePI, fileInteractions, verbose = TRUE, overwrite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filePI}{Filename of the HDF5 file with interaction scores. (input)}
  \item{fileInteractions}{Filename of the HDF5 file with adjusted p-values and interaction scores per gene pair and feature (output).}
  \item{verbose}{Prints more output on screen.}
  \item{overwrite}{If TRUE, overwrite an existing output file (fileInteractions), otherwise stops.}
}
\value{NULL is returned. As a side effect the HDF5 file 'fileInteractions' is created and adjusted p-values and pairwise interaction scores are saved to this file.}
\references{
Horn T, Sandmann T, Fischer B, Axelsson E, Huber W, Boutros M (2011). \emph{Mapping of signaling networks through synthetic genetic interaction analysis by RNAi}. Nature Methods 8: 341-346.
}
\author{
Bernd Fischer
}
\seealso{\link{DmelSGI-package}}
\examples{
  print(callInteractions)
}
\keyword{ datasets }
