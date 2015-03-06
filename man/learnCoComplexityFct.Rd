\name{learnCoComplexFct}
\alias{learnCoComplexFct}
\alias{convertCorrelations}
\title{Learning the co-complex function}
\description{Learning the co-complex function from a correlation matrix.}
\usage{
learnCoComplexFct(C, ProteinComplexes)
convertCorrelations(C, coComplexFct)
}
\arguments{
  \item{C}{A genes-by-genes correlation matrix (or matrix containing any other pairwise scores).}
  \item{ProteinComplexes}{A list of protein complexes. Each element of the list is one protein complex and contains a \code{data.frame} with at least one column \code{gene_id}.}
  \item{coComplexFct}{\code{coComplexFct} is an object containing the co-complex function as returned by \code{learnCoComplexFct}.}
}
\details{An empirical density function of the values in \code{C} is computed once for the gene pairs that are co-member of at least one protein complex and once for all other gene pairs. The definition of protein complexes that is used to learn the co-complex function is taken from \code{ProteinComplexes}.}
\value{
\code{learnCoComplexFct} will provide an object that contains the empirical density function of correlation.

\code{convertCorrelations} will provide a matrix of co-complex scores.
}
\references{A similar approach is used in Ryan, C.J., et al. (2012). Hierarchical modularity and the evolution of genetic interactomes across species. Molecular cell 46, 691-704.}
\author{Bernd Fischer}
\seealso{\link{DmelSGI-package}}
\examples{
  print(learnCoComplexFct)
  print(convertCorrelations)
}
\keyword{ datasets }
