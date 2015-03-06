\name{mymedpolish}
\alias{mymedpolish}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A variant of the medpolish function.}
\description{A variant of the R-implementation of \code{\link{medpolish}}. Fits an additive model using Tukey's median polish procedure. The variant uses negative controls to estimate the overall effect and it estimates column effects separately for each batch of rows indicated by TP (for template plate). \bold{It is highly recommended to use the original function from the R-package stats.} The function is used by \code{\link{estimatePairwiseInteractions}}.
}
\usage{
mymedpolish(x, TP, TemplateNeg, QueryNeg, eps = 1e-04, maxiter = 100, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric matrix.}
  \item{TP}{integer. Specifying the plate number for each row element (template genes). Counting starts with 1. Column effects (query genes) will be estimated for each template plate separately}
  \item{TemplateNeg}{Index of negative controls on template plate}
  \item{QueryNeg}{Index of negative controls in the columns (queries)}
  \item{eps}{real number greater than 0. A tolerance for convergence: see \code{\link{medpolish}}.}
  \item{maxiter}{the maximum number of iterations}
  \item{na.rm}{logical. Should missing values be removed?}
}
\value{
  \item{neg }{the fitted constant term.or overall effect representing the effect of negative controls}
  \item{templateMainEffect }{the fitted row effects representing the single knock down effects of the template genes}
  \item{queryMainEffect }{the fitted column effects representing the single knock down effects of the query genes. It is a matrix with dimensions query genes x template plates}
  \item{pi }{the residuals which are the pairwise interaction scores}
}
\references{
Tukey, J. W. (1977). \emph{Exploratory Data Analysis}, Reading Massachusetts: Addison-Wesley.\cr
Horn T, Sandmann T, Fischer B, Axelsson E, Huber W, Boutros M (2011). \emph{Mapping of signaling networks through synthetic genetic interaction analysis by RNAi}. Nature Methods 8: 341-346.
}
\author{
Original implementation in the R package stats (See \code{\link{medpolish}}). Changes for estimating pairwise interaction scores by Bernd Fischer.
}
\seealso{
\code{\link{medpolish}},\code{\link{estimatePairwiseInteractions}},\link{DmelSGI-package}
}
\examples{
  print(mymedpolish)
}
\keyword{ datasets }
