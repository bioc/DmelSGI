\name{Intogen}
\alias{Intogen}
\docType{data}
\title{Intogen: Interactive Onco Genomics}
\description{A list of recurrently mutated genes from Intogen.}
\usage{data("Intogen")}
\value{The Intogen dataset.}
\format{
  A data frame with 2933 observations on the following 10 variables.
  \describe{
    \item{\code{gene}}{a character vector}
    \item{\code{symbol}}{a character vector}
    \item{\code{project.name}}{a character vector}
    \item{\code{mut.freq}}{a numeric vector}
    \item{\code{MuSiC}}{a character vector}
    \item{\code{oncodriveFM}}{a character vector}
    \item{\code{oncodriveCLUST}}{a character vector}
    \item{\code{ActiveDriver}}{a character vector}
    \item{\code{MutSig}}{a character vector}
    \item{\code{driver.category}}{a character vector}
  }
}
\seealso{\link{DmelSGI-package}}
\references{
Tamborero, David, et al. "Comprehensive identification of mutational cancer driver genes across 12 tumor types." Scientific reports 3 (2013).
}
\examples{
data(Intogen)
}
\keyword{datasets}
