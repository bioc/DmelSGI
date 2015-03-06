\name{getBaseDir}
\alias{getBaseDir}
\title{Returns the base directory for the vignette.}
\description{Returns the base directory for the vignette. When knitr is applied on the main vignette, this function ensures that the subvignettes get knowledge of the base directory.}
\usage{
getBaseDir(default = ".")
}
\arguments{
  \item{default}{The default base directory.}
}
\value{Returns a character with the directory name}
\author{
Bernd Fischer
}
\seealso{\link{DmelSGI-package}}
\examples{
getBaseDir()
}
\keyword{datasets}
