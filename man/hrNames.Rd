\name{hrNames}
\alias{hrNames}

\title{Human readable feature names}
\description{Translate feature names to human readable feature
  names. Names not known to this function for conversion are returned unchanged.}
\usage{
hrNames(names)
}

\arguments{
  \item{names}{Original feature names.}
}
\value{
  A vector of translated feature names.
}
\author{
Bernd Fischer
}

\seealso{\link{DmelSGI-package}}
\examples{
hrNames(c("4x.count", "4x.ratioMitotic"))
}
\keyword{datasets}
