\name{fitepistasis}
\alias{fitepistasis}
\docType{data}
\title{Fit of the pi-score vectors as a function of main effects}
\description{The output of the linear fit of the pi-score vectors as a function of main effects. The list contains four datasets:
\itemize{
\item A 6-dimensional array of the original data with pi-scores and main effects,
\item the coefficients of the linear fit,
\item the p-values from Anova and
\item the variance explained by the main effects.
}
}
\usage{data(fitepistasis)}
\value{The output of the epistasis estimation}
\format{
  The format is:\cr
List of 2\cr
 $ Coef: num [1:3, 1:1293, 1:2, 1:72, 1:2] -0.0464 0.3702 -0.0475 -0.0794 ...\cr
  ..- attr(*, "dimnames")=List of 5\cr
  .. ..$ : chr [1:3] "const" "xt" "xq"\cr
  .. ..$ : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" "CG6833" ...\cr
  .. ..$ : chr [1:2] "1" "2"\cr
  .. ..$ : chr [1:72(1d)] "CG31156" "lilli" "Smg1" "Axn" ...\cr
  .. ..$ : chr [1:2] "1" "2"\cr
 $ Sq  : num [1:3, 1:1293, 1:2, 1:72, 1:2] 0.18515 0.00207 0.45705 0.15047 ...\cr
  ..- attr(*, "dimnames")=List of 5\cr
  .. ..$ : chr [1:3] "xt" "xq" "res"\cr
  .. ..$ : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" "CG6833" ...\cr
  .. ..$ : chr [1:2] "1" "2"\cr
  .. ..$ : chr [1:72(1d)] "CG31156" "lilli" "Smg1" "Axn" ...\cr
  .. ..$ : chr [1:2] "1" "2"
}
\seealso{\link{DmelSGI-package}}
\examples{
data(fitepistasis)
}
\keyword{datasets}
