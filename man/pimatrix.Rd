\name{pimatrix}
\alias{pimatrix}
\docType{data}
\title{Pairwise genetic interaction scores per experiment (no summary per gene pair)}
\description{An array (D) with pairwise interaction scores. The interaction scores are given per experiment are not yet summarized per gene pair. See \link{Interactions} for interactions summarized per gene pair. A contains the annotation of genes and features.}
\usage{data(pimatrix)}
\value{The matrix of pairwise interaction scores per RNAi design.}
\format{
  The format is:\cr
List of 2\cr
 $ D   : num [1:1293, 1:2, 1:72, 1:2, 1:21] -0.1086 -0.1799 -0.045 -0.0664 ...\cr
  ..- attr(*, "dimnames")=List of 5\cr
  .. ..$ target      : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" "CG6833" ...\cr
  .. ..$ targetDesign: chr [1:2] "1" "2"\cr
  .. ..$ query       : chr [1:72(1d)] "CG31156" "lilli" "Smg1" "Axn" ...\cr
  .. ..$ queryDesign : chr [1:2] "1" "2"\cr
  .. ..$ phenotype   : chr [1:21] "4x.count" "4x.ratioMitotic" ...\cr
 $ Anno:List of 5\cr
  ..$ target      :'data.frame':	1293 obs. of  6 variables:\cr
  .. ..$ TID             : chr [1:1293(1d)] "FBgn0002441" "FBgn0038760" ...\cr
  .. ..$ TargetPlate     : int [1:1293(1d)] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ group           : chr [1:1293(1d)] "sample" "sample" "sample" ...\cr
  .. ..$ Symbol          : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" ...\cr
  .. ..$ Name            : chr [1:1293(1d)] "lethal (3) malignant brain tumor"...\cr
  ..$ targetDesign:'data.frame':	2 obs. of  1 variable:\cr
  .. ..$ design: int [1:2(1d)] 1 2\cr
  ..$ query       :'data.frame':	72 obs. of  5 variables:\cr
  .. ..$ TID             : chr [1:72(1d)] "FBgn0051156" "FBgn0041111" ...\cr
  .. ..$ Batch           : int [1:72(1d)] 1 1 1 1 1 1 2 2 2 2 ...\cr
  .. ..$ Symbol          : chr [1:72(1d)] "CG31156" "lilli" "Smg1" "Axn" ...\cr
  .. ..$ Name            : chr [1:72(1d)] "-" "lilliputian" "Smg1" "Axin" ...\cr
  ..$ queryDesign :'data.frame':	2 obs. of  1 variable:\cr
  .. ..$ design: int [1:2(1d)] 1 2\cr
  ..$ phenotype   :'data.frame':	21 obs. of  1 variable:\cr
  .. ..$ phenotype: chr [1:21] "4x.count" "4x.ratioMitotic" ...
}
\seealso{\link{Interactions},\link{DmelSGI-package}}
\examples{
data(pimatrix)
str(pimatrix)
}
\keyword{datasets}
