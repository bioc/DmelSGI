\name{datamatrix}
\alias{datamatrix}
\docType{data}
\title{Pairwise perturbation screen data}
\description{An data array (D) with the feature data of the pairwise perturbation screen. A contains the annotation of genes and features}
\usage{data(datamatrix)}
\value{An array with the phenotypic data.}
\format{
  The format is:\cr
List of 2\cr
 $ D   : num [1:1293, 1:2, 1:72, 1:2, 1:21] 14.8 14.9 14.9 14.9 14.9 ...\cr
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
  .. ..$ Name            : chr [1:1293(1d)] "malignant brain tumor"  ...\cr
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
\seealso{\link{DmelSGI-package}}
\examples{
data(datamatrix)
str(datamatrix)
}
\keyword{datasets}
