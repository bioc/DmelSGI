\name{SKDdata}
\alias{SKDdata}
\docType{data}
\title{Single knock down screen}
\description{D is the single knock down screen data for the 12 negative control query genes. The annotation of each dimension of D is provided is the list A.}
\usage{data(SKDdata)}
\value{A dataset with the single knockdown data.}
\format{
  The format is:
List of 2
 $ D   : num [1:1293, 1:2, 1:12, 1:21] 14.7 15 14.8 14.8 14.8 ...\cr
  ..- attr(*, "dimnames")=List of 4\cr
  .. ..$ target      : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" "CG6833" ...\cr
  .. ..$ targetDesign: chr [1:2] "1" "2"\cr
  .. ..$ batch       : chr [1:12] "1" "2" "3" "4" ...\cr
  .. ..$ phenotype   : chr [1:21] "4x.count" "4x.ratioMitotic" ...\cr
 $ Anno:List of 4\cr
  ..$ target      :'data.frame':	1293 obs. of  6 variables:\cr
  .. ..$ TID             : chr [1:1293(1d)] "FBgn0002441" "FBgn0038760" ...\cr
  .. ..$ TargetPlate     : int [1:1293(1d)] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ group           : chr [1:1293(1d)] "sample" "sample" "sample" ...\cr
  .. ..$ Symbol          : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" ...\cr
  .. ..$ Name            : chr [1:1293(1d)] "lethal (3) malignant brain tumor"...\cr
  ..$ targetDesign:'data.frame':	2 obs. of  1 variable:\cr
  .. ..$ design: int [1:2(1d)] 1 2\cr
  ..$ batch       :'data.frame':	12 obs. of  1 variable:\cr
  .. ..$ batch: int [1:12] 1 2 3 4 5 6 7 8 9 10 ...\cr
  ..$ phenotype   :'data.frame':	21 obs. of  1 variable:\cr
  .. ..$ phenotype: chr [1:21] "4x.count" "4x.ratioMitotic" ...
}
\seealso{\link{DmelSGI-package}}
\examples{
data(SKDdata)
str(SKDdata)
}
\keyword{datasets}
