\name{mainEffects}
\alias{mainEffects}
\docType{data}
\title{Main effects (single knock down effects) estimated from the combinatorial data}
\usage{data(mainEffects)}
\description{The overall effect (effect of the negative control experiments), and the estimated main effects (single knock down effects) for template and query genes. Overall effects and template main effects are estimated separately for each batch (1:12). Query main effects are estimated separately for each template plate (1:4)}
\value{An object containing the query and target main effects.}
\format{
  The format is:\cr
List of 4\cr
 $ overall: num [1:2, 1:12, 1:21] 14.6 14.6 14 14 14.7 ...\cr
  ..- attr(*, "dimnames")=List of 3\cr
  .. ..$ targetDesign: chr [1:2] "1" "2"\cr
  .. ..$ batch       : chr [1:12] "1" "2" "3" "4" ...\cr
  .. ..$ phenotype   : chr [1:21] "4x.count" "4x.ratioMitotic" ...\cr
 $ query  : num [1:72, 1:2, 1:4, 1:21] 0.0409 0.0192 -0.587 -0.0636 -0.4709 ...\cr
  ..- attr(*, "dimnames")=List of 4\cr
  .. ..$ query        : chr [1:72(1d)] "CG31156" "lilli" "Smg1" "Axn" ...\cr
  .. ..$ queryDesign  : chr [1:2] "1" "2"\cr
  .. ..$ templatePlate: chr [1:4] "1" "2" "3" "4"\cr
  .. ..$ phenotype    : chr [1:21] "4x.count" "4x.ratioMitotic" ...\cr
 $ target : num [1:1293, 1:2, 1:12, 1:21] 0.0921 0.3444 0.1603 0.2252 0.1293 ...\cr
  ..- attr(*, "dimnames")=List of 4\cr
  .. ..$ target      : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" "CG6833" ...\cr
  .. ..$ targetDesign: chr [1:2] "1" "2"\cr
  .. ..$ batch       : chr [1:12] "1" "2" "3" "4" ...\cr
  .. ..$ phenotype   : chr [1:21] "4x.count" "4x.ratioMitotic" ...\cr
 $ Anno   :List of 7\cr
  ..$ target       :'data.frame':	1293 obs. of  6 variables:\cr
  .. ..$ TID             : chr [1:1293(1d)] "FBgn0002441" "FBgn0038760" ...\cr
  .. ..$ TargetPlate     : int [1:1293(1d)] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ group           : chr [1:1293(1d)] "sample" "sample" "sample" ...\cr
  .. ..$ Symbol          : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" ...\cr
  .. ..$ Name            : chr [1:1293(1d)] "lethal (3) malignant brain tumor"...\cr
  ..$ targetDesign :'data.frame':	2 obs. of  1 variable:\cr
  .. ..$ design: int [1:2(1d)] 1 2\cr
  ..$ query        :'data.frame':	72 obs. of  5 variables:\cr
  .. ..$ TID             : chr [1:72(1d)] "FBgn0051156" "FBgn0041111" ...\cr
  .. ..$ Batch           : int [1:72(1d)] 1 1 1 1 1 1 2 2 2 2 ...\cr
  .. ..$ Symbol          : chr [1:72(1d)] "CG31156" "lilli" "Smg1" "Axn" ...\cr
  .. ..$ Name            : chr [1:72(1d)] "-" "lilliputian" "Smg1" "Axin" ...\cr
  ..$ queryDesign  :'data.frame':	2 obs. of  1 variable:\cr
  .. ..$ design: int [1:2(1d)] 1 2\cr
  ..$ phenotype    :'data.frame':	21 obs. of  1 variable:\cr
  .. ..$ phenotype: chr [1:21] "4x.count" "4x.ratioMitotic" ...\cr
  ..$ batch        :'data.frame':	12 obs. of  1 variable:\cr
  .. ..$ batch: int [1:12] 1 2 3 4 5 6 7 8 9 10 ...\cr
  ..$ templatePlate:'data.frame':	4 obs. of  1 variable:\cr
  .. ..$ templatePlate: int [1:4] 1 2 3 4
}
\seealso{\link{DmelSGI-package}}
\examples{
data(mainEffects)
str(mainEffects)
}
\keyword{datasets}
