\name{Interactions}
\alias{Interactions}
\docType{data}
\title{Pairwise genetic interaction scores and p-values (summary per gene pairs)}
\description{Two arrays are provided in this dataset: The pairwise interaction scores summarized per gene pair and the respective adjusted p-values. p-values are computed by a moderated t-test (limma) and corrected for multiple testing by the method of Benjamini-Hochberg. The list (Anno) contains the annotation of the target genes, query genes, and features. See \link{pimatrix} for interaction scores that are not yet summarized per gene pair.
}
\usage{data(Interactions)}
\value{An object containing the pi-scores, the adjusted p-values and the annotation of the statistical genetic interactions.}
\format{
  The format is:\cr
List of 3\cr
 $ piscore: num [1:1293, 1:72, 1:21] -0.0866 -0.0924 -0.0707 -0.0878 -0.0587 ...\cr
  ..- attr(*, "dimnames")=List of 3\cr
  .. ..$ target   : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" "CG6833" ...\cr
  .. ..$ query    : chr [1:72(1d)] "CG31156" "lilli" "Smg1" "Axn" ...\cr
  .. ..$ phenotype: chr [1:21] "4x.count" "4x.ratioMitotic" ...\cr
 $ padj   : num [1:1293, 1:72, 1:21] 0.232 0.249 0.286 0.281 0.369 ...\cr
  ..- attr(*, "dimnames")=List of 3\cr
  .. ..$ target   : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" "CG6833" ...\cr
  .. ..$ query    : chr [1:72(1d)] "CG31156" "lilli" "Smg1" "Axn" ...\cr
  .. ..$ phenotype: chr [1:21] "4x.count" "4x.ratioMitotic" ...\cr
 $ Anno   :List of 3\cr
  ..$ target   :'data.frame':	1293 obs. of  6 variables:\cr
  .. ..$ TID             : chr [1:1293(1d)] "FBgn0002441" "FBgn0038760" ...\cr
  .. ..$ TargetPlate     : int [1:1293(1d)] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ group           : chr [1:1293(1d)] "sample" "sample" "sample" ...\cr
  .. ..$ Symbol          : chr [1:1293(1d)] "l(3)mbt" "MED25" "CG31156" ...\cr
  .. ..$ Name            : chr [1:1293(1d)] "lethal (3) malignant brain tumor" ...\cr
  ..$ query    :'data.frame':	72 obs. of  5 variables:\cr
  .. ..$ TID             : chr [1:72(1d)] "FBgn0051156" "FBgn0041111" ...\cr
  .. ..$ Batch           : int [1:72(1d)] 1 1 1 1 1 1 2 2 2 2 ...\cr
  .. ..$ Symbol          : chr [1:72(1d)] "CG31156" "lilli" "Smg1" "Axn" ...\cr
  .. ..$ Name            : chr [1:72(1d)] "-" "lilliputian" "Smg1" "Axin" ...\cr
  ..$ phenotype:'data.frame':	21 obs. of  1 variable:\cr
  .. ..$ phenotype: chr [1:21] "4x.count" "4x.ratioMitotic" ...
}
\seealso{
\link{pimatrix},\link{DmelSGI-package}
}
\examples{
data(Interactions)
str(Interactions)
}
\keyword{datasets}
