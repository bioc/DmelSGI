\name{qualityControlGene}
\alias{qualityControlGene}
\docType{data}
\title{Correlation of interaction profiles between independent dsRNA designs}
\description{The quality control of dsRNA designs of the template genes is described by the correlation of the interaction profile between two independent dsRNA designs. In addition the annotation of the genes  and a logical vector indicating which gene passed the quality control is given.}
\usage{data(qualityControlGene)}
\value{A dataset with the output of the gene quality control.}
\format{
  The format is:\cr
List of 3\cr
 $ correlation: num [1:1463] 0.876 0.877 0.897 0.889 0.865 ...\cr
 $ Annotation :'data.frame':	1463 obs. of  6 variables:\cr
  ..$ TID             : chr [1:1463(1d)] "FBgn0002441" "FBgn0038760" ...\cr
  ..$ TemplatePlate   : int [1:1463(1d)] 1 1 1 1 1 1 1 1 1 1 ...\cr
  ..$ group           : chr [1:1463(1d)] "sample" "sample" "sample" ...\cr
  ..$ Symbol          : chr [1:1463(1d)] "l(3)mbt" "MED25" "CG31156" ...\cr
  ..$ Name            : chr [1:1463(1d)] "lethal (3) malignant brain tumor"...\cr
 $ passed     : logi [1:1463] TRUE TRUE TRUE TRUE TRUE TRUE ...
}
\seealso{\link{DmelSGI-package}}
\examples{
data(qualityControlGene, package="DmelSGI")
str(qualityControlGene)
}
\keyword{datasets}
