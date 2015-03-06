\name{qualityControlFeature}
\alias{qualityControlFeature}
\docType{data}
\title{Correlation of features between replicates}
\description{The quality control of features is described by three vectors: The correlation between two replicates, the fraction of finite values, and a logical vector indicating which feature passed the quality control. The features are described in the dataset \link{Features}.}
\usage{data(qualityControlFeature)}
\value{A data.frame containing the output of the feature quality control.}
\format{
List of 3
 $ correlation      : num [1:328] 0.933 0.927 0.927 0.922 0.97 ...
 $ ratioFiniteValues: num [1:328] 0.999 0.999 0.999 0.999 0.999 ...
 $ passed           : logi [1:328] TRUE TRUE TRUE TRUE TRUE TRUE ...
}
\seealso{\link{DmelSGI-package}}
\examples{
data(qualityControlFeature, package="DmelSGI")
str(qualityControlFeature)
}
\keyword{datasets}
