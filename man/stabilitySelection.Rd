\name{stabilitySelection}
\alias{stabilitySelection}
\docType{data}
\title{The features selected by stability.}
\description{
The features selected by stability used in the selection process are available in this dataset. Furthermore, it contains the correlation of the residual features and the fraction of positive correlated features that is used as a stop criterion.
}
\usage{data(stabilitySelection)}
\value{An object containing the output of the feature selection.}
\format{
  The format is:\cr
List of 4\cr
 $ selected      : chr [1:50] "4x.count" "4x.ratioMitotic" ...\cr
 $ correlation   : num [1:50] 0.912 0.868 0.559 ...\cr
 $ ratioPositive : num [1:50] 1 0.957 0.95 ...\cr
 $ correlationAll:List of 50\cr
  ..$ : Named num [1:162] 0.912 0.946 0.948 ...\cr
  .. ..- attr(*, "names")= chr [1:162] "4x.count" "4x.countpH3" ...\cr
  ..$ : Named num [1:161] 0.946 0.795 0.868 ...\cr
  .. ..- attr(*, "names")= chr [1:161] "4x.countpH3" "4x.isMitotic" ...\cr
  ..$ : Named num [1:160] 0.946 0.129 0.927 ...\cr
  .. ..- attr(*, "names")= chr [1:160] "4x.countpH3" "4x.isMitotic" ...\cr
  ...\cr
}
\seealso{\link{DmelSGI-package}}
\examples{
data(stabilitySelection)
str(stabilitySelection)
}
\keyword{datasets}
