\name{TID2HUGO}
\alias{TID2HUGO}
\docType{data}
\title{Mapping of flygene names to human.}
\description{Mapping of flybase gene identifier to their human orthologues. It is a one-to-many mapping.}
\usage{data(TID2HUGO)}
\value{A list containing the conversion of gene identifier from fly to human.}
\format{
  The format is:\cr
List of 1293\cr
 $ FBgn0002441: chr [1:3] "L3MBTL1" "L3MBTL3" "L3MBTL4"\cr
 $ FBgn0038760: chr "MED25"\cr
 $ FBgn0051156: chr "SRBD1"\cr
 ... [list output truncated]\cr
 - attr(*, "dim")= int 1293\cr
 - attr(*, "dimnames")=List of 1\cr
  ..$ : chr [1:1293] "FBgn0002441" "FBgn0038760" "FBgn0051156" "FBgn0036405" ...
}
\seealso{\link{DmelSGI-package}}
\examples{
data(TID2HUGO)
}
\keyword{datasets}
