\name{DPiM}
\alias{DPiM}
\docType{data}
\title{Drosophila Protein Interaction Map (DPiM)}
\description{
The Drosophila Protein Interaction Map (DPiM) dataset contains of two parts: (1.) The experimentally identified protein interaction partners and (2.) inferred protein complexes.
}
\usage{data("DPiM")}
\value{The mass spec data and the inferred protein complexes of the DPiM dataset.}
\format{
  The format is:
List of 2
 $ interactions:'data.frame':	10969 obs. of  5 variables:
  ..$ Interactor_1     : chr [1:10969] "FBgn0036918" "FBgn0031143" "FBgn0030086" "FBgn0015019" ...
  ..$ Interactor_2     : chr [1:10969] "FBgn0037893" "FBgn0035102" "FBgn0033342" "FBgn0037632" ...
  ..$ HGSCore          : num [1:10969] 742 737 733 730 726 ...
  ..$ Evidence.in.DroID: chr [1:10969] "human_orthology; yeast_orthology" "" "human_orthology" "human_orthology; yeast_orthology" ...
  ..$ evidence         : logi [1:10969] TRUE FALSE TRUE TRUE TRUE TRUE ...
 $ complexes   :List of 556
  ..$ : chr [1:40] "FBgn0002284" "FBgn0002787" "FBgn0004066" "FBgn0010590" ...
  ..$ : chr [1:43] "FBgn0000212" "FBgn0000499" "FBgn0001276" "FBgn0001324" ...
  ..$ : chr [1:28] "FBgn0003660" "FBgn0011288" "FBgn0011708" "FBgn0013343" ...
  ..$ : chr [1:11] "FBgn0022023" "FBgn0025582" "FBgn0029629" "FBgn0033902" ...
  .. [list output truncated]
}
\seealso{\link{DmelSGI-package}}
\references{
Guruharsha, K. G., et al. "A protein complex network of Drosophila melanogaster." Cell 147.3 (2011): 690-703.
}
\examples{
data(DPiM)
}
\keyword{datasets}
