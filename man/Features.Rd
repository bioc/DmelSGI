\name{Features}
\alias{Features}
\docType{data}
\title{Description of the extracted features.}
\description{The data.frame describing all features extraced from the images. Beside the name of the feature given as the row.name of the data.frame, 7 columns describe each extracted feature:
\itemize{
\item \bold{mag} is the magnification of the image that is used (either 4x or 10x)
\item \bold{summary} Features are extracted for each single cell. Features are summarized by mean, standard deviation, quantiles, and histograms.
\item \bold{mask} The segmentation mask used to extract the features, e.g. the nuclei, the nuclei in the pH3 channel, or the cell body extracted from the a-tubulin channel
\item \bold{channel} The channel used to extract the feature
\item \bold{set} The feature category
\item \bold{type} The type of feature (area, intensity, ...)
\item \bold{param} Additional parameters of the feature extraction (quantiles, histogram bins, ...)
}
}
\usage{data(Features)}
\value{A data.frame listing all extracted features.}
\format{
  The head of the data.frame is:\cr
                mag  summary        mask            channel set type   param\cr
4x.count        "4x" "nrNuclei"     "nucleus"       "DAPI"  "M" ""     ""   \cr
4x.countpH3     "4x" "nrNuclei"     "mitoticNuclei" "pH3"   "M" ""     ""   \cr
4x.isMitotic    "4x" "nrNuclei"     "nucleus"       "DAPI"  "M" ""     ""   \cr
4x.ratioMitotic "4x" "mitoticRatio" "nucleus"       "DAPI"  "M" ""     ""   \cr
4x.areaNuc      "4x" "mean"         "nucleus"       "DAPI"  "M" "area" ""   \cr
4x.areaNucSD    "4x" "stddev"       "nucleus"       "DAPI"  "M" "area" ""   \cr
}
\seealso{\link{DmelSGI-package}}
\examples{
data(Features)
head(Features)
}
\keyword{datasets}
