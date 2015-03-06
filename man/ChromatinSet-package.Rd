\name{DmelSGI-package}

\alias{DmelSGI-package}
\alias{DmelSGI}

\docType{package}

\title{DmelSGI.}

\description{
The package contains the data and the source code to reproduce the
results and figures from the paper \emph{title TBC}.
}

\section{Package content}{

  See vignette("DmelSGI") for more detail on how to obtain the data used for specific figures. In addition this vignette contains the complete analysis and the generation of all figures.

%  The main screen can be loaded by \code{data("Dmel2PPMAPK", package="RNAinteractMAPK")}. Access to the pairwise interaction data is done via the \code{\link{getData}} function from the \code{\link{RNAinteract}}-package. See example below.

The following \bold{datasets} are provided with this package:
\itemize{
  \item \bold{Features and quality control}
    \itemize{
      \item \code{\link{Features}} Description of the extracted features.
      \item \code{\link{qualityControlFeature}} Correlation of features between replicates
      \item \code{\link{qualityControlGene}} Correlation of interaction profiles between independent dsRNA designs
    }
  \item \bold{Stability selection}
    \itemize{
      \item \code{\link{subSampleForStabilitySelection}} A subsampled dataset used to select features.
      \item \code{\link{stabilitySelection}} The features selected by stability.
    }
  \item \bold{Pairwise interaction matrix}
    \itemize{
      \item \code{\link{datamatrix}} Pairwise perturbation screen data
      \item \code{\link{pimatrix}} Pairwise genetic interaction scores per experiment (no summary per gene pair)
      \item \code{\link{Interactions}} Pairwise genetic interaction scores and p-values (summary per gene pairs)
      \item \code{\link{mainEffects}} Main effects (single knock down effects) estimated from the combinatorial data
      \item \code{\link{SKDdata}} Single knock down screen
    }
}

Functions is this package:
\itemize{
  \item \bold{stability selection}
    \itemize{
      \item \code{\link{subSampleForStabilitySelectionFct}}
      \item \code{\link{stabilitySelection}}
      \item \code{\link{applyDimensionReduction}}
    }
  \item \bold{pairwise interactions}
    \itemize{
      \item \code{\link{estimatePairwiseInteractions}}
      \item \code{\link{mymedpolish}}
      \item \code{\link{callInteractions}}
    }
}

%  Within this package a number of specialized functions is defined written for the analysis of the MAPK interaction screen and additional experiments shown in the paper. These functions are not intended to be general purpose analysis functions, but should document the steps of analysis of the paper. Therefore, these functions are not exported. A list of functions is given below.
%  A general purpose package for the analysis of genetic interaction screens is the package \code{\link{RNAinteract}}.

%  The following functions are provided within this package. 

%  Functions used for the classification: \code{MAPK.predict.classification}, \code{MAPK.cv.classifier}, \code{MAPK.getCV}, \code{MAPK.ternary.plot}, \code{MAPK.getXY}, \code{MAPK.plot.classification}.

%  Functions for the analysis of the interaction surfaces: \code{MAPK.plot.TPS.single}, \code{MAPK.plot.TPS.all}, \code{MAPK.estimate.TPS}, \code{MAPK.cv.TPS}, \code{MAPK.screen.as.array}.
  
%  A function to write the hit list: \code{MAPK.report.gene.lists.paper}.

%  A function to plot a heatmap: \code{MAPK.plot.heatmap.raster}.

%  A function to plot smaooth scatters: \code{MAPK.smooth.scatter}.

}

\details{
See vignette("DmelSGI") for details.
}

\author{
Bernd Fischer

Maintainer: Bernd Fischer <bernd.fischer@embl.de>
}
\references{
T. Horn, T. Sandmann, B. Fischer, W. Huber, M. Boutros. Mapping of Signalling Networks through Synthetic 
Genetic Interaction Analysis by RNAi. Nature Methods, 2011.
}
\keyword{ package }
\examples{
data(datamatrix, package="DmelSGI")
}



