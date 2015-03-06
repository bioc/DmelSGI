\name{RohnEtAl}
\alias{RohnEtAl}
\docType{data}
\title{Fly RNAi phenotype data}
\description{This dataset is a list of RNAi Phenotypes captures in Drosophila S2R+ cells.}
\usage{data("RohnEtAl")}
\value{A list of phenotypes. The genes showinbg the respective phenotype are listet in the vector for each phenotype.}
\format{
  A data frame with 556 observations on the following 29 variables.
  \describe{
    \item{\code{Primer}}{a character vector}
    \item{\code{Computed.Target}}{a character vector}
    \item{\code{Symbol}}{a character vector}
    \item{\code{Decreased.cell.size}}{a numeric vector}
    \item{\code{Increased.cell.size}}{a numeric vector}
    \item{\code{Cell.shape.variable}}{a numeric vector}
    \item{\code{Cell.shape.round.or.non.adherent}}{a numeric vector}
    \item{\code{Cell.shape.processes.or.spiky.or.stretchy}}{a numeric vector}
    \item{\code{Disorganised.peripheral.actin}}{a numeric vector}
    \item{\code{Increased.number.of.actin.stress.fibres}}{a numeric vector}
    \item{\code{Increased.number.of.actin.puncta.or.dots}}{a numeric vector}
    \item{\code{Asymmetric.lamellae}}{a numeric vector}
    \item{\code{Decreased.level.of.actin}}{a numeric vector}
    \item{\code{Increased.level.of.actin}}{a numeric vector}
    \item{\code{Increased.cytoplasmic.actin}}{a numeric vector}
    \item{\code{Decreased.peripheral.actin}}{a numeric vector}
    \item{\code{Increased.peripheral.actin}}{a numeric vector}
    \item{\code{Increased.nuclear.actin}}{a numeric vector}
    \item{\code{Microtubule.clumps}}{a numeric vector}
    \item{\code{Microtubules.disorganised}}{a numeric vector}
    \item{\code{Microtubule.processes}}{a numeric vector}
    \item{\code{Decreased.level.of.microtubules}}{a numeric vector}
    \item{\code{Increased.level.of.microtubules}}{a numeric vector}
    \item{\code{Increased.number.of.multinucleate.cells}}{a numeric vector}
    \item{\code{Increased.DNA.area}}{a numeric vector}
    \item{\code{No.cells}}{a numeric vector}
    \item{\code{Decreased.cell.number}}{a numeric vector}
    \item{\code{Loss.of.cell.monolayer}}{a numeric vector}
    \item{\code{Multiple.layers.of.cells}}{a numeric vector}
  }
}
\seealso{\link{DmelSGI-package}}
\references{
Rohn, Jennifer L., et al. "Comparative RNAi screening identifies a conserved core metazoan actinome by phenotype." The Journal of cell biology 194.5 (2011): 789-805.
}
\examples{
data(RohnEtAl)
}
\keyword{datasets}
