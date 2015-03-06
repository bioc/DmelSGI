
library(rhdf5)

source("applyStabilitySelection.R")
source("estimatePairwiseInteractions.R")
source("stabilitySelection.R")

## subsample the data to reduce computational time

subSampleForStabilitySelection = subSampleForStabilitySelectionFct("../h5/matrixData.h5",
  N = 3000, random.seed = 740269056)
subSampleForStabilitySelection2 = subSampleForStabilitySelectionFct("../h5/matrixData.h5",
                                                                N = 3000, random.seed = NULL)
save(subSampleForStabilitySelection, file="subSampleForStabilitySelection.rda")
save(subSampleForStabilitySelection2, file="subSampleForStabilitySelection2.rda")

# stability selection

load("subSampleForStabilitySelection.rda")
stabilitySelection = selectByStability(subSampleForStabilitySelection,
                               preselect = c("4x.count","4x.ratioMitotic","10x.meanNonmitotic.cell.0.s.area"),
                               Rdim = 50, method="PCA",
                               verbose = TRUE)
save(stabilitySelection, file="stabilitySelection.rda")

# apply stability selection to data matrix

load("stabilitySelection.rda")
applyDimensionReduction(fileMatrixData = "../h5/matrixData.h5",
                          fileNew = "../h5/matrixDataStabilitySelection.h5",
                          selected = stabilitySelection$selected[ stabilitySelection$ratioPositive > 0.5 ] )
datamatrix = h5read("DKD", file="../h5/matrixDataStabilitySelection.h5")
save(datamatrix, file="datamatrix.rda")

SKDdata = h5read(file="/Volumes/befische/experiments/MichaelBoutros/ChromatinSet/ChromatinSet10x/h5/matrixDataStabilitySelection.h5",name="SKD")
save(SKDdata, file="SKDdata.rda")

