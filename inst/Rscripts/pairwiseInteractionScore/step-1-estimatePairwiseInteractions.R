
library(ChromatinSet)

host = system("echo $HOSTNAME",intern=TRUE)
if (host != "schroedinger.embl.de") {
  stop("This script has to be executed on schroedinger")
}

load("stabilitySelection.rda")
## source("applyStabilitySelection.R")
## source("mymedpolish.R")
## source("estimatePairwiseInteractions.R")
## source("callInteractions.R")
fileMatrixData = "../h5/matrixData.h5"
fileMatrixPI = "../h5/PI.h5"
fileMatrixInteractions = "../h5/Interactions.h5"

estimatePairwiseInteractions(fileMatrixData, fileMatrixPI, useSKD=TRUE)
# estimatePairwiseInteractions(fileMatrixData, "../h5/PIfalse.h5", useSKD=FALSE)

applyDimensionReduction(fileMatrixData = "../h5/PI.h5",
                          fileNew = "../h5/PIStabilitySelection.h5",
                          selected = stabilitySelection$selected[ stabilitySelection$ratioPositive > 0.5 ] )

pimatrix = h5read("DKD", file="../h5/PIStabilitySelection.h5")
save(pimatrix, file="pimatrix.rda")
mainEffects = h5read("mainEffects", file="../h5/PIStabilitySelection.h5")
save(mainEffects, file="mainEffects.rda")

callInteractions(fileMatrixPI, fileMatrixInteractions)
selected = stabilitySelection$selected[ stabilitySelection$ratioPositive > 0.5 ]
Anno = h5read(file=fileMatrixInteractions,name="Anno")
I = match(selected,Anno$phenotype$phenotype)
piscore = h5read(file=fileMatrixInteractions, name="piscore",index=list(NULL,NULL,I))
padj = h5read(file=fileMatrixInteractions, name="padj",index=list(NULL,NULL,I))
Anno$phenotype = data.frame(phenotype = selected, stringsAsFactors=FALSE)
Interactions = list(piscore = piscore, padj = padj, Anno=Anno)
save(Interactions, file="Interactions.rda")
