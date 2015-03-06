
library(abind)
library(rhdf5)

host = system("echo $HOSTNAME",intern=TRUE)
if (host != "schroedinger.embl.de") {
  stop("This script has to be executed on schroedinger")
}

dir.create("step-5-QC")

fileMatrixDataRaw = "../h5/matrixDataRaw.h5"
fileMatrixData = "../h5/matrixData.h5"

if (file.exists(fileMatrixData)) {
  stop("file ",fileMatrixData," already exists. Remove file before running this script!\n")
}

DKD = h5read("DKD",file = fileMatrixDataRaw)
SKD = h5read("SKD",file = fileMatrixDataRaw)
plateControlDKD = h5read("plateControlDKD",file = fileMatrixDataRaw)
plateControlSKD = h5read("plateControlSKD",file = fileMatrixDataRaw)

C = rep(NA,dim(DKD$D)[5])
n1 = rep(NA,dim(DKD$D)[5])
n2 = rep(NA,dim(DKD$D)[5])
## log transform data
logtrafo <- function(x,c) {
  log2((x+sqrt(x^2+c^2))/2)
}
for (i in 1:dim(DKD$D)[5]) {
  cat("log-transform phenotype ",i, " out of ",dim(DKD$D)[5],"\r")
  a = as.vector(DKD$D[,,,,i])
  ## add a value that is small with respect to the dynamic range of the phenotype to avoid -inf
  m = quantile(DKD$D[,,,,i],probs=0.03,na.rm=TRUE)
  DKD$D[,,,,i] = logtrafo(DKD$D[,,,,i],m)
  SKD$D[,,,i] = logtrafo(SKD$D[,,,i],m)
  plateControlDKD$D[,,,,i] = logtrafo(plateControlDKD$D[,,,,i],m)
  plateControlSKD$D[,,,i] = logtrafo(plateControlSKD$D[,,,i],m)
}
cat("\n\n")

###############################
## compute correlation of phenotype
###############################
cat("compute number of finite values.\n")
rFinite = apply(is.finite(DKD$D),5,sum) / prod(dim(DKD$D)[1:4])
D = DKD$D
D = D[,1,,,] + D[,2,,,]
D = D / 2
dim(D) = c(prod(dim(D)[1:2]),dim(D)[3],dim(D)[4])
zz = 0
PhenoCor = apply(D,3,function(X) {
  cat("compute correlation phenotype ",zz," of ",dim(D)[3],"\r")
  a = X[,1]
  b = X[,2]
  F  = which(is.finite(a) & is.finite(b))
  zz <<- zz + 1
  cor(a[F],b[F])
})
cat("\n")
qualityControlFeature = list(correlation = PhenoCor,
  ratioFiniteValues = rFinite,
  passed = (PhenoCor >= 0.6) & (rFinite > 0.99))
qualityControlFeature$passed[is.na(qualityControlFeature$passed)] = FALSE
save(qualityControlFeature, file="step-5-QC/qualityControlFeature.rda")

## Select the phenotypes with correlation > 0.6
P = which((PhenoCor >= 0.6) & (rFinite > 0.99))

cat("\nnumber of selected phenotypes: ",length(P),"\n\n")


###############################
## compute correlation of dsRNA designs for template genes
###############################
## select template genes
# I = which(DKD$A$template$group == "sample")
I = 1:nrow(DKD$A$template)

D = DKD$D[I,,,,P]
for (i in 1:dim(D)[5]) {
  cat("normalize phenotype ",i, " out of ",dim(D)[5],"\r")
  D[,,,,i] = D[,,,,i] - median(D[,,,,i],na.rm=TRUE)
  D[,,,,i] = D[,,,,i] / mad(D[,,,,i],na.rm=TRUE)
}
cat("\n")
D = D[,,,1,] + D[,,,2,]
D = D / 2

GeneCor = rep(NA, length(I))
for (i in 1:length(I)) {
  cat("compute correlation gene ",i, " out of ",length(I),"\r")
  a = as.vector(D[I[i],1,,])
  b = as.vector(D[I[i],2,,])
  F  = which(is.finite(a) & is.finite(b))
  GeneCor[i] = cor(a[F],b[F])
}
cat("\n")
G = I[which((GeneCor >= 0.7) | (DKD$A$template$group[I] != "sample"))]
qualityControlGene = list(correlation = GeneCor,
  Annotation = DKD$A$template,
  passed = rep(FALSE, nrow(DKD$A$template)))
qualityControlGene$passed[G] = TRUE
save(qualityControlGene,file="step-5-QC/qualityControlGene.rda")

cat("\nnumber of selected genes: ",sum(DKD$A$template$group[G] == "sample"),"\n\n")

##############################################
## write data back
##############################################

cat("create new hdf5 file.\n")
h5createFile(fileMatrixData)

## DKD
cat("write DKD.\n")
DKD$D = DKD$D[G,,,,P]
DKD$A$phenotype = DKD$A$phenotype[P,,drop=FALSE]
DKD$A$template = DKD$A$template[G,]

h5createGroup(file=fileMatrixData, group="DKD")
h5write(DKD$A, file=fileMatrixData, name="DKD/A")
h5createDataset(file=fileMatrixData,dataset="DKD/D",
                dims=dim(DKD$D),
                storage.mode="double",level=0)
h5write(DKD$D, name="DKD/D",file=fileMatrixData)

## SKD
cat("write SKD.\n")
SKD$D = SKD$D[G,,,P]
SKD$A$phenotype = SKD$A$phenotype[P,,drop=FALSE]
SKD$A$template = SKD$A$template[G,]

h5createGroup(file=fileMatrixData, group="SKD")
h5write(SKD$A, file=fileMatrixData, name="SKD/A")
h5createDataset(file=fileMatrixData,dataset="SKD/D",
                dims=dim(SKD$D),
                storage.mode="double",level=0)
h5write(SKD$D, name="SKD/D",file=fileMatrixData)

## plateControlDKD
cat("write plateControlDKD.\n")
plateControlDKD$D = plateControlDKD$D[,,,,P]
plateControlDKD$A$phenotype = plateControlDKD$A$phenotype[P,,drop=FALSE]

h5createGroup(file=fileMatrixData, group="plateControlDKD")
h5write(plateControlDKD$A, file=fileMatrixData, name="plateControlDKD/A")
h5createDataset(file=fileMatrixData,dataset="plateControlDKD/D",
                dims=dim(plateControlDKD$D),
                storage.mode="double",level=0)
h5write(plateControlDKD$D, name="plateControlDKD/D",file=fileMatrixData)

## plateControlSKD
cat("write plateControlSKD.\n")
plateControlSKD$D = plateControlSKD$D[,,,P]
plateControlSKD$A$phenotype = plateControlSKD$A$phenotype[P,,drop=FALSE]

h5createGroup(file=fileMatrixData, group="plateControlSKD")
h5write(plateControlSKD$A, file=fileMatrixData, name="plateControlSKD/A")
h5createDataset(file=fileMatrixData,dataset="plateControlSKD/D",
                dims=dim(plateControlSKD$D),
                storage.mode="double",level=0)
h5write(plateControlSKD$D, name="plateControlSKD/D",file=fileMatrixData)

cat("Ready.\n")
