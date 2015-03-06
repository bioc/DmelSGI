
library(rhdf5)

host = system("echo $HOSTNAME",intern=TRUE)
if (host != "schroedinger.embl.de") {
  stop("This script has to be executed on schroedinger")
}

if(file.exists("../h5/screenData.h5")) {
  stop("file 'screenData.h5' already exists in dir ../h5/. Remove file before running this script!")
}

##################################
## read a list of features that are taken into account for the interaction screen
##################################

F = readLines("../h5/selectedfeatures.txt")

if(file.exists("../h5/screenData.h5")) {
  print("file 'screenData.h5' already exists")
} else {

##################################
## read the 4x data and write the annotation to the new file
##################################

  X = h5read("/",file="../h5/screenData4x.h5")
  attributes(X$Anno)$row.names = NULL
  h5createFile("../h5/screenData.h5")
  h5write(X$Anno, file="../h5/screenData.h5", name="Anno")
  h5write(X$QueryLayout, file="../h5/screenData.h5", name="QueryLayout")
  h5write(X$TemplateLayout, file="../h5/screenData.h5", name="TemplateLayout")

##################################
## read the 10x data
##################################

  X10x = h5read("/",file="../h5/screenData10x.h5")

##################################
## select the features that are considered for the interaction screen
##################################

  I = which(X10x$Phenotype$Phenotype %in% F)
  Phenotype = data.frame(Phenotype = c(X$Phenotype$Phenotype,X10x$Phenotype$Phenotype[I]), stringsAsFactors=FALSE )
  h5write(Phenotype, file="../h5/screenData.h5", name="Phenotype")

##################################
## create a new matrix and writew the 4x and 10x data to the file
##################################

  L1 = nrow(X$Phenotype)
  L2 = length(I)
  
  d = c(dim(X$D)[1],nrow(Phenotype))
  h5createDataset(dataset="Data", dims=d, file="../h5/screenData.h5",storage.mode="double",level=0)
  h5write(X$D, name="Data", file="../h5/screenData.h5", index=list(NULL,1:L1))
  h5write(X10x$D[,I], name="Data", file="../h5/screenData.h5", index=list(NULL,1:L2+L1))
}

