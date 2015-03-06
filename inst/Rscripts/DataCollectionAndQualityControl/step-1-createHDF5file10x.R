
library(rhdf5)

host = system("echo $HOSTNAME",intern=TRUE)
if (host != "schroedinger.embl.de") {
  stop("This script has to be executed on schroedinger")
}

if(file.exists("../h5/screenData10x.h5")) {
  stop("file 'screenData10x.h5' already exists in dir ../h5/. Remove file before running this script!")
} else {
  load("../version2/disk44/output/featureswell/A2004BD063-A-1-1.rda")
  P10x = paste("10x",colnames(Ftable),sep=".")

  ################################
  ## read 4x data and write annotation to 10x HDF5 file
  ################################
  
  X = h5read("/",file="../h5/screenData4x.h5")
  attributes(X$Anno)$row.names = NULL
  h5createFile("../h5/screenData10x.h5")
  h5write(X$Anno, file="../h5/screenData10x.h5", name="Anno")
  h5write(X$QueryLayout, file="../h5/screenData10x.h5", name="QueryLayout")
  h5write(X$TemplateLayout, file="../h5/screenData10x.h5", name="TemplateLayout")
  Phenotype = data.frame(Phenotype = P10x, stringsAsFactors=FALSE)
  h5write(Phenotype, file="../h5/screenData10x.h5", name="Phenotype")

  ################################
  ## create a new data matrix in the file
  ################################
  
  A = h5read("Anno", file="../h5/screenData4x.h5")
  d = c(dim(A)[1],nrow(Phenotype))
  h5createDataset(dataset="Data", dims=d, file="../h5/screenData10x.h5",storage.mode="double",level=0)
}

