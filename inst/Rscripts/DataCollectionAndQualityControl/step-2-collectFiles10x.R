
library(rhdf5)

host = system("echo $HOSTNAME",intern=TRUE)
if (host != "schroedinger.embl.de") {
  stop("This script has to be executed on schroedinger")
}

############################
## functions to compute run time
############################

startProcessing <- function(text) {
  time1 = Sys.time()
  cat(text, " ... ")
  if (nchar(text) < 25) {
    for (i in (nchar(text)+1):25) { cat(" ") }
  }
  return(time1)
}

endProcessing <- function(time1) {
  time2 = Sys.time()
  cat("finished. Time: ",format(difftime(time2, time1, tz = "",units = c("auto"))),".\n")
  return(invisible(NULL))
}

############################
## start of script, read in 4x data and location of data
############################

time0 = startProcessing("Start collection of features")
time1 = startProcessing("load libraries and screen annotation")

A = h5read("Anno", file="../h5/screenData4x.h5")

load("../COPY/diskcontent.rda")

Plates = A[[1]]
Plates = Plates[!duplicated(Plates)]

m = match(Plates,df[[1]])
df = df[m,]
endProcessing(time1)

## Batches = rep(1:78,each=16)

## for (b in 1:78) {
##   cat("batch",b,"\n")
## }

time2 = startProcessing("Allocate memory for matrix")
X = matrix(NA, nr=nrow(A), nc = 1507)
endProcessing(time2)

############################
## read features plate by plate
############################

time3 = startProcessing("Read features")
for (i in 1:length(Plates)) {
  time4 = startProcessing(sprintf("Read data from plate %d: %s",i,df[[1]][i]))
  if (!is.na(df[[1]][i])) {
    for (j in 1:16) {
      load(sprintf("../version2/%s/output/featureswell/%s-%s-1-1.rda",df[[2]][i],df[[1]][i],LETTERS[j]))
      I = 1:24+(j-1)*24 + (i-1)*384
      X[I,] = Ftable
    }
  }
  endProcessing(time4)
}
endProcessing(time3)

############################
## write 10x HDF5 file
############################

time3 = startProcessing("Write features")
h5write(X, name="Data", file="../h5/screenData10x.h5")
endProcessing(time3)

endProcessing(time0)

