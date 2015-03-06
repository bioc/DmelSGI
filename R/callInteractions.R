
callInteractions <- function(filePI, fileInteractions, verbose = TRUE, overwrite=FALSE) {
  if (file.exists(fileInteractions)) {
    if (overwrite) {
      if (verbose) { cat("file ",filePI, " is overwritten.\n") }
      file.remove(fileInteractions)
    } else {
      stop("file ",fileInteractions," already exists. Try with 'overwrite=TRUE'.")
    }
  }

  if (verbose) { cat("Read annotation\n") }
  Anno = h5read("DKD/A", file=filePI)
  P = nrow(Anno$phenotype)

  h5createFile(fileInteractions)
  h5createDataset(file=fileInteractions, dataset="piscore",
                  dims=c(sum(Anno$template$group == "sample"),nrow(Anno$query),nrow(Anno$phenotype)), 
                  storage.mode="double",level=0)
  h5createDataset(file=fileInteractions, dataset="padj",
                  dims=c(sum(Anno$template$group == "sample"),nrow(Anno$query),nrow(Anno$phenotype)), 
                  storage.mode="double",level=0)
  for (p in 1:P) {
    if (verbose) {
      cat("Phenotype ",p, " out of ",P,"\n")
    }
    Data = h5read("DKD/D", file=filePI, index=list(NULL,NULL,NULL,NULL,p))
    Data = Data[Anno$template$group == "sample",,,,,drop=FALSE]
    Data = aperm(Data, c(1,3,2,4,5))
    d = dim(Data)
    dim(Data) = c(prod(d[1:2]),prod(d[3:4]))
    fit = eBayes(lmFit(Data))
    padj = p.adjust(fit$p.value, method="BH")
    dim(padj) = d[1:2]
    piscore = fit$coefficients
    dim(piscore) = d[1:2]
    h5write(piscore, file=fileInteractions, name="piscore", index=list(NULL,NULL,p))
    h5write(padj, file=fileInteractions, name="padj", index=list(NULL,NULL,p))
  }
  Anno$template = Anno$template[Anno$template$group == "sample",]
  Anno$queryDesign = NULL
  Anno$templateDesign = NULL
  h5write(Anno, file=fileInteractions, name="Anno")

  invisible(NULL)
}

