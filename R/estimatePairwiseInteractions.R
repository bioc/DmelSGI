
estimatePairwiseInteractions <- function(fileMatrixData, filePI, verbose = TRUE, overwrite=FALSE, useSKD = TRUE) {
  if (file.exists(filePI)) {
    if (overwrite) {
#       warning("file ",filePI, " is overwritten")
    } else {
      stop("file ",filePI," already exists. Try with 'overwrite=TRUE'.")
    }
  }

  if (verbose) { cat("copy file ",fileMatrixData," to ", filePI,"\n") }
  file.copy(fileMatrixData, filePI, overwrite = overwrite)

  if (verbose) { cat("Read annotation\n") }
  Anno = h5read("DKD/A", file=fileMatrixData)

#   DKD, SKD
  TP = Anno$template$TemplatePlate

  ## Dimensions
  P = nrow(Anno$phenotype)
  B = max(Anno$query$Batch)
  TDS = max(Anno$templateDesign$design)
  QDS = max(Anno$templateDesign$design)
  MTP = max(TP)

  ## result arrays
  templateMainEffect = array(NA_real_,dim=c(nrow(Anno$template),nrow(Anno$templateDesign),B,P))
  queryMainEffect = array(NA_real_,dim=c(nrow(Anno$query),nrow(Anno$queryDesign),MTP,P))
  overallEffect = array(NA_real_,dim=c(nrow(Anno$templateDesign), B, P))

  for (p in 1:P) {
    if (verbose) {
      cat("Phenotype ",p, " out of ",P,"\n")
    }
    Data = h5read("DKD/D", file=fileMatrixData, index=list(NULL,NULL,NULL,NULL,p))
    PI = array(NA_real_,dim=dim(Data))
    DataSKD = h5read("SKD/D", file=fileMatrixData, index=list(NULL,NULL,NULL,p))
    PISKD = array(NA_real_,dim=dim(DataSKD))
    for (tds in 1:TDS) {
      for (b in 1:B) {
        IB = which(Anno$query$Batch == b)
        D = Data[,tds,IB,,1]
        dim(D) = c(dim(D)[1],prod(dim(D)[2:3]))
        D = cbind(D, DataSKD[,tds,b,1])
        if (useSKD) {
          QueryNeg = dim(D)[2]
        } else {
          QueryNeg = NULL
        }
        MD = mymedpolish(D, TP=TP, 
          TemplateNeg = which(Anno$template$group == "neg"), 
          QueryNeg = QueryNeg)
        PI[,tds,IB,,1] = MD$pi[,1:(QDS*length(IB))]
        PISKD[,tds,b,1] = MD$pi[,QDS*length(IB)+1]
        templateMainEffect[,tds,b,p] = MD$templateMainEffect
        queryMainEffect[IB,,,p] = MD$queryMainEffect[1:(QDS*length(IB)),]
        overallEffect[tds,b,p] = MD$neg
      }
    }
    h5write(PI, "DKD/D", file=filePI, index=list(NULL,NULL,NULL,NULL,p))
    h5write(PISKD, "SKD/D", file=filePI, index=list(NULL,NULL,NULL,p))
  }
  h5createGroup(file=filePI,group="mainEffects")
  h5write(templateMainEffect, "mainEffects/templateMainEffect", file=filePI)
  h5write(queryMainEffect, "mainEffects/queryMainEffect", file=filePI)
  h5write(overallEffect, "mainEffects/overallEffect", file=filePI)

  ## plateControlDKD = h5read(file=fileMatrixData, name="plateControlDKD")
  ## plateControlSKD = h5read(file=fileMatrixData, name="plateControlSKD")
  ## h5write(plateControlDKD, file=filePI, name="plateControlDKD")
  ## h5write(plateControlSKD, file=filePI, name="plateControlSKD")

  if (verbose) { cat("\n") }

  invisible(TRUE)
}



