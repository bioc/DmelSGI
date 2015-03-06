
applyDimensionReduction <- function(fileMatrixData, fileNew,
                                    selected, 
                                    verbose = TRUE, overwrite=FALSE) {

  if (verbose) { cat("create file ", fileNew,"\n") }
  if (file.exists(fileNew)) {
    if (overwrite) {
      file.remove(fileNew)
    } else {
      stop("file ",fileNew, " already exists. Try with 'overwrite=TRUE'")
    }
  }
  h5createFile(fileNew)

  L = h5ls(fileMatrixData)
  for (name in c("DKD","SKD", "plateControlDKD", "plateControlSKD")) {
    if (name %in% L$name[L$group == "/"]) {
      data = h5read(name,file = fileMatrixData)
      if (name %in% c("DKD","plateControlDKD")) {
        data$D = data$D[,,,,match(selected, data$A$phenotype$phenotype)]
      } else {
        data$D = data$D[,,,match(selected, data$A$phenotype$phenotype)]
      }
      data$A$phenotype = data.frame(phenotype = selected,stringsAsFactors=FALSE)
      h5createGroup(name, file=fileNew)
      h5write(data$A, name=sprintf("%s/A",name), file=fileNew)
      h5createDataset(file=fileNew, dataset=sprintf("%s/D",name), dims=dim(data$D), level=0)
      h5write(data$D, name=sprintf("%s/D",name), file=fileNew)
    }
  }

  if ("mainEffects" %in% L$name[L$group == "/"]) {
    A = h5read("DKD/A",file = fileMatrixData)
    data = h5read("mainEffects",file = fileMatrixData)
    data$templateMainEffect = data$templateMainEffect[,,,match(selected, A$phenotype$phenotype)]
    data$queryMainEffect = data$queryMainEffect[,,,match(selected, A$phenotype$phenotype)]
    data$overallEffect = data$overallEffect[,,match(selected, A$phenotype$phenotype)]
    h5write(data, name="mainEffects", file=fileNew)
  }

  invisible(NULL)
}

