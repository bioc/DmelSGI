
subSampleForStabilitySelectionFct <- function(fileMatrixData,
                                           N = 10000,
                                           random.seed=NULL) {
  if (!is.null(random.seed)) {
    set.seed(seed = random.seed)
  }
  DKD = h5read("DKD",file = fileMatrixData)

  n = 1
  I = phenotype = 5
  m = 3
  d = dim(DKD$D)

  D = DKD$D
  D = apply(D, 2:5, function(x) { x = (x - median(x,na.rm=TRUE)) / mad(x, na.rm=TRUE)  })

  D[!is.finite(D)] = 0.0
  dim(D) = c(prod(dim(D)[1:3]),2,dim(D)[5])
  S = sample(1:(dim(D)[1]), N)
  D = D[S,,]
  res <- list(D=D, Sample = S, phenotype = DKD$A$phenotype$phenotype)
  res
}

