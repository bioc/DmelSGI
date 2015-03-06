library(EBImage)
library(MASS)
library(sda)
## library(grid)
## library(RColorBrewer)
## library(genefilter)
library(hwriter)
## library(e1071)

zernikeRadialPolynomial <- function(rho, n, m) {
  x <- 0.0
  if ((n - m) %% 2 == 0) {
    for (s in 0:((n-abs(m))/2)) {
      x <-x + rho^(n-2*s) * ((-1)^s*factorial(n-s)) / (factorial(s)*factorial((n+abs(m))/2-s)*factorial((n-abs(m))/2-s))
    }
  } else {
    x <- 0.0
  }
  x[rho > 1] <- 0.0
  x
}

zernikePolynomial <- function(n, m, R=30) {
  r <- ceiling(R)
  D <- 2*r+1
  a <- seq(-r,r,by=1)
  M <- array(0.0,dim=c(D,D,2))
  X <- matrix(a,nr=D,nc=D)
  Y <- t(X)
  rho <- sqrt(X^2+Y^2) / R
  G = Y / (R*rho)
  G[is.na(G)] = -1
  G[G > 1.0] = 1.0
  G[G < -1.0] = -1.0
  phi <- acos(G)
  phi[X < 0] <-  2*pi-phi[X < 0]
  M[,,1] <- zernikeRadialPolynomial(rho, n, m)
  M[,,2] <- M[,,1]
  M[,,1] <- M[,,1] * cos(m * phi)
  M[,,2] <- M[,,2] * sin(-m * phi)
  M[is.na(M)] = 0.0
  M
}

zernikeBasis <- function(N=12, R=30) {
  r <- ceiling(R)
  D <- 2*r+1
  zBasis <- list()
  z <- 0
  for (n in 0:N) {
    for (m in rev(seq(n,0,by=-2))) {
      z=z+1
      zBasis[[z]] = zernikePolynomial(n,m,R)
      names(zBasis)[z] = sprintf("z.%d.%d",n,m)
    }
  }
  zBasis
}

basisTransformation <- function (x, ref, filter, m=NULL) {
  validObject(ref)
  validObject(filter)
  if (length(dim(ref)) > 2) {
    stop("this method doesn't support multi channel images")
  }
  if (colorMode(ref) == TrueColor) {
    stop("this method doesn't support the 'TrueColor' color mode")
  }
  dx = dim(ref)
  cmx = colorMode(ref)
  df = dim(filter)[1:2]
  if (any(df%%2 == 0)) {
    stop("dimensions of 'filter' matrix must be odd")
  }
  if (any(dx[1:2] < df)) {
    stop("dimensions of 'ref' must be bigger than 'filter'")
  }
  if (is.null(m)) {
    m = round(hullFeatures(x)[,c("g.x","g.y")])
  }
  A <- matrix(1:prod(dx),nr=dx[1],nc=dx[2])
  I <-A[m]
  rm(A)

  cx = dx%/%2
  cf = df%/%2

  fftx <- fft(ref)
  y <- apply(filter, 3, function(ff) {
    wf = matrix(0, nr = dx[1], nc = dx[2])
    wf[(cx[1] - cf[1]):(cx[1] + cf[1]), (cx[2] - cf[2]):(cx[2] + cf[2])] = ff
    wf = fft(wf)
    ## dim(x) = c(dx[1:2], prod(dx)/prod(dx[1:2]))
    index1 = c(cx[1]:dx[1], 1:(cx[1] - 1))
    index2 = c(cx[2]:dx[2], 1:(cx[2] - 1))
    pdx = prod(dim(ref)[1:2])
    (Re(fft(fftx * wf, inverse = T)/pdx)[index1, index2])[I]
  })
  dimnames(y) <- dimnames(filter)[[3]]
  ## dim(y) = c(dx[1:2],df[3])
  ## if (is.Image(x)) {
  ##   return(Image(y, colormode = cmx))
  ## }  else {
    return(y)
  ## }
}

extractPatches <- function (pos, ref, df, normalize = TRUE, combine = TRUE) {
  validObject(ref)
  if (length(dim(ref)) > 2) {
    stop("this method doesn't support multi channel images")
  }
  if (colorMode(ref) == TrueColor) {
    stop("this method doesn't support the 'TrueColor' color mode")
  }

  if (mode(ref) != "list") {
    ref <- list(ref)
  }
  if (is.null(names(ref))) {
    names(ref) <- sprintf("ch%d", seq_len(length(ref)))
  }
  dx = dim(ref[[1]])[1:2]

  Coord <- matrix(1:prod(dx),nr=dx[1],nc=dx[2])
  I <- Coord[pos]
  J <- Coord[1:df[1],1:df[2]]-1
  J <- J - J[round(median(seq_len(df[1]))),round(median(seq_len(df[2])))]
  I <- rep(I, each = length(J))
  J <- rep(J, times = nrow(pos))
  K <- I + J
  A <- list()
  for (i in seq_len(length(ref))) {
    A[[i]] <- ref[[i]][K]
    dim(A[[i]]) = c(prod(df),nrow(pos))

    if (normalize) {
      m <- apply(A[[i]], 2, median, na.rm=TRUE)
      A[[i]] <- A[[i]] - array(rep(m, each=prod(df)), dim = dim(A[[i]]))
      s <- apply(A[[i]], 2, mad, na.rm=TRUE)
      s[s < 1.0e-100] <- 1.0
      A[[i]] <- A[[i]] / array(rep(s, each=prod(df)), dim = dim(A[[i]]))
    }
    dim(A[[i]]) = c(df,nrow(pos))
  }
  names(A) <- names(ref)
  if (combine) {
    n <- length(A)
    Names <- names(A)
    z <- length(A)
    if (n > 1) {
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          z <- z + 1
          A[[z]] <- A[[i]] * A[[j]]
          names(A)[z] <- sprintf("%s%s", Names[i], Names[j])
        }
      }
    }
  }
  A
}

basisTransformation3 <- function (A, filterbank) {

  nchannel <- length(A)
  namesChannel <- names(A)
  if (is.null(namesChannel)) {
    namesChannel <- sprintf("ch%d", seq_len(nchannel))
  }

  nfilter <- length(filterbank)
  namesFilter <- names(filterbank)
  if (is.null(namesFilter)) {
    namesFilter <- sprintf("f%d", seq_len(nfilter))
  }

  F <- matrix(0.0, nr=dim(A[[1]])[3], nc = nchannel * nfilter)
  colnames(F) <- sprintf("%s",1:ncol(F))
  z <- 0
  for (i in 1:nchannel) {
    A1 = A[[i]]
    dim(A1) <- c(prod(dim(A[[i]])[1:2]), dim(A[[i]])[3])
    A1 = t(A1)
    for (j in 1:nfilter) {
      if (length(dim(filterbank[[j]])) > 2) {
        filter <-matrix(0.0, nr=prod(dim(filterbank[[j]])[1:2]), nc=dim(filterbank[[j]])[3])
        for (k in seq_len(dim(filterbank[[j]])[3])) {
          filter[,k] <- filterbank[[j]][,,k]
        }
        S <- A1 %*% filter
        S <- S^2
        S <- sqrt(apply(S, 1, sum))
        dim(S) <- c(dim(A1)[1],1)
      } else {
        filter <- filterbank[[j]]
        dim(filter) <- c(prod(dim(filter)),1)
        S <- A1 %*% filterbank[[j]]
      }
      z <- z + 1
      F[,z] <- S
      colnames(F)[z] <- sprintf("%s.%s", namesChannel[i], namesFilter[j])
    }
  }
}

basisTransformation2 <- function (x, ref, filterbank, m=NULL) {
  validObject(ref)
  validObject(filterbank)
  if (length(dim(ref)) > 2) {
    stop("this method doesn't support multi channel images")
  }
  if (colorMode(ref) == TrueColor) {
    stop("this method doesn't support the 'TrueColor' color mode")
  }
  if (mode(filterbank) != list) {
    filterbank <- list(filterbank)    
  }
  nfilter <- length(filterbank)
  namesFilter <- names(filterbank)
  if (is.null(namesFilter)) {
    namesFilter <- sprintf("f%d", seq_len(nfilter))
  }

  dx = dim(ref)
  cmx = colorMode(ref)
  df = dim(filter)[1:2]
  if (any(df%%2 == 0)) {
    stop("dimensions of 'filter' matrix must be odd")
  }
  if (any(dx[1:2] < df)) {
    stop("dimensions of 'ref' must be bigger than 'filter'")
  }
  if (is.null(m)) {
    m = round(hullFeatures(x)[,c("g.x","g.y")])
  }
  A <- matrix(1:prod(dx),nr=dx[1],nc=dx[2])
  I <- A[m]
  J <- A[1:df[1],1:df[2]]-1
  I <- rep(I, each = length(J))
  J <- rep(J, times = nrow(m))
  K <- I + J
  A <- ref[K]
  A[is.na(A)] = 0.0
  ## rm(A)
  ## A <- matrix(0.0, nr=length(I), nc=prod(df[1:2]))
  ## for (i in seq_len(length(I))) {
  ##   print(i)
  ##   A[i,] <- ref[J + I[i]]
  ## }
  dim(filter) = c(prod(dim(filter)[1:2]),dim(filter)[3])
  dim(A) = c(nrow(m),prod(dim(filter)[1]))
  B = A %*% filter
  dimnames(B)[[2]] <- Names
  B
}

mygetFeatures <- function (x, ref, m, nc = 32, N=9, R=15, filter) {
  ## validImage(x)
  if (colorMode(ref) == TrueColor) {
    stop("'ref' must not be in 'TrueColor' color mode")
  }
  gf = hullFeatures(x)
  mf = moments(x = x, ref = ref)
  ef = edgeFeatures(x = x, ref = ref)
  hf = haralickFeatures(x = x, ref = ref, nc = nc)
  z = basisTransformation2(x = x, ref = ref, filter = zernikeBasis(N=N, R=R), m = m)
  ei = basisTransformation2(x = x, ref = ref, filter = filter, m = m)
  features = cbind(gf, mf, ef, hf, z, ei)
  
  ## if (getNumberOfFrames(x, "total") == 1) {
  ##   impxs = match("m.pxs", colnames(mf))
  ##   mf = mf[, -impxs, drop = FALSE]
  ## } else {
  ##   impxs = match("m.pxs", colnames(mf[[1]]))
  ##   for (i in seq_along(hf)) {
  ##     mf[[i]] = mf[[i]][, -impxs, drop = FALSE]
  ##   }
  ## }
  ## if (getNumberOfFrames(x, "total") == 1) {
  ##   features = list(cbind(gf, mf, ef, hf))
  ## } else {
  ##   features = vector("list", length(gf))
  ##   for (i in seq_along(gf)) {
  ##     features[[i]] = cbind(gf[[i]], mf[[i]], ef[[i]], hf[[i]])
  ##   }
  ## }
  return(features)
}

mynormalize <- function(Img, minInt, maxInt) {
  Img[Img > maxInt] = maxInt
  Img[] = Img[] - minInt
  Img[Img < 0] = 0
  Img[] = Img[] / (maxInt - minInt)
  Img[Img > 1] = 1
  Img
}

correctFlatField <- function(Img, flatFieldModel) {
  imageData(Img) = flatFieldModel$A + flatFieldModel$B * imageData(Img)
  return(Img)
}

permLabel <- function(x, normalize = TRUE) {
  M <- max(x)
  R <- sample(M)
  x[x > 0] <- R[x[x>0]]
  if (normalize) {
    x <- normalize(x)
  }
  x
}

getGalleryImages <- function(I, x, y, Img1, Img2, nucLI, cytLI, minCh1, maxCh1, minCh2, maxCh2) {
  ## GImg = array(0.0, dim=c(length(I), 2, 81, 81, 3))
  GImg <- list()
  for (i in 1:length(I)) {
    cat("gallery ", i, "\n")
    cropX = (x[i]-40):(x[i]+40)
    cropY = (y[i]-40):(y[i]+40)
    ## nucLI2 = nucLI[cropX,cropY]
    ## nucLI2[nucLI2 != I[i]] = 0
    cytLI2 = cytLI[cropX,cropY]
    cytLI2[cytLI2 != I[i]] = 0

    Img1c = normalize(imageData(Img1)[cropX,cropY])
    Img2c = normalize(imageData(Img2)[cropX,cropY])
    ## Img1c = mynormalize(imageData(Img1)[cropX,cropY],minCh1, maxCh1)
    ## Img2c = mynormalize(imageData(Img2)[cropX,cropY],minCh2, maxCh2)

    A = array(0,dim=c(dim(Img1c),3))
    A[,,3] = Img1c
    A[,,2] = Img2c
    Imgout1 = Image(A,colormode="Color")
    ## GImg[i,1,,,] = imageData(Imgout1)
    Imgout2 = paintObjects(cytLI2, Imgout1, col='#ff0000')
    ## Imgout2 = paintObjects(cytLI2, paintObjects(nucLI2, Imgout1, col='#ff00ff'), col='#ff0000')
    ## GImg[i,2,,,] = imageData(Imgout2)
    GImg[[i]] = Imgout2
    ## Imgout2 = paintObjects(cytLI2[cropX,cropY], Imgout1, col='#ff0000')
    ## GImg[i,3,,,] = imageData(Imgout2)

    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img1c
    ## A[,,2] = Img1c
    ## A[,,3] = Img1c
    ## Imgout1 = Image(A,colormode="Color")
    ## GImg[i,4,,,] = imageData(Imgout1)
    ## Imgout2 = paintObjects(nucLI2[cropX,cropY], Imgout1, col='#ff00ff')
    ## GImg[i,5,,,] = imageData(Imgout2)
    ## Imgout2 = paintObjects(cytLI2[cropX,cropY], Imgout1, col='#ff0000')
    ## GImg[i,6,,,] = imageData(Imgout2)

    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img2c
    ## A[,,2] = Img2c
    ## A[,,3] = Img2c
    ## Imgout1 = Image(A,colormode="Color")
    ## GImg[i,7,,,] = imageData(Imgout1)
    ## Imgout2 = paintObjects(nucLI2[cropX,cropY], Imgout1, col='#ff00ff')
    ## GImg[i,8,,,] = imageData(Imgout2)
    ## Imgout2 = paintObjects(cytLI2[cropX,cropY], Imgout1, col='#ff0000')
    ## GImg[i,9,,,] = imageData(Imgout2)
  }
  return(GImg)
}

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



getNearestNeighbor2D <- function(xs,ys,xd,yd, maxDist = 10, bandwidth = maxDist) {

  X1 = matrix(c(xs,ys),nr=length(xs),nc=2)
  X2 = matrix(c(xd,yd),nr=length(xd),nc=2)
  
  rx =range(xs)
  ry =range(ys)

  d = c(min(c(xs,xd)-1.0-bandwidth),min(c(xs,xd)-1.0),seq(rx[1]+bandwidth,rx[2],by=bandwidth),max(c(xs,xd)+1.0),max(c(xs,xd)+1.0+bandwidth))
  mx = length(d)-1
  Xs = cut(xs,d,label=FALSE) 
  Xd = cut(xd,d,label=FALSE)

  d = c(min(c(ys,yd)-1.0-bandwidth),min(c(ys,yd)-1.0),seq(ry[1]+bandwidth,ry[2],by=bandwidth),max(c(ys,yd)+1.0),max(c(ys,yd)+1.0+bandwidth))
  my = length(d)-1
  Ys = cut(ys,d,label=FALSE) 
  Yd = cut(yd,d,label=FALSE) 

  ID = split(1:length(xd),factor(Xd,levels=1:mx))
  ID = lapply(ID, function(x) { split(x,factor(Yd[x],levels=1:my)) } )
  IS = split(1:length(xs),factor(Xs,levels=1:mx))
  IS = lapply(IS, function(x) { split(x,factor(Ys[x],levels=1:my)) } )

  index = rep(-1,length(xs))
  dist = rep(-1,length(xs))
  for (i in 2:(length(IS)-1)) {
    for (j in 2:(length(IS[[i]])-1)) {
      if (length(IS[[i]][[j]]) > 0) {
        ## cat("i=",i," j=",j," L=",length(IS[[i]][[j]]),"\n")
        I = c(ID[[i-1]][[j-1]],ID[[i-1]][[j]],ID[[i-1]][[j+1]],ID[[i]][[j-1]],ID[[i]][[j]],ID[[i]][[j+1]],ID[[i+1]][[j-1]],ID[[i+1]][[j]],ID[[i+1]][[j+1]])
        if (length(I) > 0) {
          A = X1[IS[[i]][[j]],,drop=FALSE]
          B = X2[I,,drop=FALSE]
          D = sqrt(matrix(A[,1]^2+A[,2]^2,nr=length(IS[[i]][[j]]),nc=length(I)) + t(matrix(B[,1]^2+B[,2]^2,nr=length(I),nc=length(IS[[i]][[j]]))) - 2 * A %*% t(B))
          J = apply(D,1,which.min)
          md = apply(D,1,min)
          J[md > maxDist] = NA
          md[md > maxDist] = NA
          index[IS[[i]][[j]]] = I[J]
          dist[IS[[i]][[j]]] = md
        } else {
          index[IS[[i]][[j]]] = NA
          dist[IS[[i]][[j]]] = NA
        }
      }
    }
  }

  res = list()
  res$index = index
  res$dist = dist

  xs = xs[is.finite(index)]
  ys = ys[is.finite(index)]
  xd = xd[index[is.finite(index)]]
  yd = yd[index[is.finite(index)]]

  xdpred = rlm(xd-xs ~ xs)$fitted + xs
  ydpred = rlm(yd-ys ~ ys)$fitted + ys
  res$disttrafo = rep(NA,length(index))
  res$disttrafoZ = rep(NA,length(index))
  res$disttrafo[is.finite(index)] = sqrt((xd - xdpred)^2 + (yd - ydpred)^2)
  res$disttrafoZ = res$disttrafo / mad(res$disttrafo,center=0.0,na.rm=TRUE)

  return(res)
}

plotTransformationField <- function(xs,ys,xd,yd,nn, mag=20, maxDistZ = NULL) {
  grid.newpage()
  pushViewport(dataViewport(c(xs,xd),c(ys,yd)))
  col = rep("red",length(xs))
  if (is.null(maxDistZ)) {
    md = 1.05 * max(nn$dist,na.rm=TRUE)
    I = which(!is.na(nn$index))
    col[I] = gray(nn$dist[I] / md)
  } else {
    I = which(!is.na(nn$disttrafoZ))
    I = I[which(nn$disttrafoZ[I] <= maxDistZ)]
    col[I] = "black"
  }
  grid.points(unit(xs,"native"),unit(ys,"native"), pch=20, gp = gpar(col=col, cex=0.4))
  grid.polyline(unit(c(xs[I],xs[I]+mag*(xd[nn$index[I]]-xs[I])),"native"),unit(c(ys[I],ys[I]+mag*(yd[nn$index[I]]-ys[I])),"native"), id = c(1:length(I),1:length(I)), gp = gpar(col=col[I]))
  popViewport()
}

localDensity <- function(x,y,bandwidth=1, scales=1:8, return.filter2D = FALSE, return.filtered = FALSE, return.mean = TRUE) {
  Img = matrix(0,nr=2048,nc=2048)
  I = matrix(0,nr=length(x),nc=2)
  I[,1] = x
  I[,2] = y
  Img[I] = 1

  Filter = list()
  for (s in scales) {
    print(s)
    Filter[[s]] = list()
    SD = bandwidth*2^(s-1)
    w = ceiling(4.0*SD)
    Filter[[s]]$bandwidth = SD
    win = dnorm(-w:w,sd=SD)
    win = win / sum(win)
    Win = matrix(win, nr=length(win), nc=length(win))
    Filter[[s]]$filter1D = win
    Filter[[s]]$filter2D = Win * t(Win)
    Filter[[s]]$filtered = filter2(Img, Filter[[s]]$filter2D)
    Filter[[s]]$LCD = Filter[[s]]$filtered[I]
    if (!return.filter2D) {
      Filter[[s]]$filter2D = NULL
    }
    if (!return.filtered) {
      Filter[[s]]$filtered = NULL
    }
  }
  if (return.mean) {
    res = c()
    for (s in scales) {
      res = c(res, mean(Filter[[s]]$LCD, na.rm=TRUE))
    }
    names(res) = sprintf("LCD%d",scales)
    res = c(res, LCDratio = res[4] / res[8])
  }
  Filter$mean = res
  Filter
}

process.images.4x <- function(fnch1, fnch2,
                              displayImages = FALSE,
                              fnFlatFieldModel = "/g/huber/users/befische/experiments/MichaelBoutros/ChromatinSet/flatFieldCorrection/flatFieldModel.rda",
                              fnControlImages = NULL,
                              fnFeatures = NULL,
                              fnFeaturesWell = NULL,
                              fnLabelImages = NULL,
                              fnGalleryImages = NULL,
                              cropX = 940:1059,
                              cropY = 950:1049,
                              minCh1 = 0.0026,
                              maxCh1 = 0.01,
                              minCh2 = 0.0017,
                              maxCh2 = 0.005 ) {
  time1 = Sys.time()
  cat("Start image processing at ", format(time1, "%a %b %d %X %Y"), "\n")

  time3 = startProcessing("Read images")
  Img1 = readImage(fnch1)
  Img2 = readImage(fnch2)
  endProcessing(time3)

  ## time3 = startProcessing("Flat field correction")
  ## load(fnFlatFieldModel)
##   load("B2.rda")
##   ff2$B = B2
  ## Img1 = correctFlatField(Img1,ff1)
##   ff2$B = ff2$B / 2
  ## Img2 = correctFlatField(Img2,ff2)
  ## endProcessing(time3)

  time3 = startProcessing("Smooth images")
  ## Img1smooth = Img1
  Img1smooth = gblur(Img1,s=1)
  Img2smooth = gblur(Img2,s=1)
  endProcessing(time3)

  time3 = startProcessing("Segmentation")
  nuct = thresh(Img1smooth, 3.5, 3.5, 0.0003)
  nucLI = bwlabel(nuct)
  cytt = thresh(Img2smooth, 3.5, 3.5, 0.0003)
  cytLI = bwlabel(cytt)

  ## minCh2 = 0.0022
  ## maxCh2 = 3 * (median(Img2[cytt > 0]) - minCh2) + minCh2
  endProcessing(time3)

  time3 = startProcessing("Feature Extraction")
  Y = matrix(1:2048,nr=2048,nc=2048)
  X = t(Y)

  xn = tapply(X[nucLI > 0], nucLI[nucLI > 0], mean)
  yn = tapply(Y[nucLI > 0], nucLI[nucLI > 0], mean)
  xt = tapply(X[cytLI > 0], cytLI[cytLI > 0], mean)
  yt = tapply(Y[cytLI > 0], cytLI[cytLI > 0], mean)
  Fnuc = matrix(0.0,nr=length(xn),nc=10)
  FH3p = matrix(0.0,nr=length(xt),nc=11)
  colnames(Fnuc) = c("x","y","area","meanint","maxint", "int", "meanint2","maxint2", "int2", "isMitotic")
  colnames(FH3p) = c("x","y","area","meanint","maxint", "int", "meanint1","maxint1", "int1", "nucindex", "dist")

  ## save(xt,yt,xn,yn,Img2,file="tmp.rda")
  nn <- getNearestNeighbor2D(xt,yt,xn,yn, maxDist=5, bandwidth = 10)
  nn2 <- getNearestNeighbor2D(xn,yn,xt,yt, maxDist=5, bandwidth = 10)
  Fnuc[,"x"] = xn
  Fnuc[,"y"] = yn
  FH3p[,"x"] = xt
  FH3p[,"y"] = yt
  FH3p[,"nucindex"] = nn$index
  FH3p[,"dist"] = nn$dist
  FH3p[,"area"] = tapply((matrix(1,nr=2048,nc=2048))[cytLI > 0], cytLI[cytLI > 0], sum)
  FH3p[,"maxint"] = tapply(imageData(Img2)[cytLI > 0], cytLI[cytLI > 0], max)
  FH3p[,"meanint"] = tapply(imageData(Img2)[cytLI > 0], cytLI[cytLI > 0], mean)
  FH3p[,"int"] = tapply(imageData(Img2)[cytLI > 0], cytLI[cytLI > 0], sum)
  Fnuc[,"area"] = tapply((matrix(1,nr=2048,nc=2048))[nucLI > 0], nucLI[nucLI > 0], sum)
  Fnuc[,"maxint"] = tapply(imageData(Img1)[nucLI > 0], nucLI[nucLI > 0], max)
  Fnuc[,"meanint"] = tapply(imageData(Img1)[nucLI > 0], nucLI[nucLI > 0], mean)
  Fnuc[,"int"] = tapply(imageData(Img1)[nucLI > 0], nucLI[nucLI > 0], sum)
  Fnuc[,"maxint2"] = tapply(imageData(Img2)[nucLI > 0], nucLI[nucLI > 0], max)
  Fnuc[,"meanint2"] = tapply(imageData(Img2)[nucLI > 0], nucLI[nucLI > 0], mean)
  Fnuc[,"int2"] = tapply(imageData(Img2)[nucLI > 0], nucLI[nucLI > 0], sum)
  FH3p[,"maxint1"] = Fnuc[nn$index,"maxint"]
  FH3p[,"meanint1"] = Fnuc[nn$index,"meanint"]
  FH3p[,"int1"] = Fnuc[nn$index,"int"]
  I = unique(nn$index[is.finite(nn$index)])
  Fnuc[I,"isMitotic"] = 1
  
  dnuc = sqrt((Fnuc[,"x"] - 1024)^2 + (Fnuc[,"y"] - 1024)^2)
  dH3p = sqrt((FH3p[,"x"] - 1024)^2 + (FH3p[,"y"] - 1024)^2)
  ## Inuc = dnuc <= 900
  Inuc = 1:length(dnuc)
  ## IH3p = (dH3p <= 900) & (nn$disttrafoZ <= 5)
  IH3p = (nn$disttrafoZ <= 5)
  LCD = localDensity(Fnuc[,"x"], Fnuc[,"y"])
  Fwell = c(count=sum(Inuc,na.rm=TRUE),
    countpH3 = sum(IH3p,na.rm=TRUE),
    isMitotic = sum(Fnuc[Inuc,"isMitotic"],na.rm=TRUE),
    ratioMitotic = sum(Fnuc[Inuc,"isMitotic"],na.rm=TRUE) / sum(Inuc,na.rm=TRUE),
    areaNuc = mean(Fnuc[,"area"],na.rm=TRUE),
    areapH3 = mean(FH3p[,"area"],na.rm=TRUE),
    intNuc = mean(Fnuc[,"meanint"],na.rm=TRUE),
    intpH3 = mean(FH3p[,"meanint"],na.rm=TRUE),
    q20intNuc = quantile(Fnuc[,"meanint"],probs=0.2,na.rm=TRUE)[[1]],
    q20intpH3 = quantile(FH3p[,"meanint"],probs=0.2,na.rm=TRUE)[[1]],
    q80intNuc = quantile(Fnuc[,"meanint"],probs=0.8,na.rm=TRUE)[[1]],
    q80intpH3 = quantile(FH3p[,"meanint"],probs=0.8,na.rm=TRUE)[[1]],
    LCD$mean
    )
  endProcessing(time3)

  if (!is.null(fnFeatures)) {
    time3 = startProcessing("Save features")
    save(Fnuc,FH3p,nn,nn2,LCD, file=fnFeatures)
    endProcessing(time3)
  }
  if (!is.null(fnFeaturesWell)) {
    time3 = startProcessing("Save features per well")
    save(Fwell,file=fnFeaturesWell)
    endProcessing(time3)
  }
  if (!is.null(fnLabelImages)) {
    time3 = startProcessing("Save label images")
    save(nucLI, cytLI, file=fnLabelImages)
    endProcessing(time3)
  }

  if (!is.null(fnControlImages)) {
    time3 = startProcessing("Save control images")

    writeImage(permLabel(nucLI)[cropX,cropY], sprintf("%s-li-A-seg.png", fnControlImages),quality=100)
    writeImage(permLabel(cytLI)[cropX,cropY], sprintf("%s-li-B-seg.png", fnControlImages),quality=100)
    Img1c = mynormalize(imageData(Img1),minCh1, maxCh1)[cropX,cropY]
    Img2c = mynormalize(imageData(Img2),minCh2, maxCh2)[cropX,cropY]

    A = array(0,dim=c(dim(Img1c),3))
    A[,,3] = 1.5*Img1c
    A[A > 1] = 1
    A[,,2] = Img2c
    Imgout1 = Image(A,colormode="Color")
    writeImage(Imgout1,sprintf("%s-color-A-img.png", fnControlImages),quality=98)
    Imgout2 = paintObjects(nucLI[cropX,cropY], Imgout1, col='#ff00ff')
    writeImage(Imgout2,sprintf("%s-color-B-nuc.png", fnControlImages),quality=98)
    Imgout2 = paintObjects(cytLI[cropX,cropY], Imgout1, col='#ff0000')
    writeImage(Imgout2,sprintf("%s-color-C-cyt.png", fnControlImages),quality=98)

    A = array(0,dim=c(dim(Img1c),3))
    A[,,1] = Img1c
    A[,,2] = Img1c
    A[,,3] = Img1c
    Imgout1 = Image(A,colormode="Color")
    writeImage(Imgout1,sprintf("%s-ch1-A-img.png", fnControlImages),quality=98)
    Imgout2 = paintObjects(nucLI[cropX,cropY], Imgout1, col='#ff00ff')
    writeImage(Imgout2,sprintf("%s-ch1-B-nuc.png", fnControlImages),quality=98)
    Imgout2 = paintObjects(cytLI[cropX,cropY], Imgout1, col='#ff0000')
    writeImage(Imgout2,sprintf("%s-ch1-C-cyt.png", fnControlImages),quality=98)

    A = array(0,dim=c(dim(Img1c),3))
    A[,,1] = Img2c
    A[,,2] = Img2c
    A[,,3] = Img2c
    Imgout1 = Image(A,colormode="Color")
    writeImage(Imgout1,sprintf("%s-ch2-A-img.png", fnControlImages),quality=98)
    Imgout2 = paintObjects(nucLI[cropX,cropY], Imgout1, col='#ff00ff')
    writeImage(Imgout2,sprintf("%s-ch2-B-nuc.png", fnControlImages),quality=98)
    Imgout2 = paintObjects(cytLI[cropX,cropY], Imgout1, col='#ff0000')
    writeImage(Imgout2,sprintf("%s-ch2-C-cyt.png", fnControlImages),quality=98)


    ## Img1 = readImage(fnch1)
    ## Img2 = readImage(fnch2)
    ## Img1c = mynormalize(imageData(Img1),minCh1, maxCh1)[cropX,cropY]
    ## Img2c = mynormalize(imageData(Img2),minCh2, maxCh2)[cropX,cropY]
    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img1c
    ## A[,,2] = Img1c
    ## A[,,3] = Img1c
    ## Imgout1 = Image(A,colormode="Color")
    ## writeImage(Imgout1,sprintf("%s-ch1-A-imgorg.png", fnControlImages),quality=98)
    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img2c
    ## A[,,2] = Img2c
    ## A[,,3] = Img2c
    ## Imgout1 = Image(A,colormode="Color")
    ## writeImage(Imgout1,sprintf("%s-ch2-A-imgorg.png", fnControlImages),quality=98)

    png(width=800,height=800,file=sprintf("%s-transformation1.png", fnControlImages))
    plotTransformationField(xt,yt,xn,yn,nn)
    dev.off()

    png(width=800,height=800,file=sprintf("%s-transformation2.png", fnControlImages))
    plotTransformationField(xt,yt,xn,yn,nn,maxDistZ = 5)
    dev.off()

    endProcessing(time3)
  }

  if (!is.null(fnGalleryImages)) {
    time3 = startProcessing("Save gallery images")
    x = round(F[,"m1.ch1.g.x"])
    y = round(F[,"m1.ch1.g.y"])
    I = which((x > 500) & (x < 1500) & (y > 500) & (y < 1500))
    if (length(I) > 20) {
      I = sample(I,20)
    }
    Fsample = F[I,]
    x=x[I]
    y=y[I]

    GImg = getGalleryImages(I, x, y, Img1, Img2, nucLI, cytLI, minCh1, maxCh1, minCh2, maxCh2)
    save(GImg, Fsample, x, y, I, file=fnGalleryImages)
    endProcessing(time3)
  }

  time2 = Sys.time()
  cat("Finished processing at ", format(time2, "%a %b %d %X %Y"), "\n")
  cat("Processing time: ",format(difftime(time2, time1, tz = "",units = c("auto"))), "\n")
  return(invisible(Fwell))
}

process.images.20x <- function(S,
                               displayImages = FALSE,
                               minCh1 = 0.00533,
                               maxCh1 = 0.04,
                               minCh2 = 0.0012,
                               maxCh2 = 0.015,
                               minCh3 = 0.0018,
                               maxCh3 = 0.010) {

  fnch1 = S$fnch1
  fnch2 = S$fnch2
  fnch3 = S$fnch3
  fnFlatFieldModel = S$fnFlatFieldModel
  fnControlImages = S$fnControlImages
  fnFeatures = S$fnFeatures
  fnFeaturesWell = S$fnFeaturesWell
  fnLabelImages = S$fnLabelImages
  fnGalleryImages = S$fnGalleryImages
  cropX = S$cropX
  cropY = S$cropY

  time1 = Sys.time()
  cat("Start image processing at ", format(time1, "%a %b %d %X %Y"), "\n")

  time3 = startProcessing("Read images")
  Img1 = readImage(fnch1)
  Img2 = readImage(fnch2)
  Img3 = readImage(fnch3)
  endProcessing(time3)

  ## time3 = startProcessing("Flat field correction")
  ## load(fnFlatFieldModel)
##   load("B2.rda")
##   ff2$B = B2
  ## Img1 = correctFlatField(Img1,ff1)
##   ff2$B = ff2$B / 2
  ## Img2 = correctFlatField(Img2,ff2)
  ## endProcessing(time3)

  ## s = 40
  ## M <- outer(dnorm(seq(-round(4*s),round(4*s),by=1),sd=s),dnorm(seq(-round(4*s),round(4*s),by=1),sd=s))
  ## M <- M / sum(M)
  ## Img1cor <- Img1 / filter2(Img1,M)
  ## Img2cor <- Img2 / filter2(Img2,M)
  ## Img3cor <- Img3 / filter2(Img3,M)

  time3 = startProcessing("Smooth images")
  Img1smooth = gblur(Img1,s=1)
  Img2smooth = gblur(Img2,s=1)
  Img3smooth = gblur(Img3,s=1)
  endProcessing(time3)

  time3 = startProcessing("Segmentation nucleus")
  nuct = thresh(Img1smooth, 10, 10, 0.0001)
  seed = bwlabel(opening(nuct,kern=makeBrush(5,shape="disc")))
  nuct2 = fillHull(thresh(Img1smooth, 20, 20, 0.00005))
  nucLI = propagate(Img1smooth, seed, mask=nuct2)
  rm(nuct); rm(seed); rm(nuct2)
  endProcessing(time3)

  time3 = startProcessing("Segmentation pH3")
  pH3nuct = thresh(Img3smooth, 10, 10, 0.0002)
  pH3seed = bwlabel(opening(pH3nuct,kern=makeBrush(5,shape="disc")))
  pH3nuct2 = fillHull(thresh(Img3smooth, 20, 20, 0.00005))
  pH3LI = propagate(Img3smooth, pH3seed, mask=pH3nuct2)
  rm(pH3nuct); rm(pH3seed); rm(pH3nuct2)
  endProcessing(time3)

  time3 = startProcessing("Segmentation tub")
  cytMask1 = opening(Img2smooth > 0.00405)
  cytMask2 = thresh(Img2smooth, 50,50, offset=0.000001)
  cytMask3 = cytMask1 + cytMask2 + nucLI
  cytMask3[cytMask3 > 1] = 1
  cytLI = propagate(Img2smooth, nucLI,lambda=1.0e-6, mask = cytMask3)
  rm(cytMask1); rm(cytMask2); rm(cytMask3)
  # cytLI = propagate(Img2smooth, nucLI,lambda=1.0, mask = filter2(nucLI,makeBrush(15,shape="disc")) > 0.00000001)

  ## SHMD = 0.008
  ## SHMD = 1.5*quantile(imageData(Img2),probs=0.05)
  ## minCh2 = quantile(imageData(Img2),probs=0.01)
  ## maxCh2 = 1.5 * quantile(imageData(Img2),probs=0.99)
  ## cytt = Img2smooth > SHMD
  ## cytLI = propagate(Img2, seeds = nucLI, mask = cytt, lambda=1e-20)
  endProcessing(time3)

  map <- tapply(imageData(nucLI)[pH3LI > 0],imageData(pH3LI)[pH3LI > 0], function(x) {
    x = x[x > 0]
    if (length(x) > 0) {
      lev <- unique(x)
      f <- factor(x, levels=lev)
      t <- table(f)
      res <- lev[which.max(t)]
    } else {
      res <- NA
    }
  } )

  time3 = startProcessing("Feature Extraction")

  ## ftnuc    = computeFeatures(nucLI,  list(DAPI=Img1cor, Tub=Img2cor, pH3=Img3cor), haralick.scales=c(1, 2, 4, 8), xname="nucleus")
  ## pftnuc  = computeFeatures(nucLI,  list(DAPI=Img1cor, Tub=Img2cor, pH3=Img3cor), haralick.scales=c(1, 2, 4, 8), xname="nucleus", properties=TRUE)
  ## ftcyt     = computeFeatures(cytLI,   list(DAPI=Img1cor, Tub=Img2cor, pH3=Img3cor), haralick.scales=c(1, 2, 4, 8), xname="cell")
  ## pftcyt   = computeFeatures(cytLI,   list(DAPI=Img1cor, Tub=Img2cor, pH3=Img3cor), haralick.scales=c(1, 2, 4, 8), xname="cell", properties=TRUE)
  ## ftpH3   = computeFeatures(pH3LI, list(DAPI=Img1cor, Tub=Img2cor, pH3=Img3cor), haralick.scales=c(1, 2, 4, 8), xname="nucleusPH3")
  ## pftpH3 = computeFeatures(pH3LI, list(DAPI=Img1cor, Tub=Img2cor, pH3=Img3cor), haralick.scales=c(1, 2, 4, 8), xname="nucleusPH3", properties=TRUE)

  ftnuc    = computeFeatures(nucLI,  list(DAPI=Img1), haralick.scales=c(1, 2), xname="nucleus")
  pftnuc  = computeFeatures(nucLI,  list(DAPI=Img1), haralick.scales=c(1, 2), xname="nucleus", properties=TRUE)
  if (is.null(dim(ftnuc))) {
    ftnuc = matrix(NA, nr=0, nc=nrow(pftnuc))
    colnames(ftnuc) = pftnuc$name
  }
  endProcessing(time3)
  time3 = startProcessing("Feature Extraction")
  ftcyt     = computeFeatures(cytLI,   list(DAPI=Img1, Tub=Img2), haralick.scales=c(1, 2), xname="cell")
  pftcyt   = computeFeatures(cytLI,   list(DAPI=Img1, Tub=Img2), haralick.scales=c(1, 2), xname="cell", properties=TRUE)
  if (is.null(dim(ftcyt))) {
    ftcyt = matrix(NA, nr=0, nc=nrow(pftcyt))
    colnames(ftcyt) = pftcyt$name
  }
  endProcessing(time3)
  time3 = startProcessing("Feature Extraction")
  ftpH3   = computeFeatures(pH3LI, list(pH3=Img3), haralick.scales=c(1, 2), xname="nucleusPH3")
  pftpH3 = computeFeatures(pH3LI, list(pH3=Img3), haralick.scales=c(1, 2), xname="nucleusPH3", properties=TRUE)
  if (is.null(dim(ftpH3))) {
    ftpH3 = matrix(NA, nr=0, nc=nrow(pftpH3))
    colnames(ftpH3) = pftpH3$name
  }
  endProcessing(time3)

  time3 = startProcessing("Feature Extraction 2")

  d <- sqrt((ftnuc[,"nucleus.DAPI.m.cx"] - 1024.5)^2 + (ftnuc[,"nucleus.DAPI.m.cy"] - 1024.5)^2)
  IpH3 <- which(d[map] <= 950)
  Inuc <- map[IpH3]
  F = cbind(ftnuc[Inuc,,drop=FALSE], ftcyt[Inuc,,drop=FALSE], ftpH3[IpH3,,drop=FALSE])
  Fwell1 = apply(F,2,mean, trim=0.1)
  Fwell2 = apply(F,2,sd)

  ## load("classifier.rda")
  ## pr = predict(model, F[,selectedFeatures])
  ## Fwell5 = table(pr$class)
  ## Fwell6 = apply(pr$posterior,2,mean)
  ## names(Fwell5) = sprintf("class.%s",c("unknown", "inter", "pro", "prometa", "meta", "anna", "telo", "dead", "oof", "earlyAnna"))
  ## names(Fwell6) = sprintf("posterior.%s",c("unknown", "inter", "pro", "prometa", "meta", "anna", "telo", "dead", "oof", "earlyAnna"))
  Inuc <- which((d <= 950) & !(1:length(d) %in% map))
  F = as.matrix(cbind(ftnuc[Inuc,,drop=FALSE], ftcyt[Inuc,,drop=FALSE]))
  Fwell3 = apply(F,2,mean, trim=0.1)
  Fwell4 = apply(F,2,sd)

  names(Fwell1) = sprintf("meanMitotic.%s",names(Fwell1))
  names(Fwell2) = sprintf("sdMitotic.%s",names(Fwell2))
  names(Fwell3) = sprintf("meanNonmitotic.%s",names(Fwell3))
  names(Fwell4) = sprintf("sdNonmitotic.%s",names(Fwell4))

  ## Fwell = c(count=length(Inuc)+length(IpH3), countpH3 = length(IpH3), ratiopH3 = length(IpH3) / (length(Inuc)+length(IpH3)), Fwell1, Fwell2, Fwell3, Fwell4, Fwell5, Fwell6)
  Fwell = c(count=length(Inuc)+length(IpH3), countpH3 = length(IpH3), ratiopH3 = length(IpH3) / (length(Inuc)+length(IpH3)), Fwell1, Fwell2, Fwell3, Fwell4)
  endProcessing(time3)

  ## d <- sqrt((ftnuc[,"nucleus.DAPI.m.cx"] - 1024.5)^2 + (ftnuc[,"nucleus.DAPI.m.cy"] - 1024.5)^2)
  ## IpH3 <- which(d[map] <= 950)
  ## if (length(IpH3) > 20) {
  ##   IpH3 <- sample(IpH3, 20)
  ## }
  ## Inuc <- map[IpH3]
  ## x <- ftnuc[Inuc,"nucleus.DAPI.m.cx"]
  ## y <- ftnuc[Inuc,"nucleus.DAPI.m.cy"]
  ## pos <- matrix(c(x,y), nr=length(x),nc=2)
  ## colnames(pos) <- c("x","y")

  ## Anno <- data.frame(imgNr = S$imgNr, IpH3 = IpH3, Inuc = Inuc, x=x, y=y)
  ## dir.create("outputm")
  ## save(Anno, file=sprintf("outputm/anno%d.rda", S$imgNr))
 
  ## FTnuc <- ftnuc[Inuc, ,drop=FALSE]
  ## FTcyt <- ftcyt[Inuc, , drop=FALSE]
  ## FTpH3 <- ftpH3[IpH3, , drop=FALSE]

  ## save(FTnuc, FTcyt, FTpH3, file=sprintf("outputm/F%d.rda", S$imgNr))

  ## A <- extractPatches(pos, list(DAPI=Img1, Tub = Img2, pH3 = Img3), normalize=TRUE, combine = TRUE, df=c(81,81))
  ## save(A, file=sprintf("outputm/A%d.rda", S$imgNr))

  ## Patch <- extractPatches(pos, list(DAPI=(Img1cor-0.5)/4, Tub = (Img2cor-0.5)/4, pH3 = (Img3cor - 0.5)/ 4), normalize=FALSE, combine = FALSE, df=c(81,81))
  ## save(Patch, file=sprintf("outputm/Patch%d.rda", S$imgNr))
  
  ## Fwell = c(count=nrow(F1),countpH3 = nrow(F3), ratiopH3 = nrow(F3) / nrow(F1), apply(Fnm,2,median))
  ## F1a = F1[map,,drop=FALSE]
  ## F2a = F2[map,,drop=FALSE]
  ## F3a = F3[,,drop=FALSE]
  ## colnames(F1a) = sprintf("MitoticNuc.%s",colnames(F1))
  ## colnames(F2a) = sprintf("MitoticTub.%s",colnames(F2))
  ## colnames(F3a) = sprintf("MitoticPH3.%s",colnames(F3))
  ## F = cbind(F1a, F2a, F3a)
  ## Fwell <- c(Fwell, apply(F,2,median))
  ## x = round(F[,"MitoticNuc.g.x"])
  ## y = round(F[,"MitoticNuc.g.y"])
  ## I = which((x > 41) & (x < 2006) & (y > 41) & (y < 2006))
  ## F <- F[I,,drop=FALSE]
  ## map <- map[I]
  ## x <- x[I]
  ## y <- y[I]
  ## endProcessing(time3)

  ## if (!is.null(fnFeatures)) {
  ##   time3 = startProcessing("Save features")
  ##   if (nrow(Fnm) > 20) {
  ##     I <- sample(1:nrow(Fnm),20)
  ##   } else {
  ##     I <- seq_len(nrow(Fnm))
  ##   }
  ##   Fnm <- Fnm[I,]
  ##   save(F,Fnm, file=fnFeatures)
  ##   endProcessing(time3)
  ## }
  if (!is.null(fnFeaturesWell)) {
    time3 = startProcessing("Save features per well")
    save(Fwell,file=fnFeaturesWell)
    endProcessing(time3)
  }
  ## if (!is.null(fnLabelImages)) {
  ##   time3 = startProcessing("Save label images")
  ##   save(nucLI, cytLI, file=fnLabelImages)
  ##   endProcessing(time3)
  ## }

  if (!is.null(fnControlImages)) {
    time3 = startProcessing("Save control images")

    Img1c = mynormalize(imageData(Img1),minCh1, maxCh1)[cropX,cropY]
    Img2c = mynormalize(imageData(Img2),minCh2, maxCh2)[cropX,cropY]
    Img3c = mynormalize(imageData(Img3),minCh3, maxCh3)[cropX,cropY]

    A = array(0,dim=c(dim(Img1c),3))
    A[,,1] = Img3c
    A[,,2] = Img2c
    A[,,3] = Img1c
    A[A > 1] = 1
    Imgout1 = Image(A,colormode="Color")
    writeImage(Imgout1,sprintf("%s-color-A-img.png", fnControlImages),quality=98)
    ## Imgout2 = paintObjects(nucLI[cropX,cropY], Imgout1, col='#ff00ff')
    ## writeImage(Imgout2,sprintf("%s-color-B-nuc.png", fnControlImages),quality=98)
    Imgout2 = paintObjects(cytLI[cropX,cropY], paintObjects(nucLI[cropX,cropY],Imgout1, col='#aa00aa'), col='#aa0000')
    writeImage(Imgout2,sprintf("%s-color-C-cyt.png", fnControlImages),quality=98)
    ## Imgout2 = paintObjects(cytLI[cropX,cropY], paintObjects(pH3LI[cropX,cropY],Imgout1, col='#ff00ff'), col='#ff0000')
    ## writeImage(Imgout2,sprintf("%s-color-D-pH3.png", fnControlImages),quality=98)

    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img1c
    ## A[,,2] = Img1c
    ## A[,,3] = Img1c
    ## A[A > 1] = 1
    ## Imgout1 = Image(A,colormode="Color")
    ## writeImage(Imgout1,sprintf("%s-ch1-A-img.png", fnControlImages),quality=98)
    ## Imgout2 = paintObjects(nucLI[cropX,cropY], Imgout1, col='#ff00ff')
    ## writeImage(Imgout2,sprintf("%s-ch1-B-nuc.png", fnControlImages),quality=98)
    ## ## Imgout2 = paintObjects(cytLI[cropX,cropY], paintObjects(pH3LI[cropX,cropY],Imgout1, col='#ff00ff'), col='#ff0000')
    ## ## writeImage(Imgout2,sprintf("%s-ch1-D-pH3.png", fnControlImages),quality=98)

    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img2c
    ## A[,,2] = Img2c
    ## A[,,3] = Img2c
    ## A[A > 1] = 1
    ## Imgout1 = Image(A,colormode="Color")
    ## writeImage(Imgout1,sprintf("%s-ch2-A-img.png", fnControlImages),quality=98)
    ## ## Imgout2 = paintObjects(nucLI[cropX,cropY], Imgout1, col='#ff00ff')
    ## ## writeImage(Imgout2,sprintf("%s-ch2-B-nuc.png", fnControlImages),quality=98)
    ## Imgout2 = paintObjects(cytLI[cropX,cropY], paintObjects(nucLI[cropX,cropY],Imgout1, col='#ff00ff'), col='#ff0000')
    ## writeImage(Imgout2,sprintf("%s-ch2-C-cyt.png", fnControlImages),quality=98)
    ## ## Imgout2 = paintObjects(cytLI[cropX,cropY], paintObjects(pH3LI[cropX,cropY],Imgout1, col='#ff00ff'), col='#ff0000')
    ## ## writeImage(Imgout2,sprintf("%s-ch2-D-pH3.png", fnControlImages),quality=98)

    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img3c
    ## A[,,2] = Img3c
    ## A[,,3] = Img3c
    ## A[A > 1] = 1
    ## Imgout1 = Image(A,colormode="Color")
    ## writeImage(Imgout1,sprintf("%s-ch3-A-img.png", fnControlImages),quality=98)
    ## ## Imgout2 = paintObjects(nucLI[cropX,cropY], Imgout1, col='#ff00ff')
    ## ## writeImage(Imgout2,sprintf("%s-ch3-B-nuc.png", fnControlImages),quality=98)
    ## ## Imgout2 = paintObjects(cytLI[cropX,cropY], paintObjects(nucLI[cropX,cropY],Imgout1, col='#ff00ff'), col='#ff0000')
    ## ## writeImage(Imgout2,sprintf("%s-ch3-C-cyt.png", fnControlImages),quality=98)
    ## Imgout2 = paintObjects(cytLI[cropX,cropY], paintObjects(pH3LI[cropX,cropY],Imgout1, col='#ff00ff'), col='#ff0000')
    ## writeImage(Imgout2,sprintf("%s-ch3-D-pH3.png", fnControlImages),quality=98)

    ## writeImage(permLabel(nucLI),sprintf("%s-li-A.png", fnControlImages),quality=100)
    ## writeImage(permLabel(cytLI),sprintf("%s-li-B.png", fnControlImages),quality=100)
    ## writeImage(permLabel(pH3LI),sprintf("%s-li-C.png", fnControlImages),quality=100)

    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img1c
    ## A[,,2] = Img1c
    ## A[,,3] = Img1c
    ## Imgout1 = Image(A,colormode="Color")
    ## writeImage(Imgout1,sprintf("%s-ch1-A-img.png", fnControlImages),quality=98)
    ## Imgout2 = paintObjects(nucLI[cropX,cropY], Imgout1, col='#ff00ff')
    ## writeImage(Imgout2,sprintf("%s-ch1-B-nuc.png", fnControlImages),quality=98)
    ## Imgout2 = paintObjects(cytLI[cropX,cropY], Imgout1, col='#ff0000')
    ## writeImage(Imgout2,sprintf("%s-ch1-C-cyt.png", fnControlImages),quality=98)

    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img2c
    ## A[,,2] = Img2c
    ## A[,,3] = Img2c
    ## Imgout1 = Image(A,colormode="Color")
    ## writeImage(Imgout1,sprintf("%s-ch2-A-img.png", fnControlImages),quality=98)
    ## Imgout2 = paintObjects(nucLI[cropX,cropY], Imgout1, col='#ff00ff')
    ## writeImage(Imgout2,sprintf("%s-ch2-B-nuc.png", fnControlImages),quality=98)
    ## Imgout2 = paintObjects(cytLI[cropX,cropY], Imgout1, col='#ff0000')
    ## writeImage(Imgout2,sprintf("%s-ch2-C-cyt.png", fnControlImages),quality=98)

    ## Img1 = readImage(fnch1)
    ## Img2 = readImage(fnch2)
    ## Img1c = mynormalize(imageData(Img1),minCh1, maxCh1)[cropX,cropY]
    ## Img2c = mynormalize(imageData(Img2),minCh2, maxCh2)[cropX,cropY]
    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img1c
    ## A[,,2] = Img1c
    ## A[,,3] = Img1c
    ## Imgout1 = Image(A,colormode="Color")
    ## writeImage(Imgout1,sprintf("%s-ch1-A-imgorg.png", fnControlImages),quality=98)
    ## A = array(0,dim=c(dim(Img1c),3))
    ## A[,,1] = Img2c
    ## A[,,2] = Img2c
    ## A[,,3] = Img2c
    ## Imgout1 = Image(A,colormode="Color")
    ## writeImage(Imgout1,sprintf("%s-ch2-A-imgorg.png", fnControlImages),quality=98)

    endProcessing(time3)
  }

  if (!is.null(fnGalleryImages)) {
    time3 = startProcessing("Save gallery images")

    ## I = which((x > 500) & (x < 1500) & (y > 500) & (y < 1500))
    ## if (length(I) > 20) {
    ##   I = sample(I,20)
    ## }
    ## Fsample = F[I,]
    ## x=x[I]
    ## y=y[I]
    if (length(map) > 0) {
      cat("write gallery images.\n")
      GImg = getGalleryImages(map, x, y, Img1, Img2, nucLI, cytLI, minCh1, maxCh1, minCh2, maxCh2)
    } else {
      cat("no gallery images written.\n")
      GImg <- NULL
    }
    save(GImg, x, y, file=fnGalleryImages)
    endProcessing(time3)
  }

  time2 = Sys.time()
  cat("Finished processing at ", format(time2, "%a %b %d %X %Y"), "\n")
  cat("Processing time: ",format(difftime(time2, time1, tz = "",units = c("auto"))), "\n")
  ## Fwell <- c(a=0.0,b=0.0,c=0.0)
  return(invisible(Fwell))
}

