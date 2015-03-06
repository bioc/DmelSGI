library(EBImage)
library(MASS)
library(grid)
library(RColorBrewer)
library(genefilter)
library(hwriter)
library(e1071)

mynormalize <- function(Img, minInt, maxInt) {
  Img[Img > maxInt] = maxInt
  Img[] = Img[] - minInt
  Img[Img < 0] = 0
  Img[] = Img[] / (maxInt - minInt)
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

local.maxima <- function(Img, Mask, mindist = 3.5, min.intensity = 0.0, verbose = FALSE) {
  if (verbose) {
    time1 = Sys.time()
    cat("Local maximum search ...")
  }
  N <- dim(Img)[1]
  M <- dim(Img)[2]

  LM1 = matrix(0,nr=N-2,nc=M-2)

  Img[!Mask[]] = 0
  Img1 = Img[2:(N-1),2:(M-1)]

  LM1[(Mask[2:(N-1),2:(M-1)] > 0) & (Img1 >= Img[1:(N-2),1:(M-2)]) & (Img1 >= Img[1:(N-2),2:(M-1)]) & (Img1 >= Img[1:(N-2),3:M])  & (Img1 >= Img[2:(N-1),1:(M-2)]) & (Img1 >= Img[2:(N-1),3:M]) & (Img1 >= Img[3:N,    1:(M-2)]) & (Img1 >= Img[3:N,    2:(M-1)]) & (Img1 >= Img[3:N,    3:M]) ] = 1
  LM = matrix(0, nr=N, nc=M)
  LM[2:(N-1),2:(M-1)] = LM1

  X = matrix(rep(1:N,times=M),nr=N,nc=M)
  Y = matrix(rep(1:M,each=N),nr=N,nc=M)
  I = which(LM > 0.5)
  if (verbose) {
    cat("\nNumber of local maxima: ",length(I)," ... ")
  }
  I = I[(X[I] > mindist) & (X[I] < N - mindist) & (Y[I] > mindist) & (Y[I] < M - mindist)]
  I = I[order(-Img[I])]
  x = X[I]
  y = Y[I]

  disk = matrix(FALSE,nr=2*mindist+1,nc=2*mindist+1)
  X2 = matrix(rep(-mindist:mindist,times=2*mindist+1),nr=2*mindist+1,nc=2*mindist+1)
  Y2 = matrix(rep(-mindist:mindist,each=2*mindist+1),nr=2*mindist+1,nc=2*mindist+1)

  disk[X2^2 + Y2^2 <= mindist^2] = TRUE
  is.neighbor = matrix(0,nr=N,nc=M)

  is.selected = rep(FALSE,length(I))
  if (length(I) > 0) {
    for (i in 1:length(I)) {
      if (is.neighbor[x[i],y[i]] < 0.5) {
        is.selected[i] = TRUE
        J = rep((x[i]-mindist):(x[i]+mindist),times=2*mindist+1) + N*rep((y[i]-mindist):(y[i]+mindist),each=2*mindist+1)
        J = J[disk]
        is.neighbor[J] = 1
      }
    }
  }

  if (verbose) {
    cat("\nnumber of local maxima after removing maxima with distance smaller ",mindist,": ", sum(is.selected)," ... ")
  }
  LM[] = 0
  LM[I[is.selected]] = 1

  LM[Img[] < min.intensity] = 0
  
  if (verbose) {
    time2 = Sys.time()
    cat("finished. Time elapsed: ",format(difftime(time2, time1, tz = "",units = c("auto"))),".\n")
  }
  return(LM)
}

label.mask <- function(LM, Img, Mask, verbose = FALSE) {
  if (verbose) {
    time1 = Sys.time()
    cat("Label image with mask ... ")
  }
  LabelImage1 = bwlabel(LM)
  LabelImage2 = propagate(Img, LabelImage1, Mask, lambda=100)
  if (verbose) {
    time2 = Sys.time()
    cat("finished. Time elapsed: ",format(difftime(time2, time1, tz = "",units = c("auto"))),".\n")
  }
  return(LabelImage2)
}

getGalleryImages <- function(I, x, y, Img1, Img2, nucLI, cytLI, minCh1, maxCh1, minCh2, maxCh2) {
  GImg = array(0.0, dim=c(length(I), 9, 81, 81, 3))

  for (i in 1:length(I)) {
    cropX = (x-40):(x+40)
    cropY = (y-40):(y+40)
    nucLI2 = nucLI
    nucLI2[nucLI2 != I[i]] = 0
    cytLI2 = cytLI
    cytLI2[cytLI2 != I[i]] = 0

    Img1c = mynormalize(imageData(Img1)[cropX,cropY],minCh1, maxCh1)
    Img2c = mynormalize(imageData(Img2)[cropX,cropY],minCh2, maxCh2)

    A = array(0,dim=c(dim(Img1c),3))
    A[,,3] = Img1c
    A[,,2] = Img2c
    Imgout1 = Image(A,colormode="Color")
    GImg[i,1,,,] = imageData(Imgout1)
    Imgout2 = paintObjects(nucLI2[cropX,cropY], Imgout1, col='#ff00ff')
    GImg[i,2,,,] = imageData(Imgout2)
    Imgout2 = paintObjects(cytLI2[cropX,cropY], Imgout1, col='#ff0000')
    GImg[i,3,,,] = imageData(Imgout2)

    A = array(0,dim=c(dim(Img1c),3))
    A[,,1] = Img1c
    A[,,2] = Img1c
    A[,,3] = Img1c
    Imgout1 = Image(A,colormode="Color")
    GImg[i,4,,,] = imageData(Imgout1)
    Imgout2 = paintObjects(nucLI2[cropX,cropY], Imgout1, col='#ff00ff')
    GImg[i,5,,,] = imageData(Imgout2)
    Imgout2 = paintObjects(cytLI2[cropX,cropY], Imgout1, col='#ff0000')
    GImg[i,6,,,] = imageData(Imgout2)

    A = array(0,dim=c(dim(Img1c),3))
    A[,,1] = Img2c
    A[,,2] = Img2c
    A[,,3] = Img2c
    Imgout1 = Image(A,colormode="Color")
    GImg[i,7,,,] = imageData(Imgout1)
    Imgout2 = paintObjects(nucLI2[cropX,cropY], Imgout1, col='#ff00ff')
    GImg[i,8,,,] = imageData(Imgout2)
    Imgout2 = paintObjects(cytLI2[cropX,cropY], Imgout1, col='#ff0000')
    GImg[i,9,,,] = imageData(Imgout2)
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
  nuct = thresh(Img1smooth, 4, 4, 0.0003)
  LM = local.maxima(Img1smooth,nuct,verbose=TRUE)
  nucLI = label.mask(LM, Img1smooth, nuct, verbose = TRUE)
  LabelImage1 = bwlabel(nuct)
  cytt = thresh(Img2smooth, 4, 4, 0.0003)
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
  Inuc = dnuc <= 900
  IH3p = (dH3p <= 900) & (nn$disttrafoZ <= 5)
  LCD = localDensity(Fnuc[,"x"], Fnuc[,"y"])
  Q1 = quantile(Fnuc[,"meanint"],probs=c(0.03,0.1,0.25,0.5,0.75,0.9,0.97),na.rm=TRUE)
  names(Q1) = sprintf("intNucQ%0.2f",probs=c(0.03,0.1,0.25,0.5,0.75,0.9,0.97))
  Q2 = quantile(FH3p[,"meanint"],probs=c(0.03,0.1,0.25,0.5,0.75,0.9,0.97),na.rm=TRUE)
  names(Q2) = sprintf("intH3pQ%0.2f",probs=c(0.03,0.1,0.25,0.5,0.75,0.9,0.97))
  Q3 = quantile(Fnuc[,"area"],probs=c(0.03,0.1,0.25,0.5,0.75,0.9,0.97),na.rm=TRUE)
  names(Q3) = sprintf("areaNucQ%0.2f",probs=c(0.03,0.1,0.25,0.5,0.75,0.9,0.97))
      # q=quantile(F[,1],na.rm=TRUE, probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95))
  H1 = hist(Fnuc[,"meanint"],breaks=c(-1,0.003165389,0.003391965,0.003584755,0.003775190,0.003983502,0.004202335,0.004467842,0.004814221,0.005061701,0.005405753,0.005954504,100), plot=FALSE)$counts
  if (sum(H1) > 0) { H1 = H1 / sum(H1) }
  names(H1) = sprintf("intNucH%d",1:12)
      # q=quantile(F[,3],na.rm=TRUE, probs=c(0.1,0.25,0.5,0.85,0.95))
  H2 = hist(FH3p[,"meanint"],breaks=c(-1,0.002631055,0.002771038,0.003014747,0.003828319,0.004562448,100), plot=FALSE)$counts
  if (sum(H2) > 0) { H2 = H2 / sum(H2) }
  names(H2) = sprintf("intH3pH%d",1:6)
      # q=quantile(F[,2],na.rm=TRUE, probs=c(0.1,0.25,0.5,0.85,0.95))
  H3 = hist(Fnuc[,"area"],breaks=c(-1,8,10,13,14,16,17,19,21,23,27,33,10000000000000), plot=FALSE)$counts
  if (sum(H3) > 0) { H3 = H3 / sum(H3) }
  names(H3) = sprintf("areaNucH%d",1:12)
  Fwell = c(count =  sum(Inuc,na.rm=TRUE),
    countpH3      =  sum(IH3p,na.rm=TRUE),
    isMitotic     =  sum(Fnuc[Inuc,"isMitotic"],na.rm=TRUE),
    ratioMitotic  =  sum(Fnuc[Inuc,"isMitotic"],na.rm=TRUE) / sum(Inuc,na.rm=TRUE),
    areaNuc       = mean(Fnuc[Inuc,"area"],na.rm=TRUE),
    areaNucSD     =   sd(Fnuc[Inuc,"area"],na.rm=TRUE),
    areapH3       = mean(FH3p[IH3p,"area"],na.rm=TRUE),
    areapH3SD     =   sd(FH3p[IH3p,"area"],na.rm=TRUE),
    intNuc        = mean(Fnuc[Inuc,"meanint"],na.rm=TRUE),
    intNucSD      =   sd(Fnuc[Inuc,"meanint"],na.rm=TRUE),
    intpH3        = mean(FH3p[IH3p,"meanint"],na.rm=TRUE),
    intpH3SD      =   sd(FH3p[IH3p,"meanint"],na.rm=TRUE),
    countAll         =  sum(Inuc,na.rm=TRUE),
    countpH3All      =  sum(IH3p,na.rm=TRUE),
    isMitoticAll     =  sum(Fnuc[,"isMitotic"],na.rm=TRUE),
    ratioMitoticAll  =  sum(Fnuc[,"isMitotic"],na.rm=TRUE) / sum(Inuc,na.rm=TRUE),
    areaNucAll       = mean(Fnuc[,"area"],na.rm=TRUE),
    areaNucSDAll     =   sd(Fnuc[,"area"],na.rm=TRUE),
    areapH3All       = mean(FH3p[,"area"],na.rm=TRUE),
    areapH3SDAll     =   sd(FH3p[,"area"],na.rm=TRUE),
    intNucAll        = mean(Fnuc[,"meanint"],na.rm=TRUE),
    intNucSDAll      =   sd(Fnuc[,"meanint"],na.rm=TRUE),
    intpH3All        = mean(FH3p[,"meanint"],na.rm=TRUE),
    intpH3SDAll      =   sd(FH3p[,"meanint"],na.rm=TRUE),
    LCD$mean, Q1, Q2, Q3, H1, H2, H3
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
    writeImage(permLabel(LabelImage1)[cropX,cropY], sprintf("%s-li-C-seg.png", fnControlImages),quality=100)
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

    ## png(width=800,height=800,file=sprintf("%s-transformation1.png", fnControlImages))
    ## plotTransformationField(xt,yt,xn,yn,nn)
    ## dev.off()

    ## png(width=800,height=800,file=sprintf("%s-transformation2.png", fnControlImages))
    ## plotTransformationField(xt,yt,xn,yn,nn,maxDistZ = 5)
    ## dev.off()

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

process.images.20x <- function(fnch1, fnch2,
                               displayImages = FALSE,
                               fnFlatFieldModel = "/g/huber/users/befische/experiments/MichaelBoutros/ChromatinSet/flatFieldCorrection/flatFieldModel.rda",
                               fnControlImages = NULL,
                               fnFeatures = NULL,
                               fnFeaturesWell = NULL,
                               fnLabelImages = NULL,
                               fnGalleryImages = NULL,
                               cropX = 940:1059,
                               cropY = 950:1049,
                               minCh1 = 0.00533,
                               maxCh1 = 0.025,
                               minCh2 = 0.008,
                               maxCh2 = 0.030 ) {
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
  Img1smooth = gblur(Img1,s=1)
  Img2smooth = gblur(Img2,s=3)
  endProcessing(time3)

  time3 = startProcessing("Segmentation")
  nuct = thresh(Img1smooth, 10, 10, 0.001)
  nucLI = bwlabel(nuct)

  SHMD = 0.008
  SHMD = 1.5*quantile(imageData(Img2),probs=0.05)
  minCh2 = quantile(imageData(Img2),probs=0.01)
  maxCh2 = 1.5 * quantile(imageData(Img2),probs=0.99)
  cytt = Img2smooth > SHMD
  cytLI = propagate(Img2, seeds = nucLI, mask = cytt, lambda=1e-20)
  endProcessing(time3)

  time3 = startProcessing("Feature Extraction")
  Fmask1ch1 = getFeatures(nucLI,Img1)[[1]]
  Fmask2ch2 = getFeatures(cytLI,Img2)[[1]]
  colnames(Fmask1ch1) = sprintf("m1.ch1.%s",colnames(Fmask1ch1))
  colnames(Fmask2ch2) = sprintf("m2.ch2.%s",colnames(Fmask2ch2))
  F = cbind(Fmask1ch1, Fmask2ch2)
  d = sqrt((F[,"m1.ch1.g.x"] - 1024)^2 + (F[,"m1.ch1.g.y"] - 1024)^2)
  Fwell = c(count=sum(d <= 900),apply(F[d <= 900,],2,median))
  endProcessing(time3)

  if (!is.null(fnFeatures)) {
    time3 = startProcessing("Save features")
    save(F,file=fnFeatures)
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


    Img1 = readImage(fnch1)
    Img2 = readImage(fnch2)
    Img1c = mynormalize(imageData(Img1),minCh1, maxCh1)[cropX,cropY]
    Img2c = mynormalize(imageData(Img2),minCh2, maxCh2)[cropX,cropY]
    A = array(0,dim=c(dim(Img1c),3))
    A[,,1] = Img1c
    A[,,2] = Img1c
    A[,,3] = Img1c
    Imgout1 = Image(A,colormode="Color")
    writeImage(Imgout1,sprintf("%s-ch1-A-imgorg.png", fnControlImages),quality=98)
    A = array(0,dim=c(dim(Img1c),3))
    A[,,1] = Img2c
    A[,,2] = Img2c
    A[,,3] = Img2c
    Imgout1 = Image(A,colormode="Color")
    writeImage(Imgout1,sprintf("%s-ch2-A-imgorg.png", fnControlImages),quality=98)


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

