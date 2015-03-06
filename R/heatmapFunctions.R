##-------------------------------------------------------------------------
## aigmoid
##-------------------------------------------------------------------------
sigmoid = function(x, theta=4, slope=1)
  sign(x) * (atan((abs(x)-theta)/slope)-atan(-theta/slope))/pi


##-------------------------------------------------------------------------
## toRaster
## x a numeric matrix, returns a matrix of same shape with colour strings
##-------------------------------------------------------------------------
toRaster = function(x, cuts, col) {
  cux =cut(x,cuts,include.lowest = TRUE, labels=FALSE)
  rv = x
  rv[] = col[cux]
  return(rv)
}

##--------------------------------------------------
## transform matrix x into a 3D array format
##--------------------------------------------------
to3Darray = function(x) {
    spl = strsplit(colnames(x), split="_")
    queryGenes  = sapply(spl, "[", 1)
    queryFeats  = sapply(spl, "[", 2)
    uQueryGenes = unique(queryGenes)
    uQueryFeats = unique(queryFeats)

    rv = array(NA_real_, dim=c(nrow(x), length(uQueryGenes), length(uQueryFeats)),
               dimnames = list(template  = rownames(x),
               queryGene = uQueryGenes,
               queryFeat = uQueryFeats))
    for(i2 in dimnames(rv)[[2]])
        for(i3 in dimnames(rv)[[3]]) {
            cn = paste(i2, i3, sep="_")
            if(cn %in% colnames(x))
                rv[, i2, i3] = x[, cn]
        }
    return(rv)
}

##--------------------------------------------------
## and the reverse... (but looses dimnames)
##--------------------------------------------------
toMatrix = function(x) {
    dim(x) = c(dim(x)[1], prod(dim(x)[-1]))
    return(x)
}

##--------------------------------------------------------------
## determine 'hclust ordering' of one dimension i of array x
##--------------------------------------------------------------
orderDim = function(x, i) {
   px = toMatrix(aperm(x, c(i, setdiff(seq(along=dim(x)), i))))
   theCorr = cor(t(px), use = "pairwise.complete.obs")
   hc = hclust(as.dist(trsf(theCorr)))
   return(hc$order)
}

##--------------------------------------------------------------
## transform a correlation coefficient (in [-1,1]) to a distance
##--------------------------------------------------------------
trsf = function(x) {
    exp(-x)-exp(-1)
}


makeDall = function(D) {
  rg = seq_len(min(10, dim(D)[3]))
  stopifnot(identical(dimnames(D)[[3]][rg],
            c("4x.count", "4x.ratioMitotic", "10x.meanNonmitotic.cell.0.s.area", "4x.countpH3", "4x.areaNucAll",
              "4x.areapH3All", "4x.areaNucH8", "4x.areaNucH3", "4x.intNucH4", "4x.areaNucH10")[rg]))

  dimnames(D)[[3]][rg] = c("cell number", "mitotic ratio", "area (non-|mitotic cells)", "number of pH3|positive cells", "nucleus area",
                       "pH3 area", "nucleus area (H8)", "nucleus area (H3)", "nucleus area (H4)", "nucleus area (H10)")[rg]

  for (i in seq_len(dim(D)[3]))
    D[,,i] = D[,,i] / mad(D[,,i])

  return(D)
}

##----------------------------------------------------------
## This function plots dim(x)[3] heatmaps of x[,,i].
## x is a 3D array (dimensions already sorted for plotting)
##----------------------------------------------------------
myHeatmap = function (x, cuts, col,
  fontsize = 18, colnames = TRUE, rownames = FALSE) {
    stopifnot(is.array(x), length(dim(x)) == 3)
    use = 0.99

    rx = toRaster(x, cuts = cuts, col = col)
    pushViewport(viewport(xscale = c(if (rownames) -1 else 0.5, dim(x)[3] + 0.5),
                          yscale = c(0, if (colnames) 1/0.8 else 1),
                          y = unit(use/2, "npc"), height = unit(use, "npc"),
                          x = unit(use/2, "npc"), width = unit(use, "npc")))
    colnm = dimnames(x)[[3]]
    if (colnames) {
        stopifnot(!is.null(colnm))
        colnm = strsplit(colnm, split = "|", fixed = TRUE)
    }
    rownm = dimnames(x)[[1]]
    if (rownames) 
        stopifnot(!is.null(rownm))
    
    for (i3 in seq_len(dim(x)[3])) {
        grid.raster(rx[, , i3], x = unit(i3, "native"), width = unit(1, "native"),
                    y = unit(0.5, "native"), height = unit(1, "native"),
                    interpolate = FALSE)
        if (colnames) {
            nLines = length(colnm[[i3]])
            for (k in seq_len(nLines))
              grid.text(paste(" ", colnm[[i3]][k], sep = ""),
                        x = unit(i3, "native") + unit((if (nLines > 1) { (k - 1)/(nLines - 1) } else 0) - 0.5, "char"),
                        y = unit(1, "native"),
                        hjust = 0,
                        vjust = 1, rot = 90,
                        gp = gpar(fontsize = fontsize))
        }
        if (rownames) {
            grid.text(paste(rownm, " ", sep = ""),
                      x = unit(0.5, "native"),
                      y = unit(seq(along = rownm)/length(rownm),"native"),
                      hjust=1,
                      vjust=1,
                      gp = gpar(fontsize = fontsize))
        }
    }
    grid.polyline(x = unit(rep(seq_len(dim(x)[3]-1)+0.5,each=2), "native"),
                  y = unit(rep(c(1.0,1.5),times=dim(x)[3]-1), "native"),
                  id = rep(seq_len(dim(x)[3]-1),each=2), gp=gpar(col="gray50") )
    popViewport()
}

##------------------------------------------------------------------------------------
## Plot a schematic of the fit matrices
##------------------------------------------------------------------------------------
fitMatrices = function(D, Fit, R, main, filename, cuts, col, fac=c(1, 1)){
  stopifnot(identical(dim(D), dim(Fit)),
            identical(dim(D), dim(R)))

  ## Layout
  nGenes = dim(D)[1]
  nCols  = prod(dim(D)[2:3])
  xsep   = 20
  whm    = nCols / (3*nCols + 2*xsep)
  wsep   = xsep  / (3*nCols + 2*xsep)
  myVP = function(n) viewport(x=whm*(n-0.5)+wsep*(n-1), width=whm, y=0.5, height=1, default.units="npc")

  png(filename, width = (3*nCols+2*xsep)*fac[1], height = nGenes*fac[2])

  pushViewport(myVP(1))
  myHeatmap(D, cuts=cuts, col=col)
  popViewport()
  pushViewport(myVP(2))
  myHeatmap(Fit, cuts=cuts, col=col)
  popViewport()
  pushViewport(myVP(3))
  myHeatmap(R, cuts=cuts, col=col)
  popViewport()

  grid.text(main, x=2*whm+1.5*wsep, y=  1, just = c(0.5, 1.5), gp=gpar(fontsize=36))
  grid.text("=",  x=1*whm+0.5*wsep, y=0.4, just = c(0.5, 0.5), gp=gpar(fontsize=36))
  grid.text("+",  x=2*whm+1.5*wsep, y=0.4, just = c(0.5, 0.5), gp=gpar(fontsize=36))

  dev.off()
  invisible(NULL)
}
##------------------------------------------------------------------------------------
## Call myHeatmap and save as PDF or PNG
## Currently I prefer PNG since the colours look crisper and no interpolation is done
##------------------------------------------------------------------------------------
saveHeatmapFile = function(x, ppiw = 0.5, ppih = 1, fac = 1,
                           filename, cuts,
                           fontsize = 1,
                           device = "png",
                           ...) {
  d = dim(x)
  for (dev in device) {
    fn = if(!missing(filename)) {
      sprintf("%s.%s",filename, dev)
    } else {
      sprintf("heatmap-%s.%s", sub("^v", "", deparse(substitute(x))), dev)
    }
    cat("Writing", fn, "\n")
    if (dev == "png") {
      png(fn, width = prod(d[-1])/ppiw*fac, height = d[1]/ppih*fac)
      myHeatmap(x, cuts=cuts, fontsize=32*fontsize, ...)
      dev.off()
    }
    if (dev == "pdf") {
      cairo_pdf(fn, width = prod(d[-1])/ppiw*fac, height = d[1]/ppih*fac)
      myHeatmap(x, cuts=cuts, fontsize=18*fontsize, ...)
      dev.off()
    }
    if (!(dev %in% c("png", "pdf"))) {
      stop("device", dev, "not supported")
    }
  }
}


mySymlink = function(from, to) {
  stopifnot(length(from)==1)
  toRemove = to[file.exists(to)]
  if(length(toRemove)>0)
    file.remove(toRemove)
  res = file.symlink(from, to)
  stopifnot(all(res))
}
