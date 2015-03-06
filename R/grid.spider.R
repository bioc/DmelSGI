orderSpiderAxis <- function(X) {
  TSP = TSP(dist(t(X)))
  TSPobj = solve_TSP(TSP,method="repetitive_nn")
  #image(TSP,T)
  P = as.integer(TSPobj)
  P
}

grid.spider <- function(v, col, col.arms="black", dlim=NULL) {
  if(is.null(nrow(v))) {
    v = matrix(v,nrow=1)
  }
  if (is.null(dlim)) {
    dlim = range(v)
    if (dlim[1] > 0.0) {
      dlim[1] = 0.0
    }
    ref = 0.0
  }
  #vp = viewport(xscale=c(0.8,ncol(X)+0.2),yscale=c(0,1))
  vp = viewport(layout=grid.layout(nrow=1,ncol=1,widths=1,heights=1,respect=TRUE))
  pushViewport(vp)
  vp = viewport(layout.pos.row=1,layout.pos.col=1,
                xscale=c(-1.08,1.08),yscale=c(-1.08,1.08))
  pushViewport(vp)
  n = ncol(v)
  v = cbind(v,v[,1])
  M = matrix(rep(1:ncol(v),each=nrow(v)))
  #  ID = matrix(rep(1:nrow(v),times=ncol(v)))
  x=as.vector(M)
  y=as.vector(v)
  x = c(x,rev(x))
  y = c(y,rep(ref,length(y)))
  x2 = (y - dlim[1])/diff(dlim) * sin(2*pi*x/n)
  y2 = (y - dlim[1])/diff(dlim) * cos(2*pi*x/n)
  x3 = 1.05*sin(2*pi*(x[1:(length(v)-1)]-0.5)/n)
  y3 = 1.05*cos(2*pi*(x[1:(length(v)-1)]-0.5)/n)
  x4 = 1.05*sin(2*pi*(x[1:(length(v)-1)]+0.5)/n)
  y4 = 1.05*cos(2*pi*(x[1:(length(v)-1)]+0.5)/n)
  #   grid.polyline(x=x2,y=y2,id=as.vector(ID),
  #                default.units="native",gp=gpar(col=col))
  for (i in 1:n) {
    grid.polygon(c(x3[i],x4[i],0.0),c(y3[i],y4[i],0.0),
                 default.units="native",gp=gpar(fill=col.arms[i],col=NA))
  }
  grid.polygon(x=x2,y=y2,
               default.units="native",gp=gpar(fill=col))
  x2 = sin(2*pi*1:n/n)
  y2 = cos(2*pi*1:n/n)
  grid.polyline(c(x2,rep(0.0,n)),c(y2,rep(0.0,n)),id=c(1:n,1:n),
                default.units="native")
#   grid.polyline(c(x2,rep(0.0,n)),c(y2,rep(0.0,n)),id=c(1:n,1:n),
#                 default.units="native",gp=gpar(col=col.arms))
#   grid.polyline(c(x3,x4),c(y3,y4),id=c(1:n,1:n),
#                 default.units="native",gp=gpar(col=col.arms,lwd=4))
  popViewport()
  popViewport()
  invisible(NULL)
}

grid.spider.legend <- function(vn, col.arms="black", dlim=NULL) {
  v = rep(1,length(vn))
  if(is.null(nrow(v))) {
    v = matrix(v,nrow=1)
  }
  if (is.null(dlim)) {
    dlim = c(0,1)
    ref = 0.0
  }
  #vp = viewport(xscale=c(0.8,ncol(X)+0.2),yscale=c(0,1))
  vp = viewport(layout=grid.layout(nrow=1,ncol=1,widths=1,heights=1,respect=TRUE))
  pushViewport(vp)
  vp = viewport(layout.pos.row=1,layout.pos.col=1,
                xscale=c(-5.03,5.03),yscale=c(-5.03,5.03))
  pushViewport(vp)
  n = ncol(v)
  v = cbind(v,v[,1])
  M = matrix(rep(1:ncol(v),each=nrow(v)))
  #  ID = matrix(rep(1:nrow(v),times=ncol(v)))
  x=as.vector(M)
  y=as.vector(v[,c(seq_along(vn),1)])
  x = c(x,rev(x))
  y = c(y,rep(ref,length(y)))
#   x2 = (y - dlim[1])/diff(dlim) * sin(2*pi*x/n)
#   y2 = (y - dlim[1])/diff(dlim) * cos(2*pi*x/n)
  #   grid.polyline(x=x2,y=y2,id=as.vector(ID),
  #                default.units="native",gp=gpar(col=col))
  #  grid.polygon(x=x2,y=y2,
  #               default.units="native",gp=gpar(fill=col))
  x2 = sin(2*pi*1:n/n)
  y2 = cos(2*pi*1:n/n)
  x3 = 1.05*sin(2*pi*(x[1:(length(v)-1)]-0.5)/n)
  y3 = 1.05*cos(2*pi*(x[1:(length(v)-1)]-0.5)/n)
  x4 = 1.05*sin(2*pi*(x[1:(length(v)-1)]+0.5)/n)
  y4 = 1.05*cos(2*pi*(x[1:(length(v)-1)]+0.5)/n)
  x5 = 1.1*sin(2*pi*1:n/n)
  y5 = 1.1*cos(2*pi*1:n/n)
  grid.polyline(c(x2,rep(0.0,n)),c(y2,rep(0.0,n)),id=c(1:n,1:n),
                default.units="native",gp=gpar(col=col.arms))
  grid.polyline(c(x3,x4),c(y3,y4),id=c(1:n,1:n),
                default.units="native",gp=gpar(col=col.arms,lwd=4))
  #  grid.text(label=vn, x2, y2,rot=seq_along(vn)/(2*pi))
  grid.text(label=vn, x5[seq_along(vn)], y5[seq_along(vn)],just=c("left","center"),
            default.units="native",rot=90+(360)*(1-seq_along(vn)/length(vn)),
            gp=gpar(col=col.arms,lwd=2))
  popViewport()
  popViewport()
  invisible(NULL)
}

