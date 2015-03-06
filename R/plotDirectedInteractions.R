
plotPIdata <- function(X,gt,gq,show="summary",...) {
  p=par(mfrow=c(1,2))
  if (show == "summary") {
    xt = apply(X[,"xt",gt,,gq,],1,mean,na.rm=TRUE)
    xq = apply(X[,"xq",gt,,gq,],1,mean,na.rm=TRUE)
    pi = apply(X[,"pi",gt,,gq,],1,mean,na.rm=TRUE)
  } else {
    xt = X[,"xt",gt,show[1],gq,show[2]]
    xq = X[,"xq",gt,show[1],gq,show[2]]
    pi = X[,"pi",gt,show[1],gq,show[2]]
  }
  plot(xt,pi,main=sprintf("%s - %s",gt,gq),
       xlab=sprintf("%s main effect",gt),ylab="pi-score",...)
  abline(v=0,h=0)
  plot(xq,pi,main=sprintf("%s - %s",gt,gq),
       xlab=sprintf("%s main effect",gq),ylab="pi-score",...)
  abline(v=0,h=0)
  par(p)
  invisible(NULL)
}


plot2Phenotypes <- function(X,gt,gq,f1,f2,length=1,...) {
  X = X[c(f1,f2),,gt,,gq,]
  dn = dimnames(X)
  dim(X) = c(dim(X)[1:2],prod(dim(X)[3:4]))
  dimnames(X) = list(DmelSGI::hrNames(dn[[1]]),dn[[2]],NULL)
  X = aperm(X,c(3,1,2))
  X = abind(X, X[,,"xt"]+X[,,"xq"], X[,,"xt"]+X[,,"xq"]+X[,,"pi"],along=3)
  dimnames(X)[[3]][4] = "exp"
  dimnames(X)[[3]][5] = "data"
  r = apply(X[,,c("xt","xq","exp","data")],2,function(x) {
    r  = range(x,finite=TRUE)
    if (prod(sign(r)) > 0) {
      if (r[1] > 0) {
        r[1] = 0
      } else {
        r[2] = 0
      }
    }
    r
  } )
  d = apply(r, 2, diff)
  a = which.min(d)
  b = c(1,2)[-a]
  d = (d[b]-d[a])/2
  r[,a] = r[,a] + c(-d,d)
  layout(matrix(c(1,2),nrow=1,ncol=2),widths=c(1.0,0.25),heights=c(1.0,0.3),respect=TRUE)
  p = par(mar=c(4,4,4,4))
  plot(NULL,xlim=r[,1],ylim=r[,2],col="orange",pch=19,xlab=dimnames(X)[[2]][1],
       ylab=dimnames(X)[[2]][2],...)
  points(X[,,"xt"],col="#F7942D",pch=19,...)
  arrows(0.0,0.0,X[,1,"xt"],X[,2,"xt"],col="#F7942D",length=length,...)
  points(X[,,"xq"],col="#418AC9",pch=19,...)
  arrows(0.0,0.0,X[,1,"xq"],X[,2,"xq"],col="#418AC9",length=length,...)
  points(X[,,"exp"],col="#939598",pch=19,...)
  segments(X[,1,"xq"],X[,2,"xq"],X[,1,"exp"],X[,2,"exp"],col="#939598",...)
  segments(X[,1,"xt"],X[,2,"xt"],X[,1,"exp"],X[,2,"exp"],col="#939598",...)
  points(X[,,"data"],col="black",pch=19,...)
  arrows(X[,1,"exp"],X[,2,"exp"],X[,1,"data"],X[,2,"data"],col="black",length=length,...)
  abline(v=0,h=0)
  par(mar=c(0,0,0,0))
  plot(NULL,xlim=c(0,20),ylim=c(-4,15),bty="n",xaxt="n",yaxt="n")
  points(c(1,1,1,1),c(1,2,3,4),pch=19,col=c("black","#939598","#418AC9","#F7942D"),...)
  text(2,1,substitute(paste(gt,"+",gq),list(gt=gt,gq=gq)),adj=c(0,0.5),...)
  text(2,2,"expected",adj=c(0,0.5),...)
  text(2,3,substitute(gq,list(gq=gq)),adj=c(0,0.5),...)
  text(2,4,substitute(gt,list(gt=gt)),adj=c(0,0.5),...)
  par(p)
  invisible(NULL)
}
