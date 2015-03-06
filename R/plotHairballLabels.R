
plotHairballLabels <- function(g, co, Labels, Col) {
  FF = rep("allOthers",length(V(g)$name))
  names(FF) = V(g)$name
  for (i in seq_along(Labels)) {
    x = Labels[[i]]
    x = x[x %in% names(FF)]
    FF[x] = names(Labels)[i]
  }
  FF = factor(FF, levels=names(Labels))
  
  P = cbind(tapply(co[,1], FF, mean), tapply(co[,2], FF, mean))
  P2 = P * 1.06 / sqrt(apply(P^2,1,sum))
  for (i in 1:nrow(P)) {
    text(P2[i,1],P2[i,2],row.names(P2)[i],adj=(cbind(-sign(P2[i,1]),-sign(P2[i,2]))+1)/2,col=Col[i])
    lines(x=c(P2[i,1],P[i,1]),c(P2[i,2],P[i,2]),col=Col[i],lwd=2)
  }
}

