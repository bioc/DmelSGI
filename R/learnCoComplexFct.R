
learnCoComplexFct <- function(C, ProteinComplexes) {
  G = matrix(FALSE,nrow=nrow(C),ncol=nrow(C))
  row.names(G) = colnames(G) = row.names(C)
  
  for (i in 1:length(ProteinComplexes)) {
    g = ProteinComplexes[[i]]$gene_id
    G[g,g] = TRUE
  }
  
  X = matrix(1:nrow(C),nrow=nrow(C),ncol=nrow(C))
  Y = t(X)
  
  Complex = C[(X > Y) & G]
  NonComplex = C[(X > Y) & !G]
  
  x = seq(-1,1,by=0.01)
  p = rep(0.0, length(x))
  for (i in 1:length(x)) {
    s1 = sum(Complex >= x[i])
    s2 = sum(NonComplex >= x[i])
    p[i] = s1 / (s1+s2)
  }
  p[is.na(p)] = 1.0
  
  coComplexFct = list(breaks = x, p=p, Complex=Complex, NonComplex=NonComplex)
  coComplexFct
}

convertCorrelations <- function(C, coComplexFct) {
  P = C
  #   P[] = approx(coComplexFct$breaks,coComplexFct$p,C[])
  P[] = cut(C[],breaks=coComplexFct$breaks,label=FALSE)
  P[] = coComplexFct$p[P[]]
  P
}
