
mymedpolish <- function (x, TP, TemplateNeg, QueryNeg, eps = 0.0001, maxiter = 100, na.rm = TRUE) {
  MTP = max(TP)
    z <- x
    nr <- nrow(z)
    nc <- ncol(z)
    t <- 0
    r <- numeric(nr)
    c <- matrix(numeric(MTP*nc),nrow=MTP)
    oldsum <- 0
    for (iter in 1L:maxiter) {
        rdelta <- apply(z, 1L, median, na.rm = na.rm)
        z <- z - matrix(rdelta, nrow = nr, ncol = nc)
        r <- r + rdelta
        cdelta <- apply(z, 2L, function(s) { tapply(s,TP,median,na.rm=na.rm) })
        z <- z - apply(cdelta,2,function(s) { s[TP] } )
        c <- c + cdelta
        if(is.null(QueryNeg)) {
          delta = median(c, na.rm = na.rm)
          if (!is.finite(delta)) {
            delta = 0.0
          }
        } else {
          delta <- mean(c[,QueryNeg],na.rm=TRUE)
          if (!is.finite(delta)) {
            delta = median(c, na.rm = na.rm)
            if (!is.finite(delta)) {
              delta = 0.0
            }
          }
        }
        c <- c - delta
        t <- t + delta
        if(is.null(TemplateNeg)) {
          delta = median(r, na.rm = na.rm)
          if (!is.finite(delta)) {
            delta = 0.0
          }
        } else {
          delta <- mean(r[TemplateNeg],na.rm=TRUE)
          if (!is.finite(delta)) {
            delta = median(r, na.rm = na.rm)
            if (!is.finite(delta)) {
              delta = 0.0
            }
          }
        }
        r <- r - delta
        t <- t + delta
        newsum <- sum(abs(z), na.rm = na.rm)
        converged <- newsum == 0 || abs(newsum - oldsum) < eps * newsum
        if (converged) 
            break
        oldsum <- newsum
    }
    if (!converged) {
      warning(gettextf("medpolish() did not converge in %d iterations", 
        maxiter), domain = NA)
    }
    ans <- list(neg = t, templateMainEffect = r, queryMainEffect = t(c), pi = z)
    ans
}


