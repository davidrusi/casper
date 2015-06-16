setMethod(qqnormGenomeWide, signature(x="ExpressionSet"), function(x, ngenes=min(1000,nrow(x)), ...) {
  qqnormGenomeWide(exprs(x))
}
)

setMethod(qqnormGenomeWide, signature(x="data.frame"), function(x, ngenes=min(1000,nrow(x)), ...) {
  qqnormGenomeWide(as.matrix(x))
}
)

setMethod(qqnormGenomeWide, signature(x="matrix"), function(x, ngenes=min(1000,nrow(x)), ...) {
  if (ngenes>nrow(x)) ngenes <- nrow(x)
  x <- x[1:ngenes,]
  z <- t(scale(t(x),center=TRUE,scale=TRUE))
  zsort <- t(apply(z,1,function(z) z[order(z)]))
  q <- qnorm(ppoints(1:ncol(z))) #theoretical quantiles
  plot(q,zsort[1,], xlab='Theoretical quantile', ylab='Observed quantile', ylim= range(z), pch='.',...)
  if (ngenes>1) {
    for (i in 2:ngenes) points(q,zsort[i,],pch='.')
    abline(0,1)
    points(q,colMeans(zsort),col=1,pch=16,cex=1.5)
  }
}
)


setMethod(qqgammaGenomeWide, signature(x="ExpressionSet"), function(x, ngenes=min(1000,nrow(x)), ...) {
  qqgammaGenomeWide(exprs(x))
}
)

setMethod(qqgammaGenomeWide, signature(x="data.frame"), function(x, ngenes=min(1000,nrow(x)), ...) {
  qqgammaGenomeWide(as.matrix(x))
}
)

setMethod(qqgammaGenomeWide, signature(x="matrix"), function(x, ngenes=min(1000,nrow(x)), ...) {
  x <- x[1:ngenes,,drop=FALSE]
  m <- rowMeans(x); v <- rowVar(x)
  lest <- m/v; aest <- m^2/v
  xsort <- t(apply(x,1,function(z) z[order(z)]))
  pp <- ppoints(1:ncol(xsort))
  q <- do.call(rbind,lapply(1:nrow(x), function(i) qgamma(pp, aest[i], lest[i]))) #theoretical quantiles
  plot(q[1,],xsort[1,], xlab='Theoretical quantile', ylab='Observed quantile', ylim= range(xsort), xlim=range(xsort), pch='.', ...)
  if (ngenes>1) { for (i in 2:ngenes) points(q[i,],xsort[i,],pch='.') }
  abline(0,1)
}
)
