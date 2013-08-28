# Compute Prob( |theta| > logfc | data), where theta is difference between group means (in log scale)
# - x is an expressionSet (w/ simulated expressions, in our case)
# - groups: variable in fData(x)
# - logfc is a biologically relevant threshold for the log-FC
# - minCount: probabilities are only computed for genes with fData(x)$readCount >= minCount
# - method: 'exact' for exact posterior probabilities (slower), 'plugin' for plug-in approximation (much faster)
probNonEquiv <- function(x, groups, logfc=log(2), minCount, method='exact') {
  #require(plyr)
  if (class(x) != 'ExpressionSet') stop('x must be of class ExpressionSet')
  if (missing(groups)) stop('groups must be specified')
  if (is.null(pData(x)[,groups])) stop(paste('Variable',groups,'not found in pData(x)'))
  if (!missing(minCount)) sel <- fData(x)$readCount >= minCount else sel <- rep(TRUE,nrow(x))
  ans <- rep(NA,nrow(x)); names(ans) <- rownames(x)
  x <- x[sel,]
  fit <- fitNNSingleHyp(x, groups=groups, B=5, trace=FALSE) 
  groups <- pData(x)[,groups]
  x <- exprs(x)  
  par <- fit$parest
  # a.sigma, b.sigma
  a.sigma <- rep(.5*(par['v0']+ncol(x)), nrow(x))
  b.sigma <- rep(.5*par['v0']*par['sigma02'], nrow(x))
  ncolsel <- integer(length(unique(groups)))
  xbar <- matrix(NA, nrow=nrow(x), ncol=length(ncolsel))
  names(ncolsel) <- colnames(xbar) <- unique(groups)
  for (j in unique(groups)) {
    colsel <- groups %in% j
    ncolsel[j] <- sum(colsel)
    xbar[,j] <- rowMeans(x[,colsel])
    b.sigma <- b.sigma + .5*rowSums((x[,colsel]-xbar[,j])^2)
  }
  # Compute integral
  if (method=='exact') { 
    f <- function(phi) {
      v <-  matrix(nrow=length(phi), ncol=length(unique(groups))); m <-  v
      k <-  1
      # m1, m2, v1, v2 for isoform i
      for (j in unique(groups)) {
        v[,k] <- 1/(ncolsel[j]/phi + 1/par['tau02']) 
        m[,k] <- (ncolsel[j]*xbar[i,j]/phi + par['mu0']/par['tau02']) * v[,k]     
        k <- k+1
      }
      (pnorm(-logfc,m[,1]-m[,2],sd=sqrt(v[,1]+v[,2])) + pnorm(logfc, m[,1]-m[,2],sd=sqrt(v[,1]+v[,2]),lower.tail=F))*dinvgammaVec(phi,a.sigma[i],b.sigma[i],logscale=F)
    }
    sel <- which(sel)
    for (i in 1:length(sel)) ans[sel[i]] <- integrate(f, 0, Inf)$value
  } else {
    phimode <- ifelse(a.sigma>1, b.sigma/(a.sigma-1), b.sigma/a.sigma)
    v <- cbind(1/(ncolsel[1]/phimode + 1/par['tau02']), 1/(ncolsel[2]/phimode + 1/par['tau02']))
    m <- (t(ncolsel*t(xbar))/phimode + par['mu0']/par['tau02']) * v
    #ans[sel] <- pnorm(-logfc,m[,1]-m[,2],sd=sqrt(v[1]+v[2])) + pnorm(logfc,m[,1]-m[,2],sd=sqrt(v[1]+v[2]),lower.tail=FALSE)
    ans[sel] <- pt((-logfc-(m[,1]-m[,2]))/sqrt(v[1]+v[2]), df=ncolsel[1]) + pt((logfc-(m[,1]-m[,2]))/sqrt(v[1]+v[2]), df=ncolsel[2], lower.tail=FALSE)    
  }
  return(ans)
} 


dinvgammaVec <- function(x, alpha, beta, logscale=TRUE) {
    ans <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
    if (!logscale) ans <- exp(ans)
    return(ans)
}



probNonEquiv2 <- function(x, groups, logfc=log(2), minCount) {
  #same as postProb, uses quantile integration
  require(gaga)
  require(plyr)
  if (missing(groups)) stop('groups must be specified')
  if (is.null(pData(x)[,groups])) stop(paste('Variable',groups,'not found in pData(x)'))
  if (!missing(minCount)) sel <- fData(x)$readCount >= minCount else sel <- rep(TRUE,nrow(x))
  ans <- double(nrow(x))
  x <- x[sel,]
  fit <- fitNNSingleHyp(x, groups=groups, B=5, trace=FALSE) 
  groups <- pData(x)[,groups]
  x <- exprs(x)  
  par <- fit$parest
  # a.sigma, b.sigma
  a.sigma <- rep(.5*(par['v0']+ncol(x)), nrow(x))
  b.sigma <- rep(.5*par['v0']*par['sigma02'], nrow(x))
  ncolsel <- integer(length(unique(groups)))
  xbar <- matrix(NA, nrow=nrow(x), ncol=length(ncolsel))
  names(ncolsel) <- colnames(xbar) <- unique(groups)
  for (j in unique(groups)) {
    colsel <- groups %in% j
    ncolsel[j] <- sum(colsel)
    xbar[,j] <- rowMeans(x[,colsel])
    b.sigma <- b.sigma + .5*rowSums((x[,colsel]-xbar[,j])^2)
  }
  # function to integrate
  f <- function(phi) {
    v <- 1/(ncolsel/phi + 1/par['tau02'])
    m <- t((ncolsel*t(xbar)/phi + par['mu0']/par['tau02']) * v)
    pp <- log(pnorm(-logfc, m[,1]-m[,2], sd=sqrt(v[1]+v[2])) + pnorm(logfc, m[,1]-m[,2], sd=sqrt(v[1]+v[2]), lower.tail=FALSE))
    #exp(pp + dinvgammaVec(phi, a.sigma, b.sigma, logscale=TRUE))
    exp(pp + dinvgammaVec(phi, a.sigma, b.sigma, logscale=TRUE) - dinvgammaVec(phi, aprop, bprop, logscale=TRUE))
  }
  B <- 1000
  aprop <- a.sigma * min(b.sigma) / mean(b.sigma); bprop <- min(b.sigma)
  qphi <- rev(1/qgamma((1:B)/(B+1), aprop, bprop))
  ans[sel] <- rowMeans(do.call(cbind, lapply(qphi, f)))
  #dprop <- dinvgammaVec(qphi, aprop, bprop, logscale=FALSE)
  #ans[sel] <- colMeans(do.call(rbind, lapply(qphi, f)) / dprop)
  return(ans)
} 
