# Compute Prob( |theta| > logfc | data), where theta is difference between group means (in log scale)
# - x is an expressionSet (w/ simulated expressions, in our case)
# - groups: variable in fData(x)
# - logfc is a biologically relevant threshold for the log-FC
# - minCount: probabilities are only computed for genes with fData(x)$readCount >= minCount
# - method: 'exact' for exact posterior probabilities (slower), 'plugin' for plug-in approximation (much faster)

setMethod("probNonEquiv", signature(x='list'), function(x, groups, logfc=log(2), minCount, method="plugin", mc.cores=1) {
  if (any(sapply(x,class) != 'ExpressionSet')) stop('All elements in x must be of class ExpressionSet')
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) {
      ans <- parallel::mclapply(x, probNonEquiv, groups=groups, logfc=logfc, minCount=minCount, method=method, mc.cores=mc.cores, mc.preschedule=TRUE)
    } else {
      stop("parallel package has not been loaded!")
    }
  } else {
    ans <- lapply(x, probNonEquiv, groups=groups, logfc=logfc, minCount=minCount, method=method)
  }
  do.call(cbind, ans)
}
)

setMethod(probNonEquiv, signature(x="ExpressionSet"), function(x, groups, logfc=log(2), minCount, method="plugin", mc.cores=1) {
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
    ans[sel] <- pnorm(-logfc,m[,1]-m[,2],sd=sqrt(v[,1]+v[,2])) + pnorm(logfc,m[,1]-m[,2],sd=sqrt(v[,1]+v[,2]),lower.tail=FALSE)
    #ans[sel] <- pt((-logfc-(m[,1]-m[,2]))/sqrt(v[,1]+v[,2]), df=ncolsel[1]) + pt((logfc-(m[,1]-m[,2]))/sqrt(v[,1]+v[,2]), df=ncolsel[2], lower.tail=FALSE)
  }
  return(ans)
} 
)

dinvgammaVec <- function(x, alpha, beta, logscale=TRUE) {
    ans <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
    if (!logscale) ans <- exp(ans)
    return(ans)
}


setMethod("pvalTreat", signature(x='list'), function(x, groups, logfc=log(2), minCount, p.adjust.method='none', mc.cores=1) {
  if (any(sapply(x,class) != 'ExpressionSet')) stop('All elements in x must be of class ExpressionSet')
  if (mc.cores>1) {
    if ('parallel' %in% loadedNamespaces()) {
      ans <- parallel::mclapply(x, pvalTreat, groups=groups, logfc=logfc, minCount=minCount, p.adjust.method=p.adjust.method, mc.cores=mc.cores, mc.preschedule=TRUE)
    } else {
      stop("parallel package has not been loaded!")
    }
  } else {
    ans <- lapply(x, pvalTreat, groups=groups, logfc=logfc, minCount=minCount, p.adjust.method=p.adjust.method)
  }
  do.call(cbind, ans)
}
)

setMethod("pvalTreat", signature(x='ExpressionSet'), function(x, groups, logfc=log(2), minCount, p.adjust.method='none', mc.cores=1) {
  if (is.null(pData(x)[,groups])) stop(paste('Variable',groups,'not found in pData(x)'))
  if (length(unique(pData(x)[,groups])) != 2) stop('pvalTreat only works for 2 groups')
  if (!missing(minCount)) sel <- fData(x)$readCount >= minCount else sel <- rep(TRUE,nrow(x))
  ans <- rep(NA,nrow(x)); names(ans) <- featureNames(x)
  fit <- contrasts.fit(lmFit(x[sel,], design=model.matrix(~ pData(x)[,'group'])), contrasts=c(0,1))
  ans[sel] <- p.adjust(treat(fit, lfc=logfc)$p.value[,1], method=p.adjust.method)
  return(ans)
}
)

