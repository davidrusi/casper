mergeBatches <- function(x, y) {
  #Merge ExpressionSets x and y, apply quantile norm and adjust for batch effects via partial residuals
  if ((class(x) != 'ExpressionSet') | (class(y) != 'ExpressionSet')) stop('x and y must be ExpressionSet objects')
  xnew <- combine(x,y)
  xnew$batch <- factor(rep(c('batch1','batch2'),c(ncol(x),ncol(y))))
  xnewnorm <- quantileNorm(xnew)
  xnewadj <- anovaAdjustment(xnewnorm, adjustmentVariable='batch')
  return(xnewadj)
}


### AUXILIZARY FUNCTIONS: QUANTILE NORMALIZATION & ANOVA ADJUSTMENT
quantileNorm <- function(x) { UseMethod("quantileNorm") }

quantileNorm.ExpressionSet <- function(x) { exprs(x) <- quantileNorm(exprs(x)); return(x) }

quantileNorm.matrix <- function(x) {
  probsx <- seq(0,1,length=nrow(x))
  qx <- quantile(x,probs=probsx)
  f <- approxfun(probsx,qx)
  for (i in 1:ncol(x)) {
    probsy <- rank(x[,i])/(nrow(x)+1)
    x[,i] <- f(probsy)
  }
  return(x)
}


anovaAdjustment <- function(x,adjustmentVariable,na.rm=FALSE) {
  #Perform ANOVA adjustment, i.e. remove effect of nuisance covariates via lmFit. This can be used to remove experimentation date or batch processing effects on gene expression.
  # - x: ExpressionSet with expression values to be adjusted
  # - adjustmentVariable: name of adjustment variable, which must be contained in pData(x).
  # Output: ExpressionSet with adjusted expression values
  if (missing(adjustmentVariable)) stop('adjustmenVariable must be specified')
  UseMethod("anovaAdjustment")
}


anovaAdjustment.ExpressionSet <- function(x,adjustmentVariable,na.rm=FALSE) {
  if (!(adjustmentVariable %in% names(pData(x)))) { stop('Specified variable not found in pData') }
  design <- model.matrix(~ -1 + as.factor(pData(x)[,adjustmentVariable]))
  beta <- exprs(x) %*% design %*% diag(1/colSums(design))
  linpred <- beta %*% t(design)
  exprs(x) <- exprs(x) - linpred + rowMeans(exprs(x),na.rm=na.rm)
  return(x)
}
