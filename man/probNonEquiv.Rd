\name{probNonEquiv}
\alias{probNonEquiv}
\title{
  Posterior probability that differences between group means are greater
  than a given threshold, i.e. that groups are non-equivalent by a
  biologically meaningful amount.
}
\description{
  The function computes v_i=P(|theta_i| > logfc | data), where theta_i is
  the difference between group means for gene i. This posterior
  probability is based on the NNGCV model from package EBarrays, which
  has a formulation similar to limma in an empirical Bayes framework.
  Notice that the null hypothesis here is that |theta_i|<logfc,
  e.g. isoforms with small fold changes are regarded as uninteresting.

  Subsequent differential expression calls are based on selecting large
  v_i. For instance, selecting v_i >= 0.95 guarantees that the posterior
  expected false discovery proportion (a Bayesian FDR analog) is below 0.05.
}
\usage{
probNonEquiv(x, groups, logfc = log(2), minCount, method = "exact")
}
\arguments{
\item{x}{ExpressionSet containing expression levels}
\item{groups}{Variable in fData(x) indicating the two groups to
  compare (the case with more than 2 groups is not implemented).}
\item{logfc}{Biologically relevant threshold for the log fold change,
  i.e. difference between groups means in log-scale}
\item{minCount}{ If specified, probabilities are only computed for rows with \code{fData(x)$readCount >= minCount}}
\item{method}{ Set to \code{'exact'} for exact posterior probabilities
  (slower), \code{'plugin'} for plug-in approximation (much faster).
  The latter is less conservative than the former, but typically controls
  the FDR below alpha when setting a threshold 1-alpha on the resulting
  posterior probabilities.}
}
\value{
  Vector with posterior probabilities (NA's are returned for rows with
  less than \code{minCount} reads).
}
\references{
Rossell D, Stephan-Otto Attolini C, Kroiss M, Stocker A. Quantifying Alternative Splicing from Paired-End
RNA-sequencing data. COBRA (http://biostats.bepress.com/cobra/art97)
}
\author{
  Victor Pena, David Rossell
}
\examples{

  #Simulate toy data
  p <- 50; n <- 10
  x <- matrix(rnorm(p*2*n),nrow=p)
  x[(p-10):p,1:n] <- x[(p-10):p,1:n] + 1.5
  x <- new("ExpressionSet",exprs=x)
  x$group <- rep(c('group1','group2'),each=n)

  #Posterior probabilities
  pp <- probNonEquiv(x, groups='group', logfc=0.5)
  d <- rowMeans(exprs(x[,1:n])) - rowMeans(exprs(x[,-1:-n]))
  plot(d,pp,xlab='Observed log-FC')
  abline(v=c(-.5,.5))

  #Check false positives
  truth <- rep(0:1,c(p-11,11))
  table(truth, pp>.5)
  table(truth, pp>.9)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ htest }
