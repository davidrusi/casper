\name{probNonEquiv}
\alias{probNonEquiv}
\alias{probNonEquiv,ExpressionSet-method}
\alias{probNonEquiv,list-method}
\alias{pvalTreat}
\alias{pvalTreat,ExpressionSet-method}
\alias{pvalTreat,list-method}
\title{
  \code{probNonEquiv} performs a Bayesian hypothesis test for equivalence between group means.
  It returns the posterior probability that |mu1-mu2|>logfc.
  
  \code{pvalTreat} is a wrapper to \code{treat} in package \code{limma},
  which returns P-values for the same hypothesis test.
}
\description{
  \code{probNonEquiv} computes v_i=P(|theta_i| > logfc | data), where theta_i is
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
probNonEquiv(x, groups, logfc = log(2), minCount, method = "plugin", mc.cores=1)

pvalTreat(x, groups, logfc = log(2), minCount, p.adjust.method='none', mc.cores = 1) 
}
\arguments{
\item{x}{ExpressionSet containing expression levels, or list of ExpressionSets}
\item{groups}{Variable in fData(x) indicating the two groups to
  compare (the case with more than 2 groups is not implemented).}
\item{logfc}{Biologically relevant threshold for the log fold change,
  i.e. difference between groups means in log-scale}
\item{minCount}{ If specified, probabilities are only computed for rows with \code{fData(x)$readCount >= minCount}}
\item{method}{ Set to \code{'exact'} for exact posterior probabilities
  (slower), \code{'plugin'} for plug-in approximation (much faster).
  Typically both give very similar results.}
\item{mc.cores}{Number of parallel processors to use. Ignored unless
  \code{x} is a list.}
\item{p.adjust.method}{P-value adjustment method, passed on to \code{p.adjust}}
}
\value{
  If \code{x} is a single \code{ExpressionSet}, \code{probNonEquiv}
  returns a vector with posterior probabilities
  (NA for rows with less than \code{minCount} reads).
  \code{pvalTreat} returns TREAT P-values instead.
  
  If \code{x} is a list of \code{ExpressionSet}, the function is applied
  to each element separately and results are returned as columns in the
  output matrix.
}
\seealso{
\code{treat} in package \code{limma}, \code{p.adjust}
}
\references{
  Rossell D, Stephan-Otto Attolini C, Kroiss M, Stocker A. Quantifying
  Alternative Splicing from Paired-End RNA-sequencing data. Annals of
  Applied Statistics, 8(1):309-330

  McCarthy DJ, Smyth GK. Testing significance relative to a fold-change
  threshold is a TREAT. Bioinformatics, 25(6):765-771
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
  truth <- rep(c(FALSE,TRUE),c(p-11,11))
  getRoc(truth, pp>.9)
  getRoc(truth, pp>.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ htest }
