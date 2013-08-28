\name{mergeBatches}
\alias{mergeBatches}
\alias{mergeBatches,ExpressionSet,ExpressionSet-method}
\alias{mergeBatches,ExpressionSet,simulatedSamples-method}
\title{
  Merge two ExpressionSet objects by doing quantile normalization
  and computing partial residuals (i.e. substracting group mean
  expression in each batch). As currently implemented the method is only
  valid for balanced designs, e.g. each batch has the same number of samples 
  per group.
}
\description{
  \code{mergeBatches} combines \code{x} and \code{y} into an
  \code{ExpressionSet},
  performs quantile normalization and adjusts for batch effects by
  subtracting the mean expression in each batch (and then adding the
  grand mean so that the mean expression per gene is unaltered).
}
\usage{
mergeBatches(x, y, mc.cores=1)
}
\arguments{
  \item{x}{\code{ExpressionSet} object with data from batch 1.}
  \item{y}{Either \code{ExpressionSet} object with data from batch 2, or
  \code{simulatedSamples} object with data from multiple simulations.}
  \item{mc.cores}{Number of processors to be used (ignored when \code{y}
  is an \code{ExpressionSet})}
}
\value{
  When \code{y} is an \code{ExpressionSet}, \code{mergeBatches} returns
  an \code{ExpressionSet} with combined expressions.
  Its \code{featureData} contains a variable \code{"batch"} indicating
  the batch that each sample corresponded to.

  When \code{y} is a \code{simulatedSamples} object, \code{mergeBatches}
  is applied to combine \code{x} with each dataset in \code{y} and a
  list of \code{ExpressionSet} objects is returned.
}
\author{
  David Rossell
}
\examples{
#Fake data from 2 batches
x <- matrix(rnorm(6),nrow=2)
colnames(x) <- paste('x',1:3,sep='')
y <- matrix(1+rnorm(6),nrow=2)
colnames(y) <- paste('y',1:3,sep='')
x <- new("ExpressionSet",exprs=x)
y <- new("ExpressionSet",exprs=y)
exprs(x)
exprs(y)

#Merge & adjust
z <- mergeBatches(x,y)
exprs(z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
