\name{quantileNorm}
\alias{quantileNorm}
\alias{quantileNorm,ExpressionSet-method}
\alias{quantileNorm,matrix-method}
\title{
  Apply quantile normalization
}
\description{
  Perform quantile normalization on the columns of a matrix or ExpressionSet
}
\usage{
quantileNorm(x)
}
\arguments{
  \item{x}{\code{ExpressionSet} or \code{matrix}}
}
\value{
  Returns \code{x} with quantile normalized columns
}
\author{
  David Rossell
}
\examples{
x <- cbind(rnorm(1000),rnorm(1000,2,4))
boxplot(x)

xnorm <- quantileNorm(x)
boxplot(xnorm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
