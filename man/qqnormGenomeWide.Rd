\name{qqnormGenomeWide}
\alias{qqnormGenomeWide}
\alias{qqnormGenomeWide,data.frame-method}
\alias{qqnormGenomeWide,ExpressionSet-method}
\alias{qqnormGenomeWide,matrix-method}
\title{
  Genome-wide qq-normal plots
}
\description{
  Overlays quantile-quantile normal plots (qqnorm) for a series of genes
  (rows in the input matrix), to provide an overall assessment of Normality.
}
\usage{
qqnormGenomeWide(x, ngenes=min(1000, nrow(x)), ...)
}
\arguments{
  \item{x}{\code{ExpressionSet}, \code{matrix} or \code{data.frame}
    with genes/isoforms in rows}
  \item{ngenes}{A qqnorm plot is produced for the first \code{ngenes}
    rows in \code{x}}
  \item{...}{Other arguments to be passed on to code{plot}}
}
\value{
  Produces a figure overlaying qq-normal plots for \code{ngenes}
  comparing observed vs. theoretical quantiles
}
\author{
  David Rossell
}
\examples{
mu <- rnorm(100)
x <- matrix(rnorm(100*5,mu),ncol=5)

qqnormGenomeWide(x)
}
\keyword{ hplot }
\keyword{distribution}
