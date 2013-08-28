\name{getRoc}
\alias{getRoc}
\alias{getRoc,logical,logical-method}
\alias{getRoc,numeric,numeric-method}
\alias{getRoc,matrix,matrix-method}
\title{
  Operating characteristics of differential expression analysis
}
\description{
  \code{getRoc} compares simulation truth and data analysis results to
  determine False Positives (FP), False Negatives (FP), True Positives
  (TP), True Negatives (TN), Positives (FP+TP), False Discovery
  Proportion (FP/P) and Power (TP/(TP+FN)). 
}
\usage{
getRoc(simTruth, decision) 
}
\arguments{
  \item{simTruth}{Binary vector or matrix indicating simulation truth
    (\code{FALSE} or 0 for non differential expression, \code{TRUE} or 1
    for differential expression)}
  \item{decision}{Binary vector or matrix with differential expression
    calls based on some data analysis.}
}
\section{Methods}{
  \describe{
    \item{\code{signature(simTruth='logical',decision='logical')}}{
      Operating characteristics are computed for a single simulation
    }
    \item{\code{signature(simTruth='numeric',decision='numeric')}}{
      Operating characteristics are computed for a single simulation
    }
    \item{\code{signature(simTruth='matrix',decision='matrix')}}{
      \code{simTruth} and \code{decision} contain truth and calls for
      several simulations (in columns).
      \code{getRoc} returns a data.frame with operating characteristics
      in each simulation.
    }
  }
}
\value{
\code{data.frame} with TP, FP, TN, FN, P, FDR and Power.
}
\author{
  David Rossell
}
\examples{
## See help(probNonEquiv) for an example
}
\keyword{ htest }% __ONLY ONE__ keyword per line

