\name{getReads}
\alias{getReads}
\alias{getReads,procBam-method}
\docType{methods}
\title{
getReads returns the reads stored in a \code{procBam} object.
}
\description{
\code{procBam} objects store reads that have been split according to
their CIGAR codes. \code{getReads} accesses these reads.
}
\usage{
getReads(x)
}
\arguments{
 \item{x}{Object of class \code{procBam}}
}
\section{Methods}{
\describe{
  \item{\code{signature(x='procBam')}}{
    Return reads stored in \code{x}.
  }
}}
\value{\code{RangedData} object with reads stored in \code{x}.}
\examples{
#See example in calcExp
}
\keyword{manip}
