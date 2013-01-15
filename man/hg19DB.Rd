\name{hg19DB}
\alias{hg19DB}
\docType{data}
\title{
Subset of human genome (UCSC hg19 version)
}
\description{
We downloaded the human genome hg19 via \code{procGenome}
and selected a few genes from chromosome 1
to use as a toy data for the vignette and examples. 
}
\usage{data(hg19DB)}
\format{
An \code{annotatedGenome} object. See \code{help(procGenome)} and
\code{help(annotatedGenome-class)} for details.
}
\examples{
data(hg19DB)
hg19DB
slotNames(hg19DB)
}
\keyword{datasets}
