\name{distrsGSE37704}
\alias{distrsGSE37704}
\docType{data}
\title{
  Estimated read start and insert size distributions from MiSeq data
  in GEO dataset GSE37704.
}
\description{
We downloaded the fastq files, aligned with TopHat and processed with
\code{wrapKnown} to obtain the estimated distributions for each of the 6
samples. \code{distrsGSE37704} is a list with the 6 corresponding elements.
The estimated distributions for HiSeq data were very similar,
hence these distributions can be used as defaults for Illumina MiSeq and
HiSeq experiments.
}
\usage{data(distrsGSE37704)}
\format{
  An \code{list} with 6 elements of class \code{readDistrs}.
  See \code{help(getDistrs)} and \code{help(readDistrs-class)} for details.
}
\examples{
data(distrsGSE37704)
distrsGSE37704
plot(distrsGSE37704[[1]],'readSt')
lines(distrsGSE37704[[2]], 'readSt', col=2)
plot(distrsGSE37704[[1]],'fragLength')
}
\keyword{datasets}
