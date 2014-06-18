\name{distrGSE37704}
\alias{distrGSE37704}
\docType{data}
\title{
  Estimated read start and insert size distributions from MiSeq data
  in GEO dataset GSE37704.
}
\description{
We downloaded the fastq files, aligned with TopHat and processed with
\code{wrapKnown} to obtain the estimated distributions for each of the 6
samples. \code{distrsMiSeq} is a list with the 6 corresponding elements.
The estimated distributions for HiSeq data were very similar,
hence these distributions can be used as defaults for Illumina MiSeq and
HiSeq experiments.
}
\usage{data(distrsMiSeq)}
\format{
  An \code{list} with 6 elements of class \code{readDistrs}.
  See \code{help(getDistrs)} and \code{help(readDistrs-class)} for details.
}
\examples{
data(distrGSE37704)
distrGSE37704
plot(distrGSE37704[[1]],'readSt')
plot(distrGSE37704[[1]],'fragLength')
}
\keyword{datasets}
