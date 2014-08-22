\name{simMAE}
\alias{simMAE}
\docType{data}
\title{
  Simulate Mean Absolute Error (MAE) in estimating isoform expression
  under various experimental settings.
}
\description{
  Simulate several future RNA-seq data under various experimental
  settings (sequencing depth, read length, insert sizes), estimate
  isoform expression and assess the MAE incurred in the estimation
  process.
  The function is a wrapper combining functions simReads and calcExp.
}
\usage{simMAE(nsim, islandid, nreads, readLength, fragLength,
burnin=1000, pc, distr, readLength.pilot, usePilot=FALSE,
retTxsError=FALSE, genomeDB, mc.cores=1, mc.cores.int=1, verbose=FALSE)}

\arguments{
  \item{nsim}{Number of RNA-seq datasets to generate (often as little as
    \code{nsim=10} suffice)}
  \item{islandid}{When specified this argument indicates to run the
  simulations only for gene islands with identifiers in
  \code{islandid}. When not specified genome-wide simulations are
  performed.}
\item{nreads}{Vector indicating the number of read pairs for each
  experimental setting}
\item{readLength}{Vector indicating the read length in each experimental
  setting}
\item{fragLength}{Vector indicating the mean insert size in each
  experimental setting}
\item{burnin}{Number of MCMC burn-in samples (passed on to \code{calcExp})}
\item{pc}{Observed path counts in pilot data}
\item{distr}{Estimated read start and insert size distributions in pilot
  data}
\item{readLength.pilot}{Read length in pilot data}
\item{usePilot}{By default \code{casper} assumes that the pilot data is from a related experiment rather than the current tissue of interest (\code{usePilot=FALSE}). Hence, the pilot data is used to simulate new RNA-seq data but not to estimate its expression. However, in some cases we may be interested in re-sequencing the pilot sample at deeper length, in which case one would want to combine the pilot data with the new data to obtain more precise estimates. This can be achieved by setting \code{usePilot=TRUE}}
\item{retTxsError}{If \code{retTxsError=TRUE}, \code{simMAE} returns posterior expected MAE
for each individual isoform.
Else the output is a \code{data.frame} with overall MAE across all isoforms}
\item{genomeDB}{\code{annotatedGenome} object, as returned by \code{procGenome}}
\item{mc.cores}{Number of cores to use in the expression estimation step, passed on to \code{calcExp}}
\item{mc.cores.int}{Number of cores to simulate RNA-seq datasets in parallel}
\item{verbose}{Set \code{verbose=TRUE} to print progress information}
}
\details{
Ideally \code{simMAE} should use pilot data from a relevant related experiment
to simulate what future data may look like for the current experiment of
interest.
The recommened way to do this is to download a .bam file from such a
related experiment and processing it in casper with function
\code{wrapKnown}.
The object output by \code{wrapKnown} is a list with elements
named \code{'pc'}, \code{'distr'} which can be given as input to \code{simMAE}.
}
\value{
If \code{retTxsError==TRUE}, \code{simMAE} returns posterior expected MAE
  for each individual isoform.
Else the output is a \code{data.frame} with overall MAE across all isoforms  
}
\seealso{wrapKnown,simReads,calcExp}
\examples{
## maybe str(simMAE) ; plot(simMAE) ...
}
\keyword{ datagen }

