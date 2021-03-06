\name{simMultSamples}
\alias{simMultSamples}
\alias{casperDesign}
\title{
  Simulate paired end reads for multiple future samples
  based on pilot data, and obtain their expression estimates via casper
}
\description{
  Simulate true expression levels and observed data (casper expression
  estimates) for future samples within each group.
  
  These simulations serve as the basis for sample size calculation:
  if one were to sequence \code{nsamples} new RNA-seq samples, what data
  would we expect to see? The simulation is posterior predictive,
  i.e. based on the current available data \code{x}.
}
\usage{
simMultSamples(nsim, nsamples, nreads, readLength, fragLength, x,
groups='group', distrs, genomeDB, model='LNNMV', verbose=TRUE, mc.cores=1)
}
\arguments{
  \item{nsim}{Number of simulations to obtain}
  \item{nsamples}{Vector indicating number of future samples per group,
    e.g. \code{nsamples=c(5,5)} to simulate 5 new samples for 2 groups.}
  \item{nreads}{Desired number of paired-end reads per
    sample. The actual number of aligned reads for any given sample
  differs from this amount, see details.}
  \item{readLength}{Read length, i.e. in an experiment with paired reads
    at 100bp each, \code{readLength=100}.}
  \item{fragLength}{Desired average insert size (size of RNA molecules
  after fragmentation). If missing, insert sizes are as obtained from \code{distrs}}
  \item{x}{\code{ExpressionSet} containing pilot data. \code{x[[group]]} indicates groups to be compared}
  \item{groups}{Name of column in \code{pData(x)} indicating the groups}
  \item{distrs}{Fragment start and length distributions. It can be
    either an object or a list of objects of class readDistrs.
  In the latter case, an element is chosen at random for each individual
  sample to consider uncertainty in these distributions.
  If not specified, it defaults to data(distrsGSE37704).}
\item{genomeDB}{annotatedGenome object}
  \item{model}{Set to \code{'LNNMV'} to simulate from log-normal normal
    with modified variance model (Yuan and Kendsiorski, 2006),
    or to \code{'GaGa'} to simulate from the
    GaGa model (Rossell, 2009). See details.}
  \item{verbose}{Set to \code{TRUE} to print progress}
  \item{mc.cores}{Number of cores to use in function.
  \code{mc.cores>1} requires package \code{parallel}}
}
\value{
  Object of class \code{simulatedSamples}, which extends a list
  of length \code{nsim}. See the class documentation for some helpful
  methods (e.g. coef, exprs, mergeBatches).
  Each element is itself a list containing an individual simulation.
  
  \item{simTruth}{\code{data.frame} indicating the mean and
  standard deviation of the Normal distribution used to generate data
  from each group}
  \item{simExpr}{\code{ExpressionSet} with Casper expression estimates,
  as returned by \code{calcExp}. \code{pData(simExpr)} indicates group
  information, and \code{fData(simExpr)} the number of simulated reads
  for each sample (in columns \code{'explCnts'}) and across all samples (in column \code{'readCount'})}
}
\details{
  The posterior predictive simulations is based on four steps:
  (1) simulate true expression for each group (mean and SD),
  (2) simulate true expression for future samples,
  (3) simulate paired reads for each future sample,
  (4) estimate expression from the reads via Casper.
  Below are some more details.

  1. Simulate true mean expression in each group and residual variance
  for each gene. If model=='LNNMV' this is based on the log-normal
  normal with modified variance model in package EBarrays
  (Yuan & Kendziorski 2006), if model=='GaGa' this is based on the GaGa
  model (Rossell, 2009).
  adapted to take into account that the expression estimates in the
  pilot data \code{x} are noisy (which is why \code{simMultSamples} requires the
  SE / posterior SD associated to \code{exprs(x)}).
  The simulated values are returned in component \code{"simTruth"} of
  the \code{simMultSamples} output.

  2. Simulate true isoform expression for each of the future samples. These are
  independent Normal draws with mean and variance generated in step
  1. True gene expression is derived from the isoform expressions.

  3. Determine the number of reads to be simulated for each gene based
  on its true expression (generated in step 2) and a Multinomial
  sampling model. For each sample:

  - The number of reads yielded by the experiment is Unif(.8*nreads,1.2*nreads)
  - A proportion of non-mappable reads is discarded using the power law in Li
  et al (2014)
  - Amongst remaining reads, we assume that a proportion Unif(0.6,0.9) were aligned
   (consistenly with reports from ENCODE project)
    
  The final number of simulated reads is reported in component \code{"simExpr"}
  of the \code{simMultSamples} output.
  
  4. Obtain expression estimates from the path counts produced in step 3
  via \code{calcExp}. These are reported in component \code{"simExpr"}
  of the \code{simMultSamples} output.
  
}
\author{
Victor Pena, David Rossell
}
\examples{
#Run casperDesign() to see full manual with examples
}
\references{
  Rossell D. (2009) GaGa: a Parsimonious and Flexible Model for
  Differential Expression Analysis. Annals of Applied Statistics, 3,
  1035-1051.

  Stephan-Otto Attolini C., Pena V., Rossell D. Bayesian designs for
  personalized alternative splicing RNA-seq studies (2015)

  Yuan, M. and Kendziorski, C. (2006). A unified approach for simultaneous gene
  clustering and differential expression identification. Biometrics, 62, 1089-1098.
}
\keyword{ datagen }
