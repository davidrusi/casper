\name{wrapKnown}
\alias{wrapKnown}
\title{
  Run all necessary steps to get expression estimates from multiple bam files with the casper pipeline.
}
\description{
  Function to analyze bam files to generate an ExpressionSet with
  expression estimates for all samples,
  read start and fragment length distributions, path counts and
  optinally processed reads.
}
\usage{
wrapKnown(bamFile, verbose=FALSE, seed=1, mc.cores.int=1,
mc.cores=1, genomeDB, readLength, rpkm=TRUE, priorq=2, priorqGeneExpr=2,
citype='none', niter=10^3, burnin=100, keep.pbam=FALSE,
keep.multihits=TRUE, chroms=NULL)
}
\arguments{
  \item{bamFile}{
    Names of bam files with the sample to analyze. These must sorted and
    indexed, and the index must be in the same directory.
}
  \item{verbose}{
    Set to \code{TRUE} to display progress information.
}
  \item{seed}{
    Set seed of random number generator.
}
  \item{mc.cores.int}{
    Number of cores to use when loading bam files. This is a memory intensive
    step, therefore number of cores must be chosen according to available RAM
    memory.
}
  \item{mc.cores}{
    Number of cores to use in expression estimation.
}
  \item{genomeDB}{
    \code{annotatedGenome} object containing annotated genome,
    as returned by the \code{procGenome} function.
}
  \item{readLength}{
    Read length in bp, e.g. in a paired-end experiment where
    75bp are sequenced on each end one would set \code{readLength=75}.
}
  \item{rpkm}{
    Set to \code{TRUE} to return reads per kilobase per million (RPKM).
    Set to \code{FALSE} to return relative expression levels.
    Important, relative expression adds up to 1 within gene island, NOT
    within gene. To get relative expressions within gene run
    \code{relexprByGene} afterwards. See details.   
}
  \item{priorq}{
    Parameter of the prior distribution on the proportion of reads
    coming from each variant. The prior is Dirichlet with prior sample
    size for each variant equal to priorq.
    We recommend \code{priorq=2} for estimation, as it pools the
    estimated expression away from 0 and 1 and returned lower
    estimation errors than \code{priorq=1} in our simulated
    experiments.
}
  \item{priorqGeneExpr}{
    Parameter for prior distribution on overall gene
    expression. Defaults to 2, which ensures non-zero estimates for
    all genes
}
  \item{citype}{
    Set to \code{"none"} to return no credibility intervals. Set to
    \code{"asymp"} to return approximate 95\% CIs (obtained via the
    delta method). Set to \code{"exact"} to obtain exact
    CIs via Monte Carlo simulation. Options \code{"asymp"}
    and especially \code{"exact"} can increase the
    computation time substantially.
}
  \item{niter}{
    Number of Monte Carlo iterations. Only used when \code{citype=="exact"}.
}
  \item{burnin}{
    Number of burnin Monte Carlo iterations. Only used when \code{citype=="exact"}.
  }
  \item{keep.pbam}{Set to \code{TRUE} to save processed bam object, as
    returned by \code{procBam}. This object can require substantial memory
    during execution and disk storage upon saving and is not needed
    for a default analysis.
}
  \item{keep.multihits}{Set to \code{FALSE} to discard reads aligned to
  multiple positions.
}
  \item{chroms}{Manually set chromosomes to be processed. By default
  only main chromosomes are considered (except 'chrM')
}
}

\details{
  The function executes the functions \code{procBam}, \code{getDistrs} and \code{pathCounts}
  in parallel for each chromosome, but is much more efficient in cpu
  speed and memory usage than running these functions separately.
  Data from multiple samples are then combined using \code{mergeExp}.
  Note that further normalization (e.g. \code{quantileNorm})
  may be needed preliminary to actual data analysis.

  When \code{rpkm} is false the function returns the estimated
  proportion of reads arising from each isoform within a gene island.
  casper groups two or more genes into a gene island whenever these
  genes share an exon (or part of an exon). Because exons are shared,
  isoform quantification must be done simultaneously for all those
  genes.

  That is, the output from \code{wrapKnown} when \code{rpkm} is \code{FALSE} are proportions that add up
  to 1 within each island. If you would like to re-normalize these
  expressions so that they add up to 1 within each gene, see the help
  for function \code{relexprByGene}.

  One last remark: casper returns the estimated proportion of reads
  generated by each isoform, which is not the same as relative
  isoform expressions. Longer isoforms tend to
  produce more reads than shorter isoforms. This is easily accounted for
  by dividing relative expressions by isoform length, see \code{relexprByGene}.
}
\value{
\item{distr}{Object of class \code{readDistrs}}
\item{pbam}{List of objects of class \code{procBam} with one element
  per chromosome}
\item{pc}{Object of class \code{pathCounts}}
\item{exp}{Object of class \code{ExpressionSet}}
}
\author{
  Camille Stephan-Otto Attolini, David Rossell
}
\references{
Rossell D, Stephan-Otto Attolini C, Kroiss M, Stocker A. Quantifying Alternative Splicing from Paired-End
RNA-sequencing data. Annals of Applied Statistics, 8(1):309-330.
}
\seealso{
\code{procGenome}, \code{relexprByGene}, \code{quantileNorm}
}
\examples{
## genDB<-makeTranscriptDbFromUCSC(genome="hg19", tablename="refGene")
## hg19DB <- procGenome(genDB, "hg19")
##  bamFile="/path_to_bam/sorted.bam"
## ans <- wrapKnown(bamFile=bamFile, mc.cores.int=4, mc.cores=3, genomeDB=hg19DB, readLength=101)
##  names(ans)
##  head(exprs(ans\$exp))
}
