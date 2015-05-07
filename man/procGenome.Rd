\name{procGenome}
\alias{procGenome}
\alias{procGenome,TxDb-method}
\alias{procGenome,GRanges-method}
\alias{procGenome,TxDb,ANY-method}
\alias{procGenome,GRanges,ANY-method}
\alias{createDenovoGenome}
\title{
  Create an annotatedGenome object that stores information about genes and transcripts
}
\description{
  \code{procGenome} processes annotations for a given transcriptome,
  either from a \code{TxDb} object created by \code{GenomicFeatures} package
  (e.g. from UCSC) or from a user-provided \code{GRanges} object (e.g. by
  importing a gtf file).
  
  \code{createDenovoGenome} creates a de novo annotated genome by
  combining UCSC annotations and observed RNA-seq data.
}
\usage{
procGenome(genDB, genome, mc.cores=1, verbose=TRUE)

createDenovoGenome(reads, DB, minLinks=2,
maxLinkDist=1e+05, maxDist=1000, minConn=2, minJunx=3, minLen=12, mc.cores=1)
}
\arguments{
  \item{genDB}{Either a \code{TxDb} object with annotations (e.g. from
    UCSC or a gtf file or a \code{GRanges} object as returned by
    \code{import} from \code{rtracklayer} package). See details.}
  \item{genome}{Character indicating genome version (e.g. "hg19",
    "dm3")}
  \item{mc.cores}{Number of cores to use in parallel processing
    (multicore package required)}
  \item{verbose}{Set to \code{TRUE} to print progress information}
  \item{DB}{\code{annotatedGenome} object, as returned by
    \code{procGenome}}
  \item{minLinks}{Minimum number of reads joining two
    exons to merge their corresponding genes}
  \item{maxLinkDist}{Maximum distance between two exons to merge
    their correspondin genes. A value of \code{0} disables this option.}
  \item{maxDist}{Maximum distance between two exons with reads joining
    them to merge their corresponding genes.}
  \item{minConn}{Minimum number of fragments connecting a new exon to an
    annotated one to add to denovo genome.}
  \item{minJunx}{Minimum number of junctions needed to redefine an annotated
    exon's end or start.}
  \item{minLen}{Minimum length of a junction to consider as a putative intron.}
  \item{reads}{Processed reads stored in a \code{RangedData}, as
    returned by \code{procBam}}
}
\section{Methods}{
\describe{
  \item{\code{signature(genDB = "transcriptDb")}}{
    \code{genDB} is usually obtained with a call to
    \code{makeTxDbFromUCSC} (package \code{GenomicFeatures}),
    e.g. genDB<-makeTxDbFromUCSC(genome="hg19", tablename="refGene")
  }
  \item{\code{signature(genDB = "GRanges")}}{
    \code{genDB} stores information about all transcripts and their
    respective exons. Chromosome, start, end and strand are stored as
    usual in GRanges objects. \code{genDB} must have a column named
    \code{"type"} taking the value \code{"transcript"} for rows
    corresponding to transcript and \code{"exon"} for rows corresponding
    to exons. It must also store transcript and gene ids. For instance, Cufflinks RABT
    module creates a gtf file with information formatted in this manner
    for known and de novo predicted isoforms.
  }
}}

\details{
  These functions create the annotation objects that are needed for subsequent
  functions. Typically these objects are created only once for a set of
  samples.

  If interested in quantifying expression for known transcripts
  only, one would typically use \code{procGenome} with a
  \code{TxDb} from the usual Bioconductor annotations,
  e.g. genDB<-makeTxDbFromUCSC(genome="hg19",tablename="refGene"),
  or imported from a gtf file
  e.g. genDB<-makeTxDbFromGFF('transcripts.gft',format='gtf').
  \code{GRanges} object (e.g. genDB <- import('transcripts.gtf')).
  Package GenomicFeatures contains more info about how to create
  \code{TxDb} objects.
  Alternatively, one can provide annotations as a \code{GRanges} object
  whith is returned when importing a gtf file with
  function \code{import} (package \code{rtracklayer}).
  
  The output from \code{procGenome} can be used in combination with
  \code{wrapKnown}, which quantifies expression for a set of known transcripts,
  or \code{wrapDenovo}, which uses Bayesian model selection methods to
  assess which transcripts are truly expressed.
  When using \code{wrapDenovo}, you should create a single \code{annotatedGenome}
  object that combines information from all samples
  (e.g. from a gtf file produced by running your favorite isoform
  prediction software jointly on all samples),
  as this increases the power to detect new exons and isoforms.
}
\value{
  Object of class \code{annotatedGenome}.
}
\author{
  Camille Stephan-Otto Attolini
}
\seealso{
  See \code{annotatedGenome-class} for a description of the class.
  To obtain the island id for an entrez gene id or a transcript id see help(getIsland).
  To extract the transcripts for a entrez or island id see
  help(transcripts).
  See \code{splitGenomeByLength} for splitting an \code{annotatedGenome}
  according to gene length.
}
\examples{
##Known transcripts from UCSC
## genDB<-makeTxDbFromUCSC(genome="hg19", tablename="refGene")
## hg19DB <- procGenome(genDB, "hg19")

##Transcripts inferred by Cufflinks-RABT module
##genDB.Cuff <- import('transcripts.gtf')
##hg19DB.Cuff <- procGenome(genDB.Cuff, genome='hg19')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~annotation }
