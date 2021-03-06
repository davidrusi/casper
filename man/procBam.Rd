\name{procBam}
\alias{procBam}
\alias{procBam-method}
\alias{procBam,ANY-method}
\alias{procBam,list,logical,integer,logical,character-method}
\alias{procBam,list,logical,integer,logical,missing-method}
\alias{procBam,list,logical,missing,logical,missing-method}
\alias{procBam,list,missing,missing,missing,missing-method}
\alias{procBam,list,missing,integer,logical,missing-method}
\alias{procBam,list,missing,integer,missing,missing-method}
\alias{procBam,list,missing,numeric,missing,missing-method}
\alias{procBam,list,logical,integer,missing,missing-method}
\alias{procBam,list,logical,missing,missing,missing-method}
\alias{procBam,list,missing,missing,logical,missing-method}
\docType{methods}
\title{
  Process BAM object
}
\description{
  Process paired-end data stored in BAM object generated by scanBam.
  Outputs GRanges objects for reads and junctions.
  }
\usage{
procBam(bam, stranded=FALSE, seed=as.integer(1), verbose=FALSE, rname='null',
keep.junx=FALSE, keep.flag=FALSE, ispaired=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam}{
    BAM object generated by \code{scanBam}
  }
  \item{stranded}{
    Set to \code{TRUE} to indicate that the RNA-seq experiment preserved
    the strand information.
  }
  \item{seed}{
    Seed for random number generator
  }
  \item{verbose}{Set to \code{TRUE} to print progress information.
  }
  \item{rname}{
   Chromosome to process be combined with the \code{which} argument in the scanBam function
  }
  \item{keep.junx}{
   Option to store junction information. Only useful for finding denovo
  exons and transcripts.
  }
  \item{keep.flag}{
      Option to store aligment flag information.
  }
  \item{ispaired}{Set to \code{TRUE} is reads are paired.}
  \item{...}{Other arguments}
}


\section{Methods}{
  \describe{
    \item{\code{signature(bam='list',stranded='logical',seed='integer',verbose='logical',
  rname='character',keep.junx='logical',keep.flag='logical')}}{
      Process paired-end data stored in BAM object generated by scanBam.
      Outputs GRanges objects for reads and (optionally) junctions.
    }
  }
}
\details{
In case of multihits with same start position for both reads but
  different insertions/deletions patterns only one alignment is chosen
  at random.
}
\value{
  An object of class \code{procBam} containing
  reads with both ends
  correctly aligned and split according to the corresponding CIGAR.
  Unique identifiers by fragment are stored. Junctions spanned by reads
  are also stored in GRanges object if the argument \'keep.junx\' is set
  to TRUE.}
\author{
  Camille Stephan-Otto Attolini
}
\seealso{
  \code{scanBam} from package \code{Rsamtools},
  \code{help("procBam-class")}, \code{getReads}.
}
\examples{
##See example in calcExp
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

