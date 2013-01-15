\name{K562.r1l1}
\alias{K562.r1l1}
\docType{data}
\title{
Toy RNA-seq data from RGASP project.
}
\description{
The paired-end RNA-seq data is from the RGASP project 
sample K562_2x75 (replicate 1, lane 1)
and was obtained at
ftp://ftp.sanger.ac.uk/pub/gencode/rgasp/RGASP1/inputdata/human_fastq.
Reads were aligned against hg19 with tophat 2.0.2 and bowtie 0.12.5, setting the insert size at -r 200,
and imported into R using \code{scanBam} from package \code{Rsamtools}.
For illustration purposes, we selected reads mapping to a few genes only
(namely, the genes that were also selected for the toy genome annotation in \code{data(hg19DB)}.
}
\usage{data(K562.r1l1)}
\format{
 A list indicating read id, chromosome, start and end locations and the position of the pair, as returned by \code{scanBam}.
}
\source{
ftp://ftp.sanger.ac.uk/pub/gencode/rgasp/RGASP1/inputdata/human_fastq
}
\references{
C Trapnell, L Pachter, SL Salzberg. TopHat: discovering splice junctions with RNA-Seq. Bioinformatics, 2009, 25, 1105-1111. doi=10.1093/bioinformatics/btp120.

B Langmead, C Trapnell, M Pop, SL Salzberg. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology, 2009, 10:R25.
}
\examples{
data(K562.r1l1)
names(K562.r1l1)
}
\keyword{datasets}
