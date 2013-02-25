\name{transcripts}
\alias{transcripts}
\alias{transcripts,character,missing,annotatedGenome-method}
\alias{transcripts,missing,character,annotatedGenome-method}
\docType{methods}
\title{
Returns all transcripts associated to a given island or
entrez id stored in an annotatedGenome object.
}
\description{
annotatedGenome objects store information regarding genes and
transcripts. When there's an overlap in exons between several genes,
these genes are grouped into gene islands.

The function transcripts retrieves all stored transcripts for a given
gene or island.
}
\usage{
transcripts(entrezid, islandid, genomeDB)
}
\arguments{
 \item{entrezid}{Character indicating single Entrez identifier. Can be
   left missing and specify islandid instead.}
 \item{islandid}{Character indicating island identifier. Can be left
   missing and specify entrezid instead}
 \item{genomeDB}{Object of class \code{annotatedGenome}}
}
\section{Methods}{
\describe{
  \item{\code{signature(x = "readDistrs")}}{
    \code{x} is an object of type \code{readDistrs}, as returned by
    \code{getDistrs}. The plot function allows to visualize the fragment
    length and read start distributions.
  }
}}
\value{
\code{IRangesList} where each element in the list corresponds to a
different transcript.
}
\examples{
data(hg19DB)
tx <- transcripts(entrezid="27",genomeDB=hg19DB)
tx
}
\keyword{manip}
