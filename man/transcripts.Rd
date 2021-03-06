\name{transcripts}
\alias{transcripts}
\alias{transcripts,annotatedGenome,missing,missing-method}
\alias{transcripts,annotatedGenome,character,missing-method}
\alias{transcripts,annotatedGenome,missing,character-method}
\alias{matchTranscripts}

\docType{methods}
\title{
  Extracts transcript information (exon start
  and ends) from an \code{annotatedGenome} object,
  either for all transcripts or only those corresponding
  to a given island or transcript.
}
\description{
annotatedGenome objects store information regarding genes and
transcripts. When there's an overlap in exons between several genes,
these genes are grouped into gene islands.

\code{transcripts} retrieves all stored transcripts for a given
transcript or island.

\code{matchTranscripts} finds transcripts in queryDB matching a
transcript in subjectDB. The best match for each transcript
in subjectDB is returned, unless difference in bp is >maxbp
}
\usage{
transcripts(genomeDB, txid, islandid)

matchTranscripts(queryDB, subjectDB, maxbp=10) 
}
\arguments{
  \item{genomeDB}{Object of class \code{annotatedGenome}}
  \item{txid}{Character indicating transcript identifier (optional)}
  \item{islandid}{Character indicating island identifier (optional)}

  \item{queryDB}{\code{annotatedGenome} with query transcripts}
  \item{subjectDB}{\code{annotatedGenome} with potentially matching
    transcripts}
  \item{maxbp}{Maximum difference in bp for transcripts to be matched}
}
\section{Methods}{
\describe{
  \item{\code{signature(genomeDB = "annotatedGenome", txid="missing", islandid="missing")}}{
    Return exons for all transcripts in \code{genomeDB}
  }
  \item{\code{signature(genomeDB = "annotatedGenome", txid="character", islandid="missing")}}{
    Return exons for transcript \code{txid}
  }
  \item{\code{signature(genomeDB = "annotatedGenome", txid="missing", islandid="character")}}{
    Return exons for all transcripts in island \code{islandid}
  }
}}
\value{
\code{IRangesList} where each element in the list corresponds to a
different transcript.
}
\seealso{
\code{genePlot} to plot the resulting transcripts
}
\examples{
data(hg19DB)
txs <- transcripts(txid="NM_005158",genomeDB=hg19DB)
txs
}
\keyword{manip}
