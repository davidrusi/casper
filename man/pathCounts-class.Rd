\name{pathCounts-class}
\Rdversion{1.1}
\docType{class}
\alias{pathCounts-class}
\alias{show,pathCounts-method}

\title{Class "pathCounts"}
\description{
  Stores exon path counts.
}
\section{Objects from the Class}{
  Objects are created with a call to \code{pathCounts}.
}
\section{Slots}{
  \describe{
    \item{counts}{List with one element per gene island. For each
      island, it contains a named vector with exon path counts. The names indicate
      the visited exons.
      
      For instance, consider that for gene '1' with 2 exons we observe 10
      reads in which the left end falls completely in exon 1 and the right
      end in exon 2. Suppose that for 5 reads the left end bridges exons 1-2
      and the right end falls in exon 2. Then \code{pc[['1']]} would contain
      \code{c(10,5)} and \code{names(pc[['1']])} would contain \code{c(".1-2.",".1.2-2.")}
    }
    \item{denovo}{Logical variable. \code{FALSE} indicates that the
      counts correspond to a known genome (i.e. created with \code{procGenome}), and \code{TRUE} to a de novo
      annotated genome (i.e. created with \code{createDenovoGenome}).
    }
    \item{stranded}{Logical variable. \code{TRUE} indicates that the
      path counts were obtained from and RNA-seq experiment where strand
      information was preserved.
    }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "pathCounts")}: Displays general information about the object. }
	 }
}
\author{
Camille Stephan-Otto Attolini
}
\examples{
showClass("pathCounts")
}
\keyword{classes}
