\name{denovoGeneExpr-class}
\Rdversion{1.1}
\docType{class}
\alias{denovoGeneExpr-class}
\alias{show,denovoGeneExpr-method}
\alias{[,denovoGeneExpr,ANY,ANY,ANY-method}
\alias{[[,denovoGeneExpr,ANY,ANY-method}
\alias{[[,denovoGeneExpr,ANY,ANY,ANY-method}
\alias{posprob}
\alias{variants}
\alias{variants<-}
\alias{posprob,denovoGeneExpr-method}
\alias{variants,denovoGeneExpr-method}
\alias{variants<-,denovoGeneExpr-method}
\alias{variants,denovoGenomeExpr-method}
\alias{names,denovoGenomeExpr-method}

\title{Class "denovoGeneExpr"}
\description{
  \code{denovoGeneExpr} stores inferred expression for de novo splicing
  variants for a single gene.
  \code{denovoGenomeExpr} stores the information for several genes
  (typically, the whole genome).
}
\section{Objects from the Class}{
  Objects are returned by \code{calcDenovo}. When running
  \code{calcDenovo} on multiple genes results are returned in a
  \code{denovoGenomeExpr} object. Results for a single gene can be
  retrieved using the [[ operator as usual, which returns a
  \code{denovoGeneExpr} object.
}
\section{Slots}{
  \describe{
    \item{posprob}{\code{data.frame} containing the posterior
      probability of each model}
    \item{expression}{\code{data.frame} with the estimated expression of
    each variant under each model}
    \item{variants}{\code{matrix} indicating the exons contained in
      each variant.}
    \item{integralSum}{Sum of the log(integrated likelihood) + log(model
    prior probability) across all considered models.}
    \item{npathDeleted}{Number of paths that had 0 probability under all
    considered variants and had to be excluded for model fitting
    purposes.}
    \item{priorq}{Input parameter to \code{calcDenovo}}
    \item{txLength}{Length of transcripts in bp (including new isoforms
      found by casper)}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "denovoGeneExpr")}: Displays
      general information about the object. }
    \item{names}{Show names (island ids)}
    \item{"["}{Selects a subset of genes}
    \item{"[["}{Selects a single gene}
    \item{posprob}{Accesses the posterior probabilities of each model
      (slot posprob)}
    \item{variants}{Accesses the variant names and their respective exons}
    \item{variants<-}{Replaces the value of the slot variants (can be
      useful for renaming variants, for instance)}
	 }
}
\author{
David Rossell
}
\seealso{
  \code{\link{calcDenovo}} to create objects from the class.
  \code{\link{denovoExpr}} to obtain expression estimates from
  \code{denovoGenomeExpr} objects.
}
\examples{
showClass("denovoGeneExpr")
}
\keyword{classes}
