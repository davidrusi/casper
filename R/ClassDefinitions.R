require(methods)



setClass("procBam", representation(pbam = "GRanges", junx="GRanges", stranded = "logical", plus="GRanges", minus="GRanges", pjunx="GRanges", mjunx="GRanges"))

setMethod("show", signature(object="procBam"), function(object) {
  cat("procBam object created from",ifelse(object@stranded,"stranded","non-stranded"),"reads\n")
  if(object@stranded) {
    cat("Contains",length(object@plus),"ranges in the positive strand corresponding to",length(unique(values(object@plus)$id)),"unique read pairs\n")
    cat("and",length(object@minus),"ranges in the negative strand corresponding to",length(unique(values(object@minus)$id)),"unique read pairs\n")
  } else cat("Contains",length(object@pbam),"ranges corresponding to",length(unique(values(object@pbam)$id)),"unique read pairs\n")
}
          )


setMethod("getReads", signature(x="procBam"), function(x) x@pbam)


###############################################################

setClass("readDistrs", representation(lenDis = "array", stDis = "function"))

setMethod("show", signature(object="readDistrs"), function(object) {
  cat("readDistrs object\n\n")
  cat("Insert size distribution (only first few shown)\n")
  show(head(object@lenDis,n=6))
  cat("...\n")
  cat("Read start cumulative distribution function in slot stDis\n")
}
)

setMethod("plot", signature(x="readDistrs"), function(x, y, ...) {
  if (missing(y)) stop("Specify type of plot: 'fragLength' or 'readSt'")
  if (y=='fragLength') {
    n <- as.numeric(names(x@lenDis))
    ylim <- c(0,max(x@lenDis))
    x2plot <- double(max(n)-min(n)+1); names(x2plot) <- min(n):max(n)
    x2plot[names(x@lenDis)] <- x@lenDis
    plot(as.numeric(names(x2plot)),x2plot,type='l',xlab='Fragment length',ylab='Counts',ylim=ylim)
  } else if (y=='readSt') {
    s <- seq(0,1,by=0.02)
    probs <- diff(x@stDis(s))
    plot(NA,NA,xlim=c(0,1),ylim=c(0,max(probs[!is.na(probs)])),xlab='Read start (relative to transcript length)',ylab='Density')
    segments(s[-length(s)],probs,s[-1])
    segments(s,c(0,probs),s,c(probs,0))
  } else {
    stop("Second argument must be either 'fragLength' or 'readSt'")
  }
}
)


###############################################################


setClass("pathCounts", representation(counts = "list", denovo = "logical", stranded="logical"))

valid_pathCounts <- function(object) {
  msg <- NULL
#validity checks go here
  if(object@stranded & length(object@counts)<2) msg <- "Stranded pathCounts must contain plus and minus counts"
  if(!(is.null(msg))) { TRUE } else { msg }
}

setValidity("pathCounts", valid_pathCounts)
setMethod("show", signature(object="pathCounts"), function(object) {
  if(object@denovo) {
    if(object@stranded) {
      cat("Stranded denovo pathCounts object with",length(object@counts[[1]])+length(object@counts[[2]]),"islands and", sum(!unlist(lapply(object@counts[['plus']], is.null))),"non zero positive islands and ", sum(!unlist(lapply(object@counts[['minus']], is.null))),"non zero negative islands\n")
          } else cat("Non-stranded denovo pathCounts object with",length(object@counts[[1]]),"islands and",sum(!unlist(lapply(object@counts[[1]], is.null))),"non zero islands.\n") 
  } else {
    if(object@stranded) {
      cat("Stranded known pathCounts object with",length(object@counts[[1]])+length(object@counts[[2]]),"islands and", sum(!unlist(lapply(object@counts[['plus']], is.null))),"non zero positive islands and ", sum(!unlist(lapply(object@counts[['minus']], is.null))),"non zero negative islands\n")
    } else cat("Non-stranded known pathCounts object with",length(object@counts[[1]]),"islands and", sum(!unlist(lapply(object@counts[[1]], is.null))),"non zero islands.\n")
  }
})



###############################################################


setClass("annotatedGenome", representation(islands = "GRangesList", transcripts = "list", exon2island = "data.frame", exonsNI="GRanges",  aliases="data.frame", genomeVersion="character", dateCreated="Date", denovo="logical"))
valid_annotatedGenome <- function(object) {
  msg <- NULL
  if (!(all(c('seqnames','start','end','width','id','island') %in% names(object@exon2island)))) msg <- "Incorrect column names in 'exon2island'"
  if(!(all(c('tx_id','tx_name','gene_id','exid','tx','island_id') %in% names(object@aliases)))) msg <- "Incorrect column names in 'aliases'"
  if (is.null(msg)) { TRUE } else { msg }
}


setValidity("annotatedGenome", valid_annotatedGenome)
setMethod("show", signature(object="annotatedGenome"), function(object) {
  if(object@denovo) {
    cat("Denovo annotatedGenome object with",length(object@islands),"gene islands,", length(unique(unlist(lapply(object@transcripts, names)))),"transcripts and",nrow(object@exon2island),"exons.\n")
  } else cat("Known annotatedGenome object with",length(object@islands),"gene islands,", length(unique(unlist(lapply(object@transcripts, names)))),"transcripts and",nrow(object@exon2island),"exons.\n")
  cat("Genome version:",object@genomeVersion,"\n")
  cat("Date created:", as.character(object@dateCreated),"\n")
}
          )
