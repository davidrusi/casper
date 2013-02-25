require(methods)

setGeneric("subsetPbam", function(pbam, strand) standardGeneric("subsetPbam"))
setMethod("subsetPbam", signature(pbam="procBam", strand="character"),
          function(pbam, strand) {
            if(strand=="+") {
              pbam@pbam <- pbam@plus
              pbam@junx <- pbam@pjunx
            } else {
              pbam@pbam <- pbam@minus
              pbam@junx <- pbam@pminus
            }
            pbam
          }
          )


buildRD<-function(reads){
  sel <- reads$end<reads$start
  st <- en <- integer(length(reads$start))
  st[sel] <- reads$end[sel]; st[!sel] <- reads$start[!sel]
  en[!sel] <- reads$end[!sel]; en[sel] <- reads$start[sel]
  if(sum(grepl("chr", unique(reads$chrom)))>0) {
    reads<-GRanges(ranges=IRanges(start=st, end=en), seqnames=reads$chrom, id=reads$key, flag=reads$flag, rid=reads$rid, strand=reads$rstrand, XS=reads$strand)
  } else {
    reads<-GRanges(ranges=IRanges(start=st, end=en), seqnames=paste("chr", reads$chrom, sep=""), id=reads$key, flag=reads$flag, rid=reads$rid, strand=reads$rstrand, XS=reads$strand)
  }
  return(reads)
}

uniquifyQname<-function(bam, seed=1){   	
  qname <- vector(mode='character', length=length(bam$qname))
  dups <- .Call("uniqQname", bam$qname, length(bam$qname), bam$pos, bam$mpos, qname)
  bam$qname <- dups[[1]]
  #if(length(dups[[2]])>0){
  #  probPos<-bam$qname %in% dups[[2]]
  #  probID<-paste(bam$qname[probPos], bam$pos[probPos])
  #  set.seed(seed)
  #  sel<-unlist(tapply(1:sum(probPos), probID, function(x) x[sample(1:length(x), 1)]))
  #  fixed<-lapply(lapply(bam, "[", probPos), "[", sel)
  #  bam1<-lapply(bam, function(x) x[!probPos])    
  #  res<-lapply(names(bam), function(x) c(bam1[[x]], fixed[[x]]))
  #  names(res)<-names(bam)
  #} else {
  #  res<-bam
  #}
  probPos <- bam$qname %in% dups[[2]]
  bam <- lapply(bam, function(x) x[!probPos])    
  qdup <- bam$qname[bam$qname %in% bam$qname[duplicated(bam$qname)]] 
  uni <- bam$qname %in% qdup
  res <- lapply(bam, "[", uni)
  names(res)<-names(bam)
  res
}

nbReads <- function(bam0) {
  tab <- table(bam0$cigar)
  count <- sapply(gregexpr("M",names(tab)),length)
  tjunx <- sum(tab*(count-1))
  treads <- sum(tab*count)
  c(tjunx=tjunx, treads=treads)
}

procB <- function(bam, strnd, seed=1, verbose=FALSE){

  lev=NULL
  if(class(bam$strand)=='factor') lev <- levels(bam$strand)
  bam$rname<-as.character(bam$rname)
  bam<-uniquifyQname(bam, seed)
  if(verbose) cat("Calculating total number of reads...\n")
  
  nreads <- nbReads(bam)
  njunx <- nreads['tjunx']
  nreads <- nreads['treads']+10

  if(verbose) cat("done.\nProcessing cigars and building read's object...\n")
  len=vector(mode="integer", length=nreads)
  flag=vector(mode="integer", length=nreads)
  strs=vector(mode="integer", length=nreads)
  rid=vector(mode="integer", length=nreads)
  key=vector(mode="character", length=nreads)
  chrom=vector(mode="character", length=nreads)
  strand=vector(mode="integer", length=nreads)
  jchrom=vector(mode="character", length=njunx)
  jstrs=vector(mode="integer", length=njunx)
  jlen=vector(mode="integer", length=njunx)
  data<-.Call("procBam", bam$qname, bam$flag,  bam$rname, bam$pos, bam$cigar, as.numeric(bam$strand), length(bam$pos), nreads, njunx, len, strs, flag, key, chrom, rid, strand, jchrom, jstrs, jlen)
  names(data) <- c('end', 'start', 'flag', 'key', 'chrom', 'rid', 'rstrand', 'jchrom', 'jstart', 'jend')
  junx <- GRanges(ranges=IRanges(start=data$jstart[data$jstart!=0], end=data$jend[data$jstart!=0]), seqnames=data$jchrom[data$jstart!=0], XS=rep(strnd, sum(data$jstart!=0)))
  if(!is.null(lev)) data[[7]] <- lev[data[[7]]]
  sel<-data[[1]]!=0 & !is.na(data[[3]])
  data<-lapply(data, "[", sel)
  data[[length(data)+1]] <- rep(strnd, length(data[[1]]))
  names(data)[length(data)] <- 'strand'
  if(verbose) cat("done.\n")
  ans<-buildRD(data)
  id <- values(ans)$id
  idnum <- rep(0,length(id))
  idnum[which(id[-length(id)] != id[-1])+1] <- 1
  idnum <- cumsum(idnum)
  values(ans)$id <- idnum
  list(pbam=ans, junx=junx)
}

setMethod("procBam", signature(bam='list',stranded='missing',seed='missing', verbose='missing') ,
          function(bam, stranded, seed, verbose) procBam(bam=bam, stranded=FALSE, seed=as.integer(1), verbose=FALSE)
          )

setMethod("procBam", signature(bam='list',stranded='missing',seed='missing', verbose='logical') ,
          function(bam, stranded, seed, verbose) procBam(bam=bam, stranded=FALSE, seed=as.integer(1), verbose=verbose)
          )

setMethod("procBam", signature(bam='list',stranded='missing',seed='integer', verbose='missing') ,
          function(bam, stranded, seed, verbose) procBam(bam=bam, stranded=FALSE, seed=seed, verbose=FALSE)
          )

setMethod("procBam", signature(bam='list',stranded='logical',seed='missing', verbose='logical') ,
          function(bam, stranded, seed, verbose) procBam(bam=bam, stranded=stranded, seed=as.integer(1), verbose=verbose)
          )

setMethod("procBam", signature(bam='list',stranded='logical',seed='integer', verbose='logical') ,
          function(bam, stranded, seed, verbose) procBam(bam=bam, stranded=stranded, seed=seed, verbose=verbose)
)

setMethod("procBam", signature(bam='list',stranded='logical',seed='integer', verbose='logical') ,
          function(bam, stranded, seed, verbose) {
            byList <- ifelse(class(bam[[1]])=='list', TRUE, FALSE)
            if(!byList) {
              ans <- procBamF(bam=bam, stranded=stranded, seed=seed, verbose=verbose)
            } else {
              ans <- lapply(bam, function(x) procBamF(bam=x, stranded=stranded, seed=seed, verbose=verbose))
              ans <- mergePbam(ans, stranded)
              lens <- seqnames(ans@pbam)@lengths
              fix <- rep(c(0, (values(ans@pbam)$id[cumsum(lens)]+1)), c(seqnames(ans@pbam)@lengths, 0))
              values(ans@pbam)$id <- values(ans@pbam)$id+fix
            }
            ans
          }
            )

mergeGRanges <- function(gr) {
   values <- do.call('rbind',unname(lapply(gr, values)))
   colnames(values) <- sub("values.", "", colnames(values))
   seqn <- unlist(RleList(unname(lapply(gr, seqnames))))
   ranges <- unlist(IRangesList(lapply(gr, ranges)))
   ans <- GRanges(ranges=ranges, seqnames=seqn)
   values(ans) <- values
   unname(ans)
 }

mergePbam <- function(ans, stranded){
  if(stranded){
    plus <- mergeGRanges(lapply(ans, function(x) slot(x, 'plus')))
    minus <- mergeGRanges(lapply(ans, function(x) slot(x, 'minus')))
    pjunx <- mergeGRanges(lapply(ans, function(x) slot(x, 'pjunx')))
    mjunx <- mergeGRanges(lapply(ans, function(x) slot(x, 'mjunx')))
    new("procBam", pbam=NULL, plus=plus, minus=minus,junx=NULL, pjunx=pjunx, mjunx=mjunx, stranded=ans$stranded)
  } else{
    pbam <- mergeGRanges(lapply(ans, function(x) slot(x, 'pbam')))
    junx <- mergeGRanges(lapply(ans, function(x) slot(x, 'junx')))
    ans <- new("procBam",pbam=pbam, junx=junx,stranded=stranded)
  }
  ans
}

    
procBamF<-function(bam, stranded=FALSE, seed=1,  verbose=FALSE){
  require(IRanges)
  if(stranded) {
    minus <- bam$tag$XS=='-'
    mjunx <- GRanges(IRanges(0,0), space=NULL)
    if(sum(minus)>0){
      minus <- lapply(bam, '[', minus)
      minus <- procB(minus, "-", seed=seed, verbose=verbose)
      mjunx <- minus$junx
      minus <- minus$pbam
    } else minus <- GRanges(IRanges(0,0), space=NULL)
    plus <- bam$tag$XS=='+'
    mplus <- GRanges(IRanges(0,0), space=NULL)
    if(sum(plus)>0){
      plus <- lapply(bam, '[', plus)
      plus <- procB(plus, "+", seed=seed, verbose=verbose)
      pjunx <- plus$junx
      plus <- plus$pbam
    } else plus <- GRanges(IRanges(0,0), space=NULL)
    ans <- list(plus=plus, pjunx=pjunx, minus=minus, mjunx=mjunx, stranded=TRUE)
  } else {
    pbam <- procB(bam, "*", seed=seed, verbose=verbose)
    ans <- list(pbam=pbam$pbam, junx=pbam$junx, stranded=FALSE)
  }
  if(!ans$stranded) ans <- new("procBam",pbam=ans$pbam, junx=ans$junx,stranded=ans$stranded)
  else ans <- new("procBam",pbam=NULL, plus=ans$plus, minus=ans$minus,junx=NULL, pjunx=ans$pjunx, mjunx=ans$mjunx, stranded=ans$stranded)
}
