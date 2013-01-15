## METHODS

setMethod("getNreads",signature(pc='pathCounts'),function(pc) {
  if(pc@stranded) {
    ans <- c(sapply(pc@counts[['plus']], sum), sapply(pc@counts[['minus']], sum))
  } else {
    ans <- sapply(pc@counts[[1]], sum)
  }
  ans
}
          )

procPaths <- function(reads, DB, mc.cores){
    cat("Finding overlaps between reads and exons\n")
    over<-findOverlaps(reads, DB@exonsNI)    
    readid<-as.character(values(reads)$id)[queryHits(over)]
    readside<-values(reads)$rid[queryHits(over)]
    exid <- names(DB@exonsNI)[subjectHits(over)]
    exst<-start(DB@exonsNI)[subjectHits(over)]
    cat("Counting paths\n")
    pCounts<-.Call("pathCounts", readid, readside, exst, exid)
    pCounts<-lapply(pCounts[1:2], function(x) x[1:pCounts[[3]]])
    names(pCounts[[2]])<-pCounts[[1]]
    pCounts<-pCounts[[2]]
    pCounts <- pCounts[grepl("-[0-9]", names(pCounts))]
    
    sel <- strsplit(names(pCounts), split='-|\\.')
    sel1 <- lapply(sel, "[", 2)
    sel1 <- unlist(sel1)
    
    nislEx <- elementLengths(DB@islands)
    nislEx <- rep(names(DB@islands), unlist(nislEx))
    islEx <- names(unlist(DB@islands))
    islEx <- sub(".*\\.","", islEx)
    names(islEx) <- nislEx
    isl <- match(sel1, islEx)
    isl <- names(islEx)[isl]
    splCounts <- split(pCounts, isl)
    splCounts <- lapply(splCounts, function(x) x[grepl("-", names(x))])
    
    if(DB@denovo){
      sel <- sapply(sel, "[", -1)
      tmp <- split(sel, isl)
      if(mc.cores>1) {
        require(multicore)
        tmp1 <- multicore::mclapply(names(tmp), function(x){
          n <- sapply(tmp[[x]], length)
          nn <- unlist(tmp[[x]])
          names(nn) <- rep(names(splCounts[[x]]), n)
          nnn <- nn %in% names(DB@islands[[x]])
          nnnn <- tapply(nnn, names(nn), all)
          splCounts[[x]][names(splCounts[[x]]) %in% names(nnnn)[nnnn]]
        }, mc.cores=mc.cores)
      } else {
        tmp1 <- lapply(names(tmp), function(x){
          n <- sapply(tmp[[x]], length)
          nn <- unlist(tmp[[x]])
          names(nn) <- rep(names(splCounts[[x]]), n)
          nnn <- nn %in% names(DB@islands[[x]])
          nnnn <- tapply(nnn, names(nn), all)
          splCounts[[x]][names(splCounts[[x]]) %in% names(nnnn)[nnnn]]
        })
      }
      names(tmp1) <- names(tmp)
      splCounts <- tmp1
    }
    
    ans <- vector(length(DB@islands), mode='list')
    names(ans) <- names(DB@islands)
    ans[names(splCounts)] <- splCounts
    ans
  }

pathCounts<-function(reads, DB, mc.cores=1) {
  if (class(reads) != 'procBam') stop('reads must be an object of class procBam')
  if(!reads@stranded) {
    counts <- procPaths(reads@pbam, DB, mc.cores)
    ans <- new("pathCounts", counts=list(counts), denovo=DB@denovo, stranded=reads@stranded)
  }
  else {
    plusDB <- genomeBystrand(DB, strand="+")
    plus <- procPaths(reads@plus, plusDB, mc.cores)
    minusDB <- genomeBystrand(DB, strand="-")
    minus <- procPaths(reads@minus, minusDB, mc.cores)
    ans <- new("pathCounts", counts=list(plus=plus, minus=minus), denovo=DB@denovo, stranded=reads@stranded)
  }
  ans
}
