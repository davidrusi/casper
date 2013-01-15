require(methods)

findNewExons <- function(pbam, junx, DB, minConn=2, minJunx=3, minLen=12){
#  if(DB@stranded){
#  } else {
## Find new putative exons by reads' islands
  cov <- coverage(pbam)
  isl <- slice(cov, lower=1,  rangesOnly=T)
  isl <- GRanges(ranges=unlist(isl), seqnames=rep(names(isl), sapply(isl, length)))
  cnt <- countOverlaps(isl, pbam)
  isl <- isl[cnt>1]
  isls <- disjoin(c(reduce(c(isl, DB@exonsNI)), DB@exonsNI))

  ## Redefine known and new exons by junctions
  isls <- lapply(unique(as.character(seqnames(DB@exonsNI)@values)), function(x){
    if(sum(seqnames(junx)==x)){
      y <- junx[seqnames(junx) == x]
      y <- y[width(y)>3]
      z <- paste(start(y), end(y), sep='.')
      zz <- table(z)
      y <- y[z %in% names(zz)[zz>=minJunx]]
      values(y) <- NULL
      y <- disjoin(c(isls[seqnames(isls)==x], y))
      ans <- unique(subsetByOverlaps(y, isls[seqnames(isls)==x]))
    } else ans <- disjoin(isls[seqnames(isls)==x])
    ans
  })
  isls <- do.call('c', isls)
  names(isls) <- 1:length(isls)

  ## Find exons with no overlap with known exons
  kex <- subsetByOverlaps(isls, DB@exonsNI)
  nisl <- isls[!(isls %in% DB@exonsNI),]
  over <- findOverlaps(nisl, DB@exonsNI, maxgap=1)
  nex <- nisl[!(1:length(nisl) %in% queryHits(over)),]

  ## Find exons to redefine (grow)
  lex <- nisl[!(names(nisl) %in% names(nex))]
  lex <- lapply(unique(as.character(seqnames(DB@exonsNI)@values)), function(x){
    junx <- junx[seqnames(junx)==x]
    scnt <- table(start(junx[width(junx)>5]))
    ecnt <- table(end(junx[width(junx)>5]))
    tmp <- lex[seqnames(lex)==x]
    cnts <- cbind(ecnt[as.character(start(tmp)-1)], scnt[as.character(end(tmp)+1)])
    cnts[is.na(cnts)] <- 0
    tmp <- tmp[cnts[,1]>=minJunx | cnts[,2]>=minJunx]
    tmp[width(tmp)>minLen]
  })
  lex <- do.call('c', lex)

  ## Find exons connected to annotated ones by reads
  over1 <- findOverlaps(pbam, nex)
  idrea <- values(pbam)$id[queryHits(over1)]
  idnex <- names(nex)[subjectHits(over1)]
  tidnex <- table(idnex)
  conn <- (values(pbam)$id[pbam %in% kex]) %in% idrea
  conn <- (values(pbam)$id[pbam %in% kex])[conn]
  conn <- table(conn)
  conn <- names(conn)[conn>minConn]
  idconn <- idnex[idrea %in% conn]
  nexconn <- nex[unique(idconn)]

  ans <- c(kex, nexconn, lex)
  names(ans) <- 1:length(ans)
  #ans
  return(ans)
}

assignExons2Gene <- function(exons, DB, reads, maxDist=1000, minLinks=2, maxLinkDist=0, stranded=FALSE, mc.cores=1){

  ### All exons were redefined. Tx structure exon names' don't coincide anymore. Remap structure before doing anything else. Fuck 

  over <- findOverlaps(exons, DB@exonsNI)
  exid <- names(exons)[queryHits(over)]
  exkey <- split(exid, names(DB@exonsNI)[subjectHits(over)])

  exlen <- sapply(exkey, length)
  otxs <- unlist(DB@transcripts, recursive=F)
  names(otxs) <- sub("[0-9]+\\.", "", names(otxs))
  newTxs <- unlist(exkey[as.character(sprintf("%d",unlist(otxs)))])
  lenkey <- exlen[as.character(sprintf("%d",unlist(otxs)))]
  lenkey <- tapply(lenkey, rep(names(otxs), sapply(otxs, length)), sum)
  lenkey <- lenkey[names(otxs)]
  newTxs <- unname(as.numeric(newTxs))
  newTxs <- split(newTxs, rep(names(lenkey), lenkey))
  newTxs <- newTxs[names(otxs)]
  newTxs <- split(newTxs, rep(names(DB@transcripts), sapply(DB@transcripts, length)))
  
  exs <- exons[!(exons %in% DB@exonsNI)]
  rea <- subsetByOverlaps(reads, exs)
  
#Select read ids in these exons
  area <- reads[values(reads)$id %in% values(rea)$id]
  over <- findOverlaps(exons, area)
  exs1 <- names(exons[queryHits(over)])
  exs1 <- sub("\\..*","",exs1)
  rea1 <- values(area[subjectHits(over)])$id
  len <- length(unique(rea1))
  if(length(unique(exs1))>1){
    ans <- .Call("joinExons", as.integer(exs1), rea1, len)
    junx <- ans[[1]][ans[[2]]>=minLinks]
    junx <- strsplit(junx, split=".", fixed=T)
    names(junx) <- 1:length(junx)
    nalljunx <- rep(1:length(junx), unlist(lapply(junx, length)))
    alljunx <- unlist(junx)
    tmpex <- as.data.frame(exons)
    rownames(tmpex) <- names(exons)
    tmpex <- tmpex[alljunx,]
    alljunx <- as.integer(alljunx)
    names(alljunx) <- nalljunx
    pos <- (tmpex$end+tmpex$start)/2
    dis <- tapply(pos, names(alljunx), function(x) (max(x)-min(x) < maxLinkDist ))
    alljunx <- alljunx[names(alljunx) %in% as.character(names(dis)[dis])]
    if(sum(!(as.numeric(names(exs)) %in% unique(alljunx)))>0) {
      sing <- as.numeric(names(exs))[!(as.numeric(names(exs)) %in% unique(alljunx))]
      mname <- max(as.numeric(names(alljunx)))
      names(sing) <- (mname+1):(mname+length(sing)) 
      alljunx <- c(alljunx, sing)
    }
  } else {
    alljunx <- unique(exs1)
    names(alljunx) <- '1'
  }

#Build gene islands
  oldexs <- unlist(newTxs, recursive=F)
  txids <- sub("[0-9]+\\.", "", names(oldexs))
  txids <- rep(txids, unlist(lapply(oldexs, length)))
  oldexs <- unlist(oldexs)
  names(oldexs) <- txids
  allexs <- c(oldexs, alljunx)
  islands <- makeIslands(allexs) 
  values(exons)$island <- islands[names(exons)] 
  exon2gene <- as.data.frame(exons)
  exon2gene$id <- names(exons)
  exon2gene$island <- islands[as.character(exon2gene$id)]
  rownames(exon2gene) <- exon2gene$id
  islands <- GRanges(IRanges(exon2gene$start, exon2gene$end), seqnames=exon2gene$seqnames, id=exon2gene$id)
  names(islands) <- exon2gene$id
  exon2island <- as.data.frame(exons)  
  exon2island$id <- names(exons)

  if(any(strand(DB@islands@unlistData)=='+')) {
    strand(islands)[exon2gene$island %in% names(newTxs)] <- '+'
    islands <- split(islands, exon2gene$island[match(exon2gene$id, values(islands)$id)])
  } else {
    strand(islands)[exon2gene$island %in% names(newTxs)] <- '-'
    islands <- split(rev(islands), exon2gene$island[match(exon2gene$id, rev(values(islands)$id))])
  }

  extxs <- unlist(newTxs, recursive=F)
  names(extxs) <- sub("[0-9]+\\.", "", names(extxs))
  sel <- unlist(lapply(extxs, "[", 1))
  sel <- match(as.character(sel), rownames(exon2island))
  tx2gene <- exon2gene$island[sel]
  transcripts <- vector(length=length(islands), mode="list")
  names(transcripts) <- names(islands)
  tmp <- split(extxs, tx2gene)
  transcripts[names(tmp)] <- tmp
  
  ans <- new("annotatedGenome", islands=islands, transcripts=transcripts, exon2island=exon2island, exonsNI=exons, aliases=DB@aliases, genomeVersion=DB@genomeVersion, dateCreated=Sys.Date(), denovo=TRUE)
  ans
}

mergeStrDenovo <- function(plus, minus){  
  nullplus <- sapply(plus@transcripts, is.null)
  newplus <- by(names(plus@islands[nullplus]@unlistData), rep(1:sum(nullplus), elementLengths(plus@islands[nullplus])), function(x) paste(x, collapse='.'))
  names(minus@islands) <- 1:length(minus@islands)
  nullminus <- unlist(lapply(minus@transcripts, is.null))
  newminus <- by(names(minus@islands[nullminus]@unlistData), rep(1:sum(nullminus), elementLengths(minus@islands[nullminus])), function(x) paste(sort(x), collapse='.'))
  names(newminus) <- names(minus@islands)[nullminus]
  common <- names(newminus)[!(newminus %in% newplus)]
  allislands <- c(plus@islands, minus@islands[!(names(minus@islands) %in% common)])
  names(allislands) <- 1:length(allislands)
  alltrans <- c(plus@transcripts, minus@transcripts[!(names(minus@islands) %in% common)])
  names(alltrans) <- 1:length(alltrans)
  allexonsNI <- unique(c(plus@exonsNI, minus@exonsNI))
  ex2is <- as.data.frame(allexonsNI)
  ex2is$id <- names(allexonsNI)

  ans <- new("annotatedGenome", aliases=plus@aliases, denovo=TRUE, exonsNI=allexonsNI, transcripts=alltrans, exon2island=ex2is, dateCreated=Sys.Date(), genomeVersion=plus@genomeVersion, islands=allislands)
  ans
}

createDenovoGenome <- function(reads, DB, stranded=FALSE,  minLinks=2, maxLinkDist=1e+05, maxDist=1000, minConn=2, minJunx=3, minLen=12, mc.cores=1){
  cat("Finding new exons\n")
  somex <- NULL
  DBplus <- NULL
  DBminus <- NULL
  if(any(strand(DB@islands@unlistData)=='+')) DBplus <- genomeBystrand(DB, "+")
  if(any(strand(DB@islands@unlistData)=='-')) DBminus <- genomeBystrand(DB, "-")
  if(reads@stranded){
    newexplus <- NULL
    newexminus <- NULL
    if(!is.null(DBplus)) newexplus <- findNewExons(reads@plus, reads@pjunx, DBplus, minConn=minConn, minJunx=minJunx, minLen=minLen)
    if(!is.null(DBminus)) newexminus <- findNewExons(reads@minus, reads@mjunx, DBminus, minConn=minConn, minJunx=minJunx, minLen=minLen)
    if(!is.null(newexplus) | !is.null(newexminus)) somex <- 1
  } else {
    newex <- findNewExons(reads@pbam, reads@junx, DB, minConn=minConn, minJunx=minJunx, minLen=minLen)
    if(!is.null(newex)) somex <- 1
  }
  
  if(!is.null(somex)){
    if(!is.null(DBplus)){
      cat("Done...\nCreating denovo genome for positive strand\n")
      if(reads@stranded) denovoplus <- assignExons2Gene(newexplus, DBplus, reads@plus, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
      else {
        newexplus <- newex
        if(!is.null(DBminus)) newexplus <- newex[!(newex %in% DBminus@exonsNI[!(DBminus@exonsNI %in% DBplus@exonsNI)])]
        denovoplus <- assignExons2Gene(newexplus, DBplus, reads@pbam, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
      }
    }
    if(!is.null(DBminus)){
      cat("Done...\nCreating denovo genome for negative strand\n")
      if(reads@stranded) denovominus <- assignExons2Gene(newexminus, DBminus, reads@minus, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
      else {
        newexminus <- newex
        if(!is.null(DBplus)) newexminus <- newex[!(newex %in% DBplus@exonsNI[!(DBplus@exonsNI %in% DBminus@exonsNI)])]
        denovominus <- assignExons2Gene(newexminus, DBminus, reads@pbam, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
      }
    }
    if(!is.null(DBminus) & !is.null(DBplus)){
      cat("Done...\nMerging denovo genome\n")
      denovo <- mergeStrDenovo(denovoplus, denovominus)
    } else {
      if(!is.null(DBplus)) denovo <- denovoplus
      if(!is.null(DBminus)) denovo <- denovominus
    }
  } else {
    denovo <- DB
    denovo@denovo <- TRUE
  }
  denovo

}

 
