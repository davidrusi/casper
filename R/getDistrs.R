transDistr <- function(distr, newmean, readLength){
  lprobs <- distr@lenDis
  totF <- sum(lprobs)
  nmean <-  mean(rep(as.numeric(names(lprobs)), lprobs))
  nmean <- newmean - nmean
  names(lprobs) <- as.numeric(names(lprobs)) + round(nmean)
  nprobs <- names(lprobs)
  lprobs <- lprobs[as.numeric(nprobs)>0]
  nprobs <- names(lprobs)
  distr@lenDis <- as.array(lprobs)
  if(min(as.numeric(names(lprobs)))>=readLength) {
    probs <- as.integer(ceiling(dpois(readLength:(min(as.numeric(names(lprobs)))-1), lambda=newmean)*totF))
    names(probs) <- as.character(readLength:(min(as.numeric(names(lprobs)))-1))
   #Combine into new prob vector
    nprobs <- c(names(probs), nprobs)
    lprobs <- c(probs, lprobs)
    distr@lenDis <- as.array(lprobs)
    names(distr@lenDis) <- nprobs
  }
  distr
}

startDist <- function(st,fragLength,txLength, nreads=NULL) {
                                        # Estimate relative start distribution under left??truncation (st < 1 ?? fragLength/txLength)
                                        # ?? st: relative start (i.e. start/txLength)
                                        # ?? fragLength: fragment length
                                        # ?? txLength: transcript length
                                        # Output: cumulative probability function (actually, a linear interpolation)
  if(!is.null(nreads)){
    if(nreads<length(st)){
      ran <- sample(1:length(st), size=nreads, replace=T)
      st <- st[ran]
      fragLength <- fragLength[ran]
      txLength <- txLength[ran]
    }
  }
  
  trunc <- 1-fragLength/txLength
  sel <- trunc>(st+1e-10)
  trunc <- 1-trunc[sel]
  st <- 1-st[sel]

  fit <- summary(survfit(Surv(time=trunc,time2=st,event=rep( TRUE,length(st))) ~ 1))
  s <- 1-fit$time
  pcum <- fit$surv
  if(length(s)>1){
    f <- approxfun(s,pcum)
    sseq <- seq(0,1,.001)
    startcdf <- f(sseq)
    startcdf[1] <- 0; startcdf[length(startcdf)] <- 1
    f <- approxfun(sseq[!is.na(startcdf)], startcdf[!is.na(startcdf)])
  } else f <- approxfun(c(0,1), c(0,0))
  return(f)
}

firstBamReads <- function(bam, nreads) {
  if (nreads < length(bam$qname)) {
    id <- unique(bam$qname[1:nreads])
    sel <- bam$qname %in% id
    bam <- lapply(bam,function(z) z[sel])
  }
  return(bam)
}


getDistrsSplit <- function(DBsplit, pbamSplit) {
  #Estimate distributions from split pbam object. Runs getDistrs on each chromosome separately and then combines
  distrsplit <- ans <- vector("list",length(DBsplit))
  names(ans) <- names(DBsplit)
  for (i in 1:length(distrsplit)) {
    distrsplit[[i]] <- vector("list",length(pbamSplit))
    for (j in 1:length(pbamSplit)) { distrsplit[[i]][[j]] <- getDistrs(DB=DBsplit[[i]], pbam=pbamSplit[[j]], verbose=FALSE); cat('.') }
    cat('\n')
  }
  for (i in 1:length(ans)) ans[[i]] <- mergeDisWr(distrsplit[[i]])
  return(ans)
}


getDistrs <- function(DB, bam, pbam, islandid=NULL, verbose=FALSE, nreads=4*10^6, readLength, min.gt.freq = NULL, tgroups=5, mc.cores=1){    

    if('gene_type' %in% colnames(DB@aliases) & !'gene_type_merged' %in% colnames(DB@aliases)){
        types <- table(DB@aliases$gene_type)
        if(!is.null(min.gt.freq)){
            freqs <- types/sum(types)
            mer <- names(freqs[freqs < min.gt.freq])
        } else {
            freqs <- sort(types, decreasing=T)
            mer <- names(freqs[tgroups:length(freqs)])
        }
        newgt <- DB@aliases$gene_type
        newgt[newgt %in% mer] <- 'merged'
        DB@aliases$gene_type_merged <- as.character(newgt) 
    }
    
    if (missing(pbam)) {
        if('gene_type_merged' %in% colnames(DB@aliases)){
            types <- unique(DB@aliases$gene_type_merged)
            if(mc.cores > 1) {
                if ('parallel' %in% loadedNamespaces()) {
                    ans <- parallel::mclapply(types, function(x) {
                        isl <- unique(as.character(DB@aliases$island_id[DB@aliases$gene_type_merged==x]))
                        getDistrsFromBam(DB=DB, bam=bam, islandid=islandid, verbose=verbose, nreads=nreads, readLength=readLength, selislands=isl)
                    }, mc.cores=min(mc.cores, length(types)))
                } else stop('parallel library has not been loaded!')
            } else {
                ans <- lapply(types, function(x) {
                    isl <- unique(as.character(DB@aliases$island_id[DB@aliases$gene_type_merged==x]))
                    getDistrsFromBam(DB=DB, bam=bam, islandid=islandid, verbose=verbose, nreads=nreads, selislands=isl)
                })
            }
            ld <- lapply(ans, truncLenDis)
            names(ld) <- types
            stDis <- lapply(ans, slot, "stDis")
            names(stDis) <- types
            ans <- new("readDistrsList",lenDis=ld,stDis=stDis)
        } else {
            ans <- getDistrsFromBam(DB=DB, bam=bam, islandid=islandid, verbose=verbose, nreads=nreads, readLength=readLength)
            ans@lenDis <- truncLenDis(ans)
        }
    } else {
        if('gene_type_merged' %in% colnames(DB@aliases)){
            types <- unique(DB@aliases$gene_type_merged)
            if(mc.cores > 1) {
                if ('parallel' %in% loadedNamespaces()) {  
                    ans <- parallel::mclapply(types, function(x) {
                        isl <- unique(as.character(DB@aliases$island_id[DB@aliases$gene_type_merged==x]))
                        getDistrsFrompBam(DB=DB, pbam=pbam, islandid=islandid, verbose=verbose, nreads=nreads, selislands=isl)
                    }, mc.cores=min(mc.cores, length(types)))
                } else stop('parallel library has not been loaded!')   
            } else {
                ans <- lapply(types, function(x) {
                    isl <- unique(as.character(DB@aliases$island_id[DB@aliases$gene_type_merged==x]))
                    getDistrsFrompBam(DB=DB, pbam=pbam, islandid=islandid, verbose=verbose, nreads=nreads, selislands=isl)
                })
            }                
            ld <- lapply(ans, truncLenDis)
            names(ld) <- types
            stDis <- lapply(ans, slot, "stDis")
            names(stDis) <- types
            ans <- new("readDistrsList",lenDis=ld,stDis=stDis)
        } else {
            ans <- getDistrsFrompBam(DB=DB, pbam=pbam, islandid=islandid, verbose=verbose, nreads=nreads)
            ans@lenDis <- truncLenDis(ans)
        }
    }
    return(ans)
}

truncLenDis <- function(ans){
    #Truncate length distr to (median - 3*IQR, median + 6*IQR)
    if(length(ans@lenDis)>1) {
        w <- ans@lenDis/max(ans@lenDis)
        wcum <- cumsum(w/sum(w))
        q <- as.numeric(c(names(wcum)[which(wcum>.25)[1]],names(wcum)[which(wcum>.75)[1]]))
        iqr <- q[2]-q[1]
        sel <- (as.numeric(names(wcum))>= q[1]-3*iqr) & (as.numeric(names(wcum))<= q[2]+6*iqr)
        ans@lenDis <- as.array(ans@lenDis[sel])
    }
    return(ans@lenDis)
}


getDistrsFrompBam <- function(DB, pbam, islandid=NULL, verbose=FALSE, nreads=4*10^6, selislands=NULL){

  if (class(pbam) != 'procBam') stop('Argument pbam must be of class procBam')
  
  if(!is.null(selislands)) {
      exonsRD <- DB@exonsNI[names(DB@islands[selislands]@unlistData)]
  } else exonsRD <- DB@exonsNI
  if(!any(unique(seqnames(pbam@pbam)) %in% levels(seqnames(exonsRD)))) stop('Different chromosome names in pbam and genome')

  #Find fragment length distribution for fragments aligning to exons larger than 1000 bases  
  if(verbose) cat("Calculating fragment length distribution\n")

  #Remove reads with >2 appearances
  sel <- c(TRUE, pbam@pbam$names[-1]!=pbam@pbam$names[-length(pbam@pbam)] | (pbam@pbam$rid[-1] != pbam@pbam$rid[-length(pbam@pbam)]))
  pbam@pbam <- pbam@pbam[sel,]

  sel <- pbam@pbam$rid==1
  sel1 <- pbam@pbam$rid==2
  sel2 <- as.character(seqnames(pbam@pbam[sel])) == as.character(seqnames(pbam@pbam[sel1]))
  
  frags <- GRanges(IRanges(start(pbam@pbam)[sel][sel2], end(pbam@pbam)[sel1][sel2]), seqnames=seqnames(pbam@pbam)[sel][sel2])
  n <- levels(seqnames(frags))[levels(seqnames(frags)) %in% levels(seqnames(exonsRD))]
  fragsL<-frags[levels(seqnames(frags)) %in% n]
  over <- suppressWarnings(findOverlaps(fragsL, subset(exonsRD, width(exonsRD)>1000), type="within"))
  if (length(subjectHits(over))<10) {
    #over <- suppressWarnings(findOverlaps(fragsL, subset(exonsRD, width(exonsRD)>500), type="within"))
    ld <- array(0); names(ld) <- '300'
  } else {
    ld <- table(width(fragsL)[queryHits(over)])
    ld <- ld[ld/sum(ld) > 0.0001]
    ns <- names(ld)
    ld <- array(ld)
    names(ld) <- ns
  }

  #Find fragment start distribution for fragments aligning to transcripts in genes with only one anno tated tx
  if(verbose) cat("Calculating fragment start distribution\n")
  if(is.null(islandid) & is.null(selislands)){
    sel <- unlist(lapply(DB@transcripts, length))==1
  } else sel= islandid
  if(!is.null(selislands)){
      sel <- unlist(lapply(DB@transcripts, length))==1 & names(DB@transcripts) %in% selislands
  }
  oneTx <- DB@transcripts[sel]
  oneTx <- unlist(oneTx, recursive=F)
  names(oneTx) <- sub("[0-9]+\\.", "", names(oneTx))
  oneExons <- exonsRD[names(exonsRD) %in% unlist(oneTx)]
  oneExons <- oneExons[as.character(unlist(oneTx))]
  
  islandStrand <- as.character(strand(DB@islands@unlistData))[cumsum(c(1, elementNROWS(DB@islands)[-length(DB@islands)]))]
  names(islandStrand) <- names(DB@islands)
  

  exon_rank <- unlist(lapply(oneTx, function(x) 1:length(x)))
  exon_strand <- islandStrand[as.character(DB@exon2island$island[match(names(oneExons), rownames(DB@exon2island))])]
  strand(oneExons) <- exon_strand
  values(oneExons)$exon_rank <- exon_rank
  txid <- cumsum( values(oneExons)$exon_rank==1 )
  wid <- end(oneExons) - start(oneExons) + 1
  values(oneExons)$exstnogap <- unlist(tapply(wid, txid, function(z) c(0, cumsum(z[-length(z)]))))
  values(oneExons)$txlength <- unlist(tapply((end(oneExons) - start(oneExons) + 1),INDEX=txid,function(z) rep(sum(z),length(z))))
  
  frags <- frags[width(frags)<max(width(oneExons))]
  over <- suppressWarnings(findOverlaps(frags, oneExons, type="within"))

  if (length(over)>0) {
    exstnogap <- values(oneExons)$exstnogap[subjectHits(over)]
    txlength <- values(oneExons)$txlength[subjectHits(over)]
    exst <- start(oneExons)[subjectHits(over)]
    exen <- end(oneExons)[subjectHits(over)]
    readst <- start(frags)[queryHits(over)]
    readen <- end(frags)[queryHits(over)]
    str <- as.character(strand(oneExons))[subjectHits(over)]
    frlen <- width(frags)[queryHits(over)]
     
    stDis <- double(length(readst))
    sel <- str=='+';  stDis[sel] <- (exstnogap[sel]+readst[sel]-exst[sel])/txlength[sel]
    sel <- str=='-'; stDis[sel] <- (exstnogap[sel]+exen[sel]-readen[sel])/txlength[sel]
    stDis <- startDist(stDis, frlen, txlength)
  } else {
    stDis <- function(z) return(0)
  }
  
  new("readDistrs",lenDis=ld,stDis=stDis)
}



getDistrsFromBam <- function(DB, bam, islandid=NULL, verbose=FALSE, nreads=4*10^6, readLength, selislands=NULL){

  if (!all(c('qname','rname','pos','mpos') %in% names(bam))) stop('bam must contain elements qname, rname, pos, mpos')

  #Select a sample of reads
  if (length(bam[['qname']] > nreads)) bam <- firstBamReads(bam, nreads=nreads)
  #frags <- GRanges(IRanges(start=bam$pos),seqnames=bam$rname)
  #islands <- as.data.frame(DB@islands)
  #islandrg <- sqldf("select element, seqnames, min(start), max(end) from islands group by element")
  #islandrg <- GRanges(IRanges(islandrg[,3],islandrg[,4]), seqnames=islandrg[,2])
  #sel <- frags %over% islandrg
  
  #Find fragment length distribution for fragments aligning to exons larger than 1000 bases  
  if(verbose) cat("Calculating fragment length distribution\n")
  if(!is.null(selislands)) {
      exonsRD <- DB@exonsNI[names(DB@islands[selislands]@unlistData)]
  } else exonsRD <- DB@exonsNI
  d <- bam$mpos - bam$pos
  sel <- d<0; n <- bam$qname[sel]; sp <- bam$rname[sel]; names(sp) <- n
  en <- bam$pos[sel]+readLength-1; names(en) <- n  #faster than bam$pos[sel]+bam$qwidth[sel]-1
  sel <- d>0; st <- bam$pos[sel]; names(st) <- bam$qname[sel];
  sel <- match(n, names(st))
  st <- st[sel[!is.na(sel)]]
  en <- en[names(st)]; sp <- sp[names(st)]
  sel <- st<en; st <- st[sel]; en <- en[sel]; sp <- sp[sel]
  
  if(!any(unique(sp) %in% names(exonsRD))) {
    if(any(grepl('chr', names(exonsRD)))) {
      sp <- as.factor(paste('chr', as.character(sp), sep=''))
    } else if(any(grepl('chr', unique(sp)))) sp <- sub('chr', '', sp)
  }

  if(!any(unique(sp) %in% levels(seqnames((exonsRD))))) {
    if(any(grepl('chr', levels(seqnames(exonsRD))))) {
      sp <- as.factor(paste('chr', as.character(sp), sep=''))
    } else if(any(grepl('chr', unique(sp)))) sp <- sub('chr', '', sp)
  }
  
  if(!any(unique(sp) %in% levels(seqnames(exonsRD)))) stop('Different chromosome names in bam and genome')
  frags <- GRanges(IRanges(start=st,end=en),seqnames=sp)
  n <- levels(seqnames(frags))[levels(seqnames(frags)) %in% levels(seqnames(exonsRD))]
  fragsL<-frags[levels(seqnames(frags)) %in% n]
  over <- suppressWarnings(findOverlaps(fragsL, subset(exonsRD, width(exonsRD)>1000), type="within"))
  if(length(subjectHits(over))==0) over <- suppressWarnings(findOverlaps(fragsL, subset(exonsRD, width(exonsRD)>500), type="within"))
  ld<-table(width(fragsL)[queryHits(over)])
  ld <- ld[ld/sum(ld) > 0.0001]
  ns <- names(ld)
  ld <- array(ld)
  names(ld) <- ns
  

  #Find fragment start distribution for fragments aligning to transcripts in genes with only one annotated tx

  if(verbose) cat("Calculating fragment start distribution\n")
  if(is.null(islandid)){
    sel <- unlist(lapply(DB@transcripts, length))==1
  } else sel= islandid
  oneTx <- DB@transcripts[sel]
  oneTx <- unlist(oneTx, recursive=F)
  names(oneTx) <- sub("[0-9]+\\.", "", names(oneTx))
  oneExons <- exonsRD[names(exonsRD) %in% unlist(oneTx)]
  oneExons <- oneExons[as.character(unlist(oneTx))]
  
  islandStrand <- as.character(strand(DB@islands@unlistData))[cumsum(c(1, elementNROWS(DB@islands)[-length(DB@islands)]))]
  names(islandStrand) <- names(DB@islands)
  

  exon_rank <- unlist(lapply(oneTx, function(x) 1:length(x)))
  exon_strand <- islandStrand[as.character(DB@exon2island$island[match(names(oneExons), rownames(DB@exon2island))])]
  strand(oneExons) <- exon_strand
  values(oneExons)$exon_rank <- exon_rank
  txid <- cumsum( values(oneExons)$exon_rank==1 )
  wid <- end(oneExons) - start(oneExons) + 1
  values(oneExons)$exstnogap <- unlist(tapply(wid, txid, function(z) c(0, cumsum(z[-length(z)]))))
  values(oneExons)$txlength <- unlist(tapply((end(oneExons) - start(oneExons) + 1),INDEX=txid,function(z) rep(sum(z),length(z))))
  
  frags <- frags[width(frags)<max(width(oneExons))]
  over <- suppressWarnings(findOverlaps(frags, oneExons, type="within"))

  exstnogap <- values(oneExons)$exstnogap[subjectHits(over)]
  txlength <- values(oneExons)$txlength[subjectHits(over)]
  exst <- start(oneExons)[subjectHits(over)]
  exen <- end(oneExons)[subjectHits(over)]
  readst <- start(frags)[queryHits(over)]
  readen <- end(frags)[queryHits(over)]
  str <- as.character(strand(oneExons))[subjectHits(over)]
  frlen <- width(frags)[queryHits(over)]
  
  stDis <- double(length(readst))
  sel <- str=='+';  stDis[sel] <- (exstnogap[sel]+readst[sel]-exst[sel])/txlength[sel]
  sel <- str=='-'; stDis[sel] <- (exstnogap[sel]+exen[sel]-readen[sel])/txlength[sel]

  stDis <- startDist(stDis, frlen, txlength)
  lenDis <- new("array", ld); names(lenDis)= names(ld)
  new("readDistrs",lenDis=lenDis,stDis=stDis)
}

