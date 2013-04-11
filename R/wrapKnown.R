
mergePCWr <- function(x, genomeDB){
  pc <- lapply(x, '[[', 'pc')
  tmp <- vector(mode='list', length=length(genomeDB@transcripts))
  names(tmp) <- names(genomeDB@transcripts)
  pcs <- unlist(lapply(pc, function(y) y@counts[[1]]), recursive=F)
  tmp[names(pcs)] <- pcs
  new("pathCounts", counts=list(tmp), stranded=x[[1]]$pc@stranded, denovo=FALSE)
}

mergeDisWr <- function(dis){
  minlen <- min(unlist(lapply(dis, function(x) min(as.numeric(names(x@lenDis))))))
  maxlen <- max(unlist(lapply(dis, function(x) max(as.numeric(names(x@lenDis))))))
  tmp <- lapply(dis, function(x){
    y <- x@lenDis
    tmp <- vector(mode='numeric', length=maxlen-minlen+1)
    names(tmp) <- minlen:maxlen
    tmp[names(y)] <- y
    tmp
  } )
  tmp <- do.call(cbind, tmp)
  tmp <- as.array(rowSums(tmp))
  distr <- new('readDistrs', lenDis=tmp)
                                        #stDis
  th <- seq(0,1,length=10000)
  tmp <- lapply(dis, function(x){
    probs <- x@stDis(th)
    ll <- sample(th, prob=probs, replace=TRUE, size=10000)
    ll
  })
  tmp <- unlist(tmp)
  tmp <- density(tmp)
  tmp <- approxfun(tmp$x, tmp$y)
  distr@stDis <- tmp
  distr
}


wrapKnown <- function(bamFile, verbose=FALSE, seed=1, mc.cores.int=1, mc.cores=1, genomeDB, readLength, rpkm=TRUE, priorq=2, priorqGeneExpr=2, citype='none', niter=10^3, burnin=100) {
  what <- scanBamWhat(); what <- what[!(what %in% c('seq','qual', 'flag'))]
  t <- scanBamHeader(bamFile)[[1]][["targets"]]
  which <- GRanges(names(t), IRanges(1, unname(t)))
  which <- which[!grepl("_",as.character(seqnames(which)))]
  which <- which[!as.character(seqnames(which))=='chrM']
  flag <- scanBamFlag(isPaired=TRUE,hasUnmappedMate=FALSE)
  if (mc.cores>1 ){
    if ('multicore' %in% loadedNamespaces()) {
      ans <- multicore::mclapply(1:length(which), function(i){
        param <- ScanBamParam(flag=flag,what=what, which=which[i], tag='XS')
        cat("Processing chromosome: ", as.character(seqnames(which[i])), "\n")
        bam <- scanBam(file=bamFile,param=param)
        pbam <- procBam(bam=bam, stranded=FALSE, seed=as.integer(seed), verbose=verbose)[[1]]
        distr <- getDistrs(DB=genomeDB, bam=bam[[1]], verbose=verbose)
        rm(bam)
        gc()
        pc <- pathCounts(reads=pbam, DB=genomeDB, mc.cores=mc.cores, verbose=verbose)
        list(pbam=pbam, distr=distr, pc=pc)
      }, mc.cores=mc.cores.int)
    } else stop('multicore library has not been loaded!')
  } else {
    ans <- lapply(1:length(which), function(i){
      param <- ScanBamParam(flag=flag,what=what, which=which[i], tag='XS')
      cat("Processing chromosome: ", as.character(seqnames(which[i])), "\n")
      bam <- scanBam(file=bamFile,param=param)
      pbam <- procBam(bam=bam, stranded=FALSE, seed=as.integer(seed), verbose=verbose)[[1]]
      distr <- getDistrs(DB=genomeDB, bam=bam[[1]], verbose=verbose)
      rm(bam)
      gc()
      pc <- pathCounts(reads=pbam, DB=genomeDB, mc.cores=mc.cores, verbose=verbose)
      list(pbam=pbam, distr=distr, pc=pc)
    })
  }
  allpbam <- lapply(ans, '[[', 'pbam')
  allpc <- mergePCWr(ans, genomeDB)
  alldistr <- suppressWarnings(mergeDisWr(lapply(ans, '[[', 'distr')))
  exp <- calcExp(distrs=alldistr, genomeDB=genomeDB, pc=allpc, readLength=readLength, rpkm=rpkm, priorq=priorq, priorqGeneExpr=priorqGeneExpr, citype=citype, niter=niter, burnin=burnin, mc.cores=mc.cores, verbose=verbose)
  list(pc=allpc, distr=alldistr, exp=exp, pbam=allpbam)
}

