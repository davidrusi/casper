procExp <- function(distrs, genomeDB, pc, readLength, islandid, rpkm=FALSE, priorq=2, citype='none', niter=10^3, burnin=100, mc.cores=1, verbose=FALSE, totReads) {
  #Format input
  startcdf <- distrs@stDis(seq(0,1,.001))
  lendis <- as.double(distrs@lenDis/sum(distrs@lenDis))
  lenvals <- as.integer(names(distrs@lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  citype <- as.integer(switch(citype, none=0, asymp=1, exact=2))
  if (length(citype)==0) stop("citype should take the value 'none', 'asymp', or 'exact'")
  niter <- as.integer(niter)
  burnin <- as.integer(burnin)
  verbose <- as.integer(verbose)
  if ((niter<=burnin) & citype==2) stop("Too many burnin iterations specified. Decrease burnin or increase niter")
  if (missing(islandid)) islandid <- names(genomeDB@islands)[elementLengths(genomeDB@islands)>1]

  exons <- names(unlist(genomeDB@islands))
  exons <- as.integer(sub(".*\\.", "", exons))
  names(exons) <- rep(names(genomeDB@islands), elementLengths(genomeDB@islands))
  exons <- split(unname(exons), names(exons))
  exonwidth <- unlist(width(genomeDB@islands))
  names(exonwidth) <- rep(names(genomeDB@islands), elementLengths(genomeDB@islands))
  exonwidth <- split(unname(exonwidth), names(exonwidth))
  strand <- as.character(strand(genomeDB@islands@unlistData))[cumsum(c(1, elementLengths(genomeDB@islands)[-length(genomeDB@islands)]))]
  names(strand) <- names(genomeDB@islands)
  
  if (!all(islandid %in% names(exons))) stop('islandid not found in genomeDB@islands')
  if (!all(islandid %in% names(pc))) stop('islandid not found in pc')
  if (!all(islandid %in% names(genomeDB@transcripts))) stop('islandid not found in genomeDB@transcripts')
  
      #Define basic function
  
  f <- function(z) {
    islandid <- as.integer(z)
    exons <- exons[z]
    exonwidth <- exonwidth[z]
    transcripts <- genomeDB@transcripts[z]
    strand <- as.list(as.integer(ifelse(strand[z]=='+', 1, -1)))
    pc <- pc[z]
    ans <- calcKnownMultiple(exons=exons,exonwidth=exonwidth,transcripts=transcripts,islandid=as.list(islandid),pc=pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,priorq=priorq, strand=strand, citype=citype, niter=niter, burnin=burnin, verbose=verbose)
    if (citype==0) {
      ans <- lapply(ans, function(z) { res=vector("list",1); res[[1]]= z[[1]]; names(res[[1]])= z[[2]]; res })
    } else if (citype==1) {
      ans <- lapply(ans, function(z) { res=vector("list",2); res[[1]]= z[[1]]; res[[2]]= z[[3]]; names(res[[1]])= names(res[[2]])= z[[2]]; res })
    } else if (citype==2) {
      ans <- lapply(ans, function(z) { res=vector("list",3); res[[1]]= z[[1]]; res[[2]]= z[[3]]; res[[3]]= matrix(z[[4]],nrow=niter-burnin); names(res[[1]])= names(res[[2]])= z[[2]]; res })
    }
    ans
  }

      #Run
    sel <- !sapply(pc[islandid], is.null)
    all <- islandid
    islandid <- islandid[sel]
  if (verbose) cat("Obtaining expression estimates...\n")
    if (mc.cores>1 && length(islandid)>mc.cores) {
      if ('multicore' %in% loadedNamespaces()) {
                    #split into smaller jobs
        islandid <- split(islandid, cut(1:length(islandid), mc.cores))
        ans <- multicore::mclapply(islandid,f,mc.cores=mc.cores)
        ans <- unlist(ans, recursive=F)
        names(ans) <- unlist(islandid)
      } else stop('multicore library has not been loaded!')
    } else {
      ans <- f(islandid)
      names(ans) <- islandid
    }

    miss <- lapply(genomeDB@transcripts[all[!(all %in% names(ans))]], names)  
    #Format as ExpressionSet

    if (verbose) cat("Formatting output...\n")
    transcript <- c(unlist(lapply(ans, function(z) names(z[[1]]))), unname(unlist(miss)))
    gene <- c(rep(names(ans),sapply(ans,function(z) length(z[[1]]))), rep(names(miss), sapply(miss, length)))
    fdata <- data.frame(transcript=transcript, gene=gene)
    misse <- rep(NA, length(unlist(miss)))
    names(misse) <- unlist(miss)
    exprsx <- matrix(c(unlist(lapply(ans,'[[',1)), misse) ,ncol=1)
    rownames(exprsx) <- c(unlist(lapply(ans, function(z) names(z[[1]]))), names(misse))
    exprsx[exprsx==-1] <- NA
    exprsx <- round(exprsx,10)
    if (citype==1) {
      se <- c(unlist(lapply(ans,'[[',2)), misse)
      se[se==-1] <- NA
      alpha <- .05
      #strategy 1
      #fdata$ci95.high <- fdata$ci95.low <- rep(NA,nrow(fdata))
      #sel <- !is.na(se) & se!=0
      #fdata$ci95.high[sel] <- mapply(function(m,s) qtnorm(1-alpha/2,mean=m,sd=s,lower=0,upper=1), as.list(exprsx[sel,1]), as.list(se[sel]))
      #fdata$ci95.low[sel] <- mapply(function(m,s) qtnorm(alpha/2,mean=m,sd=s,lower=0,upper=1), as.list(exprsx[sel,1]), as.list(se[sel]))
      #sel <- !is.na(se) & se==0; fdata$ci95.high[sel] <- fdata$ci95.low[sel] <- exprsx[sel,1]
      #sel <- fdata$ci95.high < exprsx[,1] & !is.na(fdata$ci95.high); fdata$ci95.high[sel] <- exprsx[sel,1]
      #sel <- fdata$ci95.low > exprsx[,1] & !is.na(fdata$ci95.low); fdata$ci95.low[sel] <- exprsx[sel,1]
      #strategy 2
      fdata$ci95.low <- qnorm(alpha/2,mean=exprsx,sd=se)
      fdata$ci95.high <- qnorm(1-alpha/2,mean=exprsx,sd=se)
      sel <- round(fdata$ci95.low,10)<0 & !is.na(se); fdata$ci95.low[sel] <- 0
      fdata$ci95.high[sel] <- unlist(mapply(function(m,s) qtnorm(1-alpha/2,mean=m,sd=s,lower=0,upper=1), as.list(exprsx[sel,1]), as.list(se[sel])))
      sel <- round(fdata$ci95.high,10)>1 & !is.na(se); fdata$ci95.high[sel] <- 1
      fdata$ci95.low[sel] <- unlist(mapply(function(m,s) qtnorm(alpha/2,mean=m,sd=s,lower=0,upper=1), as.list(exprsx[sel,1]), as.list(se[sel])))
    }
    if (citype==2) {
      if(sum(unlist(lapply(ans, function(x) is.na(x[[3]]))))==0){ 
        q <- lapply(ans,function(z) apply(z[[3]],2,quantile,probs=c(.025,.975)))
        q <- t(do.call(cbind,q))
        q <- rbind(q, matrix(NA, nrow=length(misse), ncol=2))
        fdata$ci95.low <- q[,1]; fdata$ci95.high <- q[,2]
        tmp <- cbind(exprsx, fdata[,3:4])
        tmp <- as.data.frame(t(apply(tmp, 1, function(x){
          if(!any(is.na(x))){
            y <- x
            if(x[2]>x[1]) y[2]=y[1]
            if(x[3]<x[1]) y[3]=y[1]
          } else y <- x
          y
          })))
        fdata$ci95.low <- tmp$ci95.low
        fdata$ci95.high <- tmp$ci95.high
      }
    }
      
    rownames(exprsx) <- rownames(fdata) <- fdata$transcript
    fdata <- new("AnnotatedDataFrame",fdata)
    ans <- new("ExpressionSet",exprs=exprsx,featureData=fdata)

    #Return absolute expression levels
  if (rpkm) {
    nreads <- sapply(pc[unique(as.character(ans@featureData$gene))],sum)
    geneLength <- sapply(genomeDB@islands[unique(as.character(ans@featureData$gene))], function(x) sum(width(x)))
    nreads <- (nreads/totReads)/(geneLength/10^9)
    exprs(ans) <- exprs(ans)*nreads[as.character(ans@featureData$gene)]
  }
  ans
}

mergeExp <- function(minus, plus){
  exp <- rbind(exprs(plus), exprs(minus))
  fdata <- new("AnnotatedDataFrame", rbind(fData(plus), fData(minus)))
  ans <- new("ExpressionSet",exprs=exp,featureData=fdata)
  ans
}

calcExp <- function(distrs, genomeDB, pc, readLength, islandid, rpkm=FALSE, priorq=2, citype='none', niter=10^3, burnin=100, mc.cores=1, verbose=FALSE) {
 totReads <- sum(getNreads(pc))
  
  if (missing(readLength)) stop("readLength must be specified")
  if (genomeDB@denovo) stop("genomeDB must be a known genome")
  if (pc@denovo) stop("pc must be a pathCounts object from known genome")
  if(pc@stranded){
    if(missing(islandid)) islandid <- c(names(pc@counts$plus), names(pc@counts$minus))
    is <- as.character(strand(genomeDB@islands@unlistData))[cumsum(c(1, elementLengths(genomeDB@islands)[-length(genomeDB@islands)]))]
    names(is) <- names(genomeDB@islands)
    plusGI <- islandid[is[islandid]=="+"]
    plus <- NULL
    if(length(plusGI)>0) {
      plusDB <- genomeBystrand(genomeDB, "+")
      plus <- procExp(distrs, plusDB, pc=pc@counts$plus, readLength=readLength, islandid=plusGI, rpkm=rpkm, priorq=priorq, niter=niter, burnin=burnin, mc.cores=mc.cores, citype=citype, totReads=totReads)
    }
    minusGI <- islandid[is[islandid]=="-"]
    minus <- NULL
    if(length(minusGI)>0) {
      minusDB <- genomeBystrand(genomeDB, "-")
      minus <- procExp(distrs, minusDB, pc=pc@counts$minus, readLength=readLength, islandid=minusGI, rpkm=rpkm, priorq=priorq, niter=niter, burnin=burnin, mc.cores=mc.cores, citype=citype, totReads=totReads)
    } 
    if(is.null(plus) & is.null(minus)) stop("No counts in islandid genes")
    if(!(is.null(plus) | is.null(minus))) { ans <- mergeExp(plus, minus) }
    else { if(is.null(plus)) {ans <- minus } else ans <- plus }
  } else {
    if(missing(islandid)) islandid <- names(pc@counts[[1]])
    if(sum(!unlist(lapply(pc@counts[islandid], is.null)))) stop("No counts in islandid genes")
    ans <- procExp(distrs, genomeDB, pc=pc@counts[[1]], readLength=readLength, islandid=islandid, rpkm=rpkm, priorq=priorq, niter=niter, burnin=burnin, mc.cores=mc.cores, citype=citype, totReads=totReads)
  }

  ans
}
    
calcKnownMultiple <- function(exons, exonwidth, transcripts, islandid, pc, startcdf, lendis, lenvals, readLength, priorq, strand, citype, niter, burnin, verbose) {
  ans <- .Call("calcKnownMultiple",exons,exonwidth,transcripts,islandid,pc,startcdf, lendis, lenvals, readLength, priorq, strand, citype, niter, burnin, verbose)
  return(ans)
}


lhoodGrid <- function(pc, distrs, genomeDB, readLength, islandid, grid, priorq=2) {
  if (length(islandid)>1) stop("Only 1 gene is allowed, please specify a single gene in argument 'islandid'")
  if (missing(readLength)) stop("readLength must be specified")
  if (genomeDB@denovo) stop("genomeDB must be a known genome")
  if (pc@denovo) stop("pc must be a pathCounts object from known genome")

  #Format input
  startcdf <- distrs@stDis(seq(0,1,.001))
  lendis <- as.double(distrs@lenDis/sum(distrs@lenDis))
  lenvals <- as.integer(names(distrs@lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)

  exons <- as.integer(names(genomeDB@islands[[islandid]]))
  exonwidth <- width(genomeDB@islands[[islandid]])
  strand <- as.character(strand(genomeDB@islands@unlistData))[cumsum(c(1, elementLengths(genomeDB@islands)[-length(genomeDB@islands)]))]
  names(strand) <- names(genomeDB@islands)
  strand   <- strand[islandid]
  transcripts <- genomeDB@transcripts[[islandid]]
  if (length(transcripts)==1) stop("Single transcript specified, estimation not run")
  strand <- as.integer(ifelse(strand=='+', 1, -1))
  if (pc@stranded) {
     if (strand==1) pc <- pc@counts[['plus']][[islandid]] else pc <- pc@counts[['minus']][[islandid]]
  } else {
     pc <- pc@counts[[1]][[islandid]]
  }
  islandid <- as.integer(islandid)
  
  if (missing(grid)) {
    grid <- lapply(1:(length(transcripts)-1), function(z) seq(-10,10,length=10^(5/(length(transcripts)-1))))
    #grid <- lapply(1:length(transcripts), function(z) seq(0.0001,.9999,length=10^(5/length(transcripts))))
  } else {
    if (!is.list(grid)) stop("grid must be of class list")
    if (length(grid) != length(transcripts)-1) stop("grid must have one row less than the number of known transcripts")
  }
  milogit <- function(theta) cbind(rep(1,nrow(theta)),exp(theta))/(1+rowSums(exp(theta))) #matrix conversion
  gridmat <- t(milogit(expand.grid(grid)))
  
  ans <- .Call("lhoodGrid",gridmat,exons,exonwidth,transcripts,pc,startcdf,lendis,lenvals,readLength,priorq,strand)
  names(ans[[2]]) <- ans[[3]]; ans[[3]] <- NULL
  names(ans) <- c('logpos','piest','logpos.piest','S','pathprob','marginalLhood')
  rownames(ans$pathprob) <- names(ans$piest)
  ans$pathprob <- t(ans$pathprob)
  names(ans$marginalLhood) <- c('laplace','IS')
  ans$grid <- grid
  return(ans)  
}


