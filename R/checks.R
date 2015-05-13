simMAEcheck <- function(nsim, islandid, burnin=1000, pc, distr, readLength.pilot, eset.pilot, usePilot=FALSE, retTxsError=FALSE, genomeDB, mc.cores=1, mc.cores.int=1, verbose=FALSE) {
  if (missing(islandid)) { 
    islandid <- names(genomeDB@transcripts) 
  }
  distr.pilot <- distr
  reads = getNreads(pc)
  n = sum(reads[islandid])
  r = readLength.pilot
  f = mean(rep(as.numeric(names(distr@lenDis)), distr@lenDis))
  readLength = r
  U <- NULL
  txe <- list()
  pcs <- list()
  disArray <- vector("list", length(f))
  for (i in 1:length(disArray)) disArray[[i]] <- setfragLength(distr, fragLength=f[i])
  if (missing(pc)) {
    if (retTxsError) stop("Option retTxsError==TRUE is only available when argument pc is specified")
    if (verbose) cat("Simulating pilot data...\n")
    th <- exprs(eset.pilot)[sample(1:nrow(eset.pilot),length(islandid),replace=TRUE),sample(1:ncol(eset.pilot),1)]
    th <- 2^(th-max(th)); th <- th/sum(th)
    nSimReads <- rmultinom(n=1, size=10^7, prob=th)[,1]
    names(nSimReads) <- islandid
    simpisDirichlet <- function(z) { if (length(z)==1) { z <- 1 } else { z <- rgamma(length(z), shape=1/length(z)); z <- z/sum(z) }; return(z) }
    txs <- unlist(lapply(genomeDB@transcripts[islandid], names))
    pis <- unlist(lapply(genomeDB@transcripts[islandid], simpisDirichlet))
    names(pis) <- txs
    pc <- simReads(islandid,nSimReads=nSimReads,pis=pis,rl=readLength.pilot,distrs=distr.pilot,genomeDB=genomeDB,writeBam=FALSE,verbose=FALSE,mc.cores=mc.cores,seed=1)
  }
  for (j in 1:length(n)) {
    if (length(disArray[[j]])==1) { distrs <- disArray[[j]] } else { distrs <- disArray[[j]][[sample(1:length(disArray[[j]]),1)]] }
    if(verbose) cat(paste("Generating posterior samples j =",j, "\n"))
    pis <- simPost(islandid=islandid, nsim=nsim, distrs=distr.pilot, genomeDB=genomeDB, pc=pc, readLength=readLength.pilot, mc.cores=mc.cores.int*mc.cores, verbose=verbose)
    if(verbose) cat(paste("Running simulations for j =",j, "\n"))
    if ((mc.cores.int>1) && requireNamespace("parallel", quietly=TRUE)) {
        res <- parallel::mclapply(1:nsim, function(i){
          readYield <- runif(1,0.8,1.2) #actual reads +/- 20% within target
          pmapped <- (1-probNonMappable(readLength)) * runif(1,0.6,0.9)  #mapped reads 60%-90% of mappable reads
          N <- round(n[j]*readYield*pmapped)
          sim.pc <-  simPostPred(islandid=islandid, nreads=N, pis=pis[i,], pc=pc, distrs=distrs, rl=r[j], genomeDB=genomeDB, verbose=verbose)
          if(usePilot) sim.pc$pc <- mergePCs(pcs=list(sim.pc$pc,pc), genomeDB=genomeDB)
          sim.exp <- exprs(calcExp(islandid=islandid, distrs=distrs, genomeDB=genomeDB, pc=sim.pc$pc, readLength=r[j], rpkm=FALSE, mc.cores=mc.cores))
          list(maes=abs(sim.exp[colnames(pis),]-pis[i,]), pc=sim.pc$pc)
        }, mc.cores=mc.cores.int, mc.preschedule=FALSE)
    }
    else {
      res <- lapply(1:nsim, function(i){
        readYield <- runif(1,0.8,1.2) #actual reads +/- 20% within target
        pmapped <- (1-probNonMappable(readLength)) * runif(1,0.6,0.9)  #mapped reads 60%-90% of mappable reads
        N <- round(n[j]*readYield*pmapped)
        sim.pc <-  simPostPred(islandid=islandid, nreads=N, pis=pis[i,], pc=pc, distrs=distrs, rl=r[j], genomeDB=genomeDB, verbose=verbose)
        if(usePilot) sim.pc$pc <- mergePCs(pcs=list(sim.pc$pc,pc), genomeDB=genomeDB)
        sim.exp <- exprs(calcExp(islandid=islandid, distrs=distrs, genomeDB=genomeDB, pc=sim.pc$pc, readLength=r[j], rpkm=FALSE))
        list(maes=abs(sim.exp[colnames(pis),]-pis[i,]), pc=sim.pc$pc)
      })
    }
    Un <- do.call(cbind, lapply(res, "[[", "maes"))
    df <- data.frame(MAE= colMeans(Un), Nreads=rep(n[j], nsim), ReadLength=rep(r[j], nsim), frLength=rep(f[j], nsim))
    U <- rbind(U, df)
    
    if(retTxsError) {
      rownames(Un) <- colnames(pis)
      txe[[j]] <- Un
      pcs[[j]] <- lapply(res, '[[', 'pc')
    }
    
    if (j==length(n)) { 
      l1 = res
      nrvector = lapply(l1, function(x) getNreads(x$pc)[islandid])
      nrvector = lapply(nrvector, function(x) x[sort(names(x))])
      rown  = names(nrvector[[1]])
      df = data.frame(matrix(unlist(nrvector), nrow=length(nrvector[[1]])))
      rownames(df) = rown
      minim = apply(df, 1, min)
      maxim = apply(df, 1, max)
      nrpilot = getNreads(pc)[rown]
      s1 = sum(nrpilot >= minim & nrpilot <= maxim)
    }
    
  }
  
  if(retTxsError){
    names(txe) <- paste(n, r, f, sep='-')
    names(pcs) <- paste(n, r, f, sep='-')
    ans <- list(txe=txe, U=U, pcs=pcs, pis=pis)
  } 
  
  else {
    expected = round(length(islandid)*(nsim-1)/(nsim+1))
    ans <- list(U=U, mod_check=data.frame(Expected=expected, Observed=s1))
  }
  ans
}
