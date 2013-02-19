mergeExp <- function(..., sampleNames,  keep=c('transcript','gene')) {
  esets <- list(...)
  if (length(unique(sapply(esets,nrow))) != 1) stop('Number of rows do not match')
  if (any(sapply(esets,ncol) != 1)) stop("All ExpressionSet objects in '...' should have 1 column")
  if (missing(sampleNames)) sampleNames <- paste('Sample',1:length(esets),sep='')
  if (length(esets)>1) {
    sel <- !(keep %in% c('transcript','gene'))
    fData(esets[[1]]) <- fData(esets[[1]])[,keep]
    fvarLabels(esets[[1]])[sel] <- paste(sampleNames[1],fvarLabels(esets[[1]])[sel],sep='.')
    sampleNames(esets[[1]]) <- sampleNames[1]
    n <- featureNames(esets[[1]])
    for (i in 2:length(esets)) {
      if (any(!(keep %in% fvarLabels(esets[[i]])))) stop("Some variables in argument 'keep' were not found")
      fData(esets[[i]]) <- fData(esets[[i]])[,keep]
      sel <- !(keep %in% c('transcript','gene'))
      fvarLabels(esets[[i]])[sel] <- paste(sampleNames[i],fvarLabels(esets[[i]])[sel],sep='.')
      esets[[i]] <- esets[[i]][n,]
      sampleNames(esets[[i]]) <- sampleNames[i]
    }
    ans <- do.call("combine",esets)
  } else {
    ans <- esets[[1]]
  }
}
