setMethod("getRoc", signature(simTruth='matrix',decision='matrix'), function(simTruth, decision) {
  tp <- colSums(simTruth & decision,na.rm=TRUE)
  fp <- colSums(!simTruth & decision,na.rm=TRUE)
  tn <- colSums(!simTruth & !decision,na.rm=TRUE)
  fn <- colSums(simTruth & !decision,na.rm=TRUE)
  data.frame(tp=tp, fp=fp, tn=tn, fn=fn, p=tp+fp, fdr=fp/(fp+tp), pow=tp/(tp+fn))
}
)

setMethod("getRoc", signature(simTruth='numeric',decision='numeric'), function(simTruth, decision) {
  getRoc(as.logical(simTruth), as.logical(decision))
}
)

setMethod("getRoc", signature(simTruth='logical',decision='logical'), function(simTruth, decision) {
  roc <- data.frame(tp=sum(simTruth & decision,na.rm=TRUE), fp=sum(!simTruth & decision,na.rm=TRUE), tn=sum(!simTruth & !decision,na.rm=TRUE), fn=sum(simTruth & !decision,na.rm=TRUE))
  roc$p <- roc$tp+roc$fp; roc$fdr <- roc$fp/(roc$fp+roc$tp); roc$pow <- roc$tp/(roc$tp+roc$fn)
  return(roc)
}
)
