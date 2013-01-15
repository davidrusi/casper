#Remove short inserts
rmShortInserts <- function(bam, isizeMin=100) {
  sel <- abs(bam$isize)>isizeMin
  xs <- bam$tag$XS[sel]
  ans <- lapply(bam[names(bam)!='tag'],function(z) z[sel])
  ans$tag <- list(XS=xs)
  return(ans)
}
