MDSplot <- function(rf, fac, k=2, ...) {
  if(!inherits(rf, "randomForest")) {
    stop(deparse(substitute(rf)), " must be a randomForest object")
  }
  if(is.null(rf$proximity)) {
    stop(deparse(substitute(rf)), " does not contain a proximity matrix")
  }
  op <- par(pty="s")
  on.exit(par(op))
  Rver <- substring(paste(R.version[c("major", "minor")], collapse="."), 1, 3)
  cmdscale <- if (as.numeric(Rver) >= 1.9) stats:::cmdscale else mva:::cmdscale
  rf.mds <- cmdscale(1 - rf$proximity, eig=TRUE, k=k)
  colnames(rf.mds$points) <- paste("Dim", 1:k)
  nlevs <- length(levels(fac))
  if( require(RColorBrewer) && nlevs < 12)
    pal <- brewer.pal(nlevs,"Set1")
  else
    pal <- rainbow(nlevs)
  if(k <= 2) {
    plot(rf.mds$points, col=pal[as.numeric(fac)], pch=20, ...)
  } else {
    pairs(rf.mds$points, col=pal[as.numeric(fac)], pch=20, ...)
  }
  invisible(rf.mds)
}
