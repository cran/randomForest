var.imp.plot <- function(x, sort=TRUE,
                         n.var=min(30, if(is.null(dim(x$importance)))
                           length(x$importance) else nrow(x$importance)),
                         main=deparse(substitute(x)), ...) {
  if(!inherits(x, "randomForest"))
    stop("This function only works for objects of class `randomForest'")
  if(is.null(dim(x$importance))) {
    if(sort) {
      ord <- order(x$importance, decreasing=TRUE)[1:n.var]
      dotchart(rev(x$importance[ord]), xlab="Importance", pch=16,
           ylab="", main=main, xlim=c(0, max(x$importance[ord])), ...)
    } else {
      dotchart(x$importance, xlab="Importance", pch=16,
           ylab="", main=main, xlim=c(0, x$importance), ...)
    }
  } else {
    n.meas <- ncol(x$importance)
    op <- par(mfrow=c(1, 2), mar=c(4, 5, 4, 1), mgp=c(2, .8, 0),
              oma=c(0,0,2,0))
    on.exit(par(op))
  
    for(i in 1:n.meas) {
      if(sort) {
        ord <- order(x$importance[,i], decreasing=TRUE)[1:n.var]
        maximp <- max(x$importance[ord, i])
        dotchart(rev(x$importance[ord, i]), pch=16, xlab="Importance",
                 ylab="", main=colnames(x$importance)[i], xlim=c(0, maximp),
                 ...)
      } else {
        dotchart(x$importance[,i], pch=16, xlab="Importance", ylab="",
                 main=(dimnames(x$importance)[[2]])[i],
                 xlim=c(0, max(x$importance[,i])), ...)
      }
    }
    mtext(outer=TRUE, side=3, text=main, cex=1.1)
  }
  invisible(x$importance)
}
