var.imp.plot <- function(x, sort=TRUE, label=TRUE,
                         n.var=min(50, if(is.null(dim(x$importance)))
                           length(x$importance) else nrow(x$importance)),
                         main=deparse(substitute(x)), ...) {
  if(!inherits(x, "randomForest"))
    stop("This function only works for objects of class `randomForest'")
  if(is.null(dim(x$importance))) {
    if(sort) {
      ord <- order(x$importance, decreasing=TRUE)[1:n.var]
      maximp <- max(x$importance)
      if(label) maximp <- maximp*1.08  # make room for labels
      plot(x$importance[ord], ylim=c(0, maximp), ylab="Importance",
           xlab="", axes=FALSE, type="h", main=main, ...)
      axis(2)
      box(bty="l")
      if(label) {
        text(1:n.var, x$importance[ord] + 0.05*maximp,
           names(x$importance[ord]), cex=.8)
      }
    } else {
      plot(x$importance, ylim=c(0, max(x$importance)), ylab="Importance",
           xlab="", axes=FALSE, type="h", main=main, ...)
      axis(2)
      box(bty="l")
    }
  } else {
    n.meas <- ncol(x$importance)
    op <- par(mfrow=c(n.meas/2, 2), mar=c(2, 4, 4, 1), mgp=c(2, .8, 0),
              oma=c(0,0,2,0))
    on.exit(par(op))
  
    for(i in 1:n.meas) {
      if(sort) {
        ord <- order(x$importance[,i], decreasing=TRUE)[1:n.var]
        maximp <- max(x$importance[,i])
        if(label) maximp <- maximp*1.07  # make room for labels
##        print(x$importance[ord,i])
        plot(seq(along=ord), x$importance[ord, i], ylim=c(0, maximp),
             ylab="Importance", xlab="", axes=FALSE, type="h",
             main=colnames(x$importance)[i], ...)
        axis(2)
        box(bty="l")
        if(label) {
          text(1:n.var, x$importance[ord, i] + 0.04*maximp,
               rownames(x$importance)[ord], cex=.8)
        }
      } else {
        plot(x$importance[,i], ylim=c(0, max(x$importance[,i])),
             ylab="Importance", xlab="", axes=FALSE, type="h",
             main=(dimnames(x$importance)[[2]])[i], ...)
        axis(2)
        box(bty="l")
      }
    }
    mtext(outer=TRUE, side=3, text=main, cex=1.05)
  }
  invisible(x$importance)
}

