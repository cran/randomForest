var.imp.plot <- function(x, sort=TRUE,
                         n.var=min(30, if(is.null(dim(x$importance)))
                           length(x$importance) else nrow(x$importance)),
                         class = NULL,
                         main=deparse(substitute(x)), ...) {
  if (!inherits(x, "randomForest"))
    stop("This function only works for objects of class `randomForest'")
  if (is.null(dim(x$importance)) || !is.null(class)) {
    if (is.null(class)) {
      imp <- x$importance
    } else {
      imp <- x$importance[,class]
    }
    if (sort) {
      ord <- order(imp, decreasing=TRUE)[1:n.var]
      imp <- imp[ord]
      dotchart(rev(imp), xlab="Importance", pch=16,
           ylab="", main=main, ...)
    } else {
      dotchart(imp, xlab="Importance", pch=16,
           ylab="", main=main, ...)
    }
  } else {
    ## Extract the last two columns.
    impmat <- x$importance[,(ncol(x$importance)-1):ncol(x$importance)]
    imp <- vector(2, mode="list")
    names(imp) <- colnames(impmat)
    op <- par(mfrow=c(1, 2), mar=c(4, 5, 4, 1), mgp=c(2, .8, 0),
              oma=c(0,0,2,0))
    on.exit(par(op))
  
    for (i in 1:2) {
      if(sort) {
        ord <- order(impmat[,i], decreasing=TRUE)[1:n.var]
        imp[[i]] <- impmat[ord, i]
        maximp <- max(imp[[i]])
        dotchart(rev(imp[[i]]), pch=16, xlab="Importance",
                 ylab="", main=names(imp)[i], xlim=c(0, maximp),
                 ...)
      } else {
        imp[[i]] <- impmat[, i]
        dotchart(imp[[i]], pch=16, xlab="Importance", ylab="",
                 main=colnames(imp)[i], ...)
      }
    }
    mtext(outer=TRUE, side=3, text=main, cex=1.1)
  }
  invisible(imp)
}
