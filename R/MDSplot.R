MDSplot <- function(rf, fac, k=2, palette=NULL, pch=20, ...) {
    if (!inherits(rf, "randomForest")) 
        stop(deparse(substitute(rf)), " must be a randomForest object")
    if(is.null(rf$proximity)) 
        stop(deparse(substitute(rf)), " does not contain a proximity matrix")
    op <- par(pty="s")
    on.exit(par(op))
    Rver <- as.numeric(substring(paste(R.version[c("major", "minor")],
                                       collapse="."), 1, 3))
    cmdscale <- if (Rver >= 1.9) stats:::cmdscale else mva:::cmdscale
    rf.mds <- cmdscale(1 - rf$proximity, eig=TRUE, k=k)
    colnames(rf.mds$points) <- paste("Dim", 1:k)
    nlevs <- nlevels(fac)
    if (is.null(palette)) {
        palette <- if (require(RColorBrewer) && nlevs < 12)
            brewer.pal(nlevs, "Set1") else rainbow(nlevs)
    }
    if (k <= 2) {
        plot(rf.mds$points, col=palette[as.numeric(fac)], pch=pch, ...)
    } else {
        pairs(rf.mds$points, col=palette[as.numeric(fac)], pch=pch, ...)
    }
    invisible(rf.mds)
}
