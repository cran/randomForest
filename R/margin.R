margin <- function(rf, observed) {
    if( !inherits(rf, "randomForest") ) {
        stop("margin defined for Random Forests")
    }
    if( is.null(rf$votes) ) {
        stop("margin is only defined if votes are present")
    }
    if( !is.factor(observed) ) {
        stop(deparse(substitute(observed)), " is not a factor")
    }
    augD <- rf$votes
    if( any(augD > 1) ) {
        augD <- sweep(augD, 1, rowSums(augD), "/")
    }
    augD <- data.frame(augD, observed)
    names(augD) <- c(dimnames(rf$votes)[[2]], "observed")
    nlev <- length(levels(observed))
    
    ans<- apply(augD, 1, function(x) { pos <- match(x[nlev+1], names(x));
                                       t1 <- as.numeric(x[pos]);
                                       t2 <- max(as.numeric(x[-c(pos, nlev+1)]));
                                       t1 - t2 }
                )
    names(ans) <- observed
    class(ans) <- "margin"
    ans
}

plot.margin <- function(x, sort=TRUE, ...) {
    if (sort) x <- sort(x)
    nF <- factor(names(x))
    nlevs <- length(levels(nF))
    if ( require(RColorBrewer) && nlevs < 12) {
        pal <- brewer.pal(nlevs,"Set1")
    } else {
        pal <- rainbow(nlevs)
    }
    plot.default(x, col=pal[as.numeric(nF)], pch=20, ... )
}
