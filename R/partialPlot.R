partialPlot <- function(x, ...) UseMethod("partialPlot")

partialPlot.default <- function(x, ...)
    stop("partial dependence plot not implemented for this class of objects.\n")

partialPlot.randomForest <-
    function (x, pred.data, x.var, which.class, add = FALSE,
              n.pt = min(length(unique(pred.data[, xname])), 51), rug = TRUE,
              xlab=deparse(substitute(x.var)), ylab="",
              main=paste("Partial Dependence on", deparse(substitute(x.var))),
              ...) 
{
    classRF <- x$type != "regression"
    if (is.null(x$forest)) 
        stop("The randomForest object must contain the forest.\n")
    xname <- if (is.name(substitute(x.var))) deparse(substitute(x.var)) else eval(substitute(x.var))
    xv <- pred.data[, xname]
    n <- nrow(pred.data)
    if (classRF) {
        if (missing(which.class)) {
            focus <- 1
        }
        else {
            focus <- charmatch(which.class, colnames(x$votes))
            if (is.na(focus)) 
                stop(which.class, "is not one of the class labels.")
        }
    }
    if (is.factor(xv) && !is.ordered(xv)) {
        x.pt <- levels(xv)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- pred.data
            x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
            if (classRF) {
                pr <- predict(x, x.data, type = "prob")
                y.pt[i] <- mean(log(ifelse(pr[, focus] > 0,
                                           pr[, focus], 1)) -
                                rowMeans(log(ifelse(pr > 0, pr, 1))))
            } else {
                y.pt[i] <- mean(predict(x, x.data))
            }
        }
        if (add) {
            points(1:length(x.pt), y.pt, type = "h", lwd = 2, 
                   ...)
        } else {
            barplot(y.pt, width=rep(1, length(y.pt)), col="blue", xlab = xlab, 
                    ylab = ylab, main=main, names.arg=x.pt, ...)
        }
    } else {
        if (is.ordered(xv)) 
            xv <- as.numeric(xv)
        x.pt <- seq(min(xv), max(xv), length = n.pt)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- pred.data
            x.data[, xname] <- rep(x.pt[i], n)
            if (classRF) {
                pr <- predict(x, x.data, type = "prob")
                y.pt[i] <- mean(log(ifelse(pr[, focus] == 0, 1, pr[, focus]))
                                - rowMeans(log(ifelse(pr == 0, 1, pr))))
            } else {
                y.pt[i] <- mean(predict(x, x.data))
            }
        }
        if (add) {
            lines(x.pt, y.pt, ...)
        } else {
            plot(x.pt, y.pt, type = "l", 
                 xlab=xlab, ylab=ylab, main = main, ...)
        }
        if (rug) {
            if (n.pt > 10) {
                rug(quantile(xv, seq(0.1, 0.9, by = 0.1)), side = 1)
            } else {
                rug(unique(xv, side = 1))
            }
        }
    }
    invisible(list(x = x.pt, y = y.pt))
}
