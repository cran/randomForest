partial.plot <- function(x, ...) UseMethod("partial.plot")

partial.plot.default <- function(x, ...)
  stop("partial dependence plot not implemented for this class of objects.\n")

partial.plot.randomForest <-
  function(x, pred.data, x.var, which.class, add=FALSE,
           n.pt=min(length(unique(pred.data[,deparse(substitute(x.var))])),
             51), rug=TRUE, ...) {
    classRF <- x$type != "regression"
    if(is.null(x$forest))
      stop("The randomForest object must contain the forest.\n")
    xname <- deparse(substitute(x.var))
    xv <- pred.data[, xname]
    n <- nrow(pred.data)
    if(classRF) {
      if(missing(which.class)) {
        focus <- 1
      } else {
        focus <- charmatch(which.class, colnames(x$votes))
        if(is.na(focus)) stop(which.class, "is not one of the class labels.")
      }
    }
    
    if(is.factor(xv)) {
      x.pt <- levels(xv)
      y.pt <- numeric(length(x.pt))
      for(i in seq(along=x.pt)) {
        x.data <- pred.data
        x.data[,xname] <- factor(rep(x.pt[i], n), levels = x.pt)
        if(classRF) {
##          pr <- predict(x, x.data, type="prob")
##          y.pt[i] <- mean(log(pr[,focus]) - rowMeans(log(pr)))
          y.pt[i] <- mean(predict(x, x.data, type="prob")[,focus])
        } else {
          y.pt[i] <- mean(predict(x, x.data))
        }
      }
      if(add) {
        points(1:length(x.pt),  y.pt, type="h", lwd=2, ...)
      } else {
        plot(1:length(x.pt), y.pt, type="h", lwd=2, xlab=xname, ylab="",
             main = paste("Partial Dependence of", xname))
      }
    } else {
      x.pt <- seq(min(xv), max(xv), length=n.pt)
      y.pt <- numeric(length(x.pt))
      for(i in seq(along=x.pt)) {
        x.data <- pred.data
        x.data[, xname] <- rep(x.pt[i], n)
        if(classRF) {
##          pr <- predict(x, x.data, type="prob")
##          y.pt[i] <- mean(log(pr[,focus]) - rowMeans(log(pr)))
          y.pt[i] <- mean(predict(x, x.data, type="prob")[,focus])
        } else {
          y.pt[i] <- mean(predict(x, x.data))
        }
      }

      if(add) {
        lines(x.pt, y.pt, ...)
      } else {
        plot(x.pt, y.pt, xlab=xname, ylab="", type="l",
             main=paste("Partial Dependence of", xname))
      }
      if(rug & n.pt > 10) rug(quantile(xv, seq(.1, .9, by=.1)), side=1)
    }
    invisible(list(x=x.pt, y=y.pt))
  }
