plot.randomForest <- function(x, type="l", main=deparse(substitute(x)), ...) {
  if(x$type == "unsupervised")
    stop("No plot for unsupervised randomForest.")
  test <- !is.null(x$test)
  if(x$type == "regression") {
    err <- x$mse
    if(test) err <- cbind(err, x$test$mse)
  } else {
    err <- x$err.rate
    if(test) err <- cbind(err, x$test$err.rate)
  }
  if(test) {
    colnames(err) <- c("OOB", "Test")
    matplot(1:x$ntree, err, type = type, xlab="trees", ylab="Error",
            main=main, ...)
  } else {
    plot(err, type = type, xlab="trees", ylab="Error", main=main, ...)
  }
  invisible(err)
}

  
