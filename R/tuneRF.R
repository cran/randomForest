tuneRF <- function(x, y, mtryStart=if(is.factor(y)) floor(sqrt(ncol(x))) else
                   floor(ncol(x)/3), ntreeTry=50, stepFactor=2,
                   improve=0.05, trace=TRUE, plot=TRUE, doBest=FALSE, ...) {
  if (improve < 0) stop ("improve must be non-negative.")
  classRF <- is.factor(y)
  errorOld <- if (classRF) {
    randomForest(x, y, mtry=mtryStart, ntree=ntreeTry,
                 keep.forest=FALSE, ...)$err.rate[ntreeTry,1]
  } else {
    randomForest(x, y, mtry=mtryStart, ntree=ntreeTry,
                 keep.forest=FALSE, ...)$mse[ntreeTry]
  }
  if (trace) {
    cat("mtry =", mtryStart, " OOB error =",
        if (classRF) paste(100*round(errorOld, 4), "%", sep="") else
        errorOld, "\n")
  }

  oobError <- list()
  oobError[[as.character(mtryStart)]] <- errorOld  
  
  for (direction in c("left", "right")) {
    if (trace) cat("Searching", direction, "...\n")
    Improve <- 1.1*improve
    mtryBest <- mtryStart
    mtryCur <- mtryStart
    while (Improve >= improve) {
      mtryOld <- mtryCur
      mtryCur <- if (direction == "left") {
        max(1, ceiling(mtryCur / stepFactor))
      } else {
        min(ncol(x), floor(mtryCur * stepFactor))
      }
      if (mtryCur == mtryOld) break
      errorCur <- if (classRF) {
        randomForest(x, y, mtry=mtryCur, ntree=ntreeTry,
                     keep.forest=FALSE, ...)$err.rate[ntreeTry,"OOB"]
      } else {
        randomForest(x, y, mtry=mtryCur, ntree=ntreeTry,
                     keep.forest=FALSE, ...)$mse[ntreeTry]
      }
      if (trace) {
        cat("mtry =",mtryCur, "\tOOB error =",
            if (classRF) paste(100*round(errorCur, 4), "%", sep="") else
            errorCur, "\n")
      }
      oobError[[as.character(mtryCur)]] <- errorCur
      Improve <- 1 - errorCur/errorOld
      cat(Improve, improve, "\n")
      if (Improve > improve) {
        errorOld <- errorCur
        mtryBest <- mtryCur
      }
    }
  }
  res <- unlist(oobError)[order(as.numeric(names(oobError)))]
  res <- cbind(mtry=as.numeric(names(res)), OOBError=res)

  if (plot) {
    plot(res, xlab=expression(m[try]), ylab="OOB Error", type="o", log="x",
         xaxt="n")
    axis(1, at=res[,"mtry"])
  }

  if (doBest) 
    res <- randomForest(x, y, mtry=res[which.min(res[,2]), 1], ...)
  
  res
}
