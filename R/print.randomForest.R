"print.randomForest" <-
function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n")
  cat("                        Number of trees: ", x$ntree, "\n",sep="")
  cat("Number of variables tried at each split: ", x$mtry, "\n\n", sep="")
  cat("Final error rate: ", round(x$err.rate*100, dig=2), "%\n", sep="")
  if(x$supervised) {
    cat("Confusion matrix:\n")
    print(x$confusion)
  }
}
