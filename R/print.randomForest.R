"print.randomForest" <-
function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n")
  cat("               Type of random forest: ", x$type, "\n", sep="")
  cat("                     Number of trees: ", x$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", x$mtry, "\n\n", sep="")
  if(x$type == "classification") {
    if(!is.null(x$confusion)) {
      cat("        OOB estimate of  error rate: ",
          round(x$err.rate[length(x$err.rate)]*100, dig=2), "%\n", sep="")
      cat("Confusion matrix:\n")
      print(x$confusion)
      if(!is.null(x$test)) {
        cat("                Test set error rate: ",
            round(x$test$err.rate[length(x$test$err.rate)]*100, dig=2), "%\n",
            sep="")
        cat("Confusion matrix:\n")
        print(x$test$confusion)
      }
    }
  }
  if(x$type == "regression") {
    if(!is.null(x$mse)) {
      cat("          Mean of squared residuals: ", x$mse[length(x$mse)],
          "\n", sep="")
      cat("                    % Var explained: ",
          round(100*x$rsq[length(x$rsq)], dig=2), "\n", sep="")
      if(!is.null(x$test)) {
        cat("                      Test set rMSE: ",
            round(x$test$mse[length(x$test$mse)]*100, dig=2), "%\n", sep="")
        cat("                    % Var explained: ",
            round(100*x$test$rsq[length(x$test$rsq)], dig=2), "\n", sep="")
      }
    }
  }
}
