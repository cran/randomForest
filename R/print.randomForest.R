"print.randomForest" <-
function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n")
  cat("                  Type of random forest: ", x$type, "\n", sep="")
  cat("                        Number of trees: ", x$ntree, "\n",sep="")
  cat("Number of variables tried at each split: ", x$mtry, "\n\n", sep="")
  if(x$type == "classification") {
    cat("            OOB estimate of  error rate: ",
        round(x$err.rate*100, dig=2), "%\n", sep="")
    cat("Confusion matrix:\n")
    print(x$confusion)
  }
  if(x$type == "regression") {
    cat("             Mean of squared residuals: ", x$mse, "\n", sep="")
    cat("                       % Var explained: ", round(100*x$rsq, dig=2),
        "\n", sep="")
  }
}
