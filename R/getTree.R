getTree <- function(rfobj, k=1) {
  if (is.null(rfobj$forest)) {
    stop("No forest component in ", deparse(substitute(rfobj)))
  }
  if (k > rfobj$ntree) {
    stop("There are fewer than ", k, "trees in the forest")
  }
  if (rfobj$type == "regression") {
      tree <- cbind(rfobj$forest$leftDaughter[,k],
                    rfobj$forest$rightDaughter[,k],
                    rfobj$forest$bestvar[,k],
                    rfobj$forest$xbestsplit[,k],
                    rfobj$forest$nodestatus[,k],
                    rfobj$forest$nodepred[,k])[1:rfobj$forest$ndbigtree[k],]
  } else {
      tree <- cbind(rfobj$forest$treemap[,,k],
                    rfobj$forest$rightDaughter[,k],
                    rfobj$forest$bestvar[,k],
                    rfobj$forest$xbestsplit[,k],
                    rfobj$forest$nodestatus[,k],
                    rfobj$forest$nodepred[,k])[1:rfobj$forest$ndbigtree[k],]
  }      
  colnames(tree) <- c("left daughter", "right daughter", "split var",
                      "split point", "status", "prediction")
  tree
}

