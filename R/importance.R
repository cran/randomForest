importance <- function(x, ...)  UseMethod("importance")

importance.default <- function(x, ...)
  stop("No method implemented for this class of object")

importance.randomForest <- function(x, type=2, class, ...) {
  if (!inherits(x, "randomForest"))
    stop("x is not of class randomForest")
  if (!(type %in% 1:2)) stop("Wrong type specified")
  hasImp <- !is.null(dim(x$importance))
  if (type == 2) {
    if (hasImp) {
      return(x$importance[,ncol(x$importance)])
    } else {
      return(x$importance)
    }
  } else {
    if (!hasImp) stop("that measure was not computed")
    imp <- if (missing(class)) {
      x$importance[,ncol(x$importance)-1]
    } else {
      x$importance[,class]
    }
    names(imp) <- rownames(x$importance)
    return(imp)
  }
}
