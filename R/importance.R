importance <- function(x, ...)  UseMethod("importance")

importance.default <- function(x, ...)
    stop("No method implemented for this class of object")

importance.randomForest <- function(x, type=NULL, class=NULL, scale=TRUE,
                                    ...) {
    if (!inherits(x, "randomForest"))
        stop("x is not of class randomForest")
    classRF <- x$type != "regression"
    hasImp <- !is.null(dim(x$importance))
    hasType <- !is.null(type)
    allImp <- is.null(type) && hasImp
    if (hasType) {
        if (!(type %in% 1:2)) stop("Wrong type specified")
        if (type == 1 && !hasImp) stop("That measure was not computed")
        if (type == 2 && !is.null(class))
            stop("No class-specific measure for that type")
    }
    
    imp <- x$importance
    if (hasType && type == 2) {
        if (hasImp) imp <- imp[, ncol(imp)]
    } else {
        if (hasImp) {
            if (scale) {
                SD <- x$importanceSD
                imp[,-ncol(imp)] <-
                    imp[,-ncol(imp)] / ifelse(SD < .Machine$double.eps, 1, SD)
            }
            if (!allImp) {
                if (is.null(class)) {
                    ## The average decrease in accuracy measure:
                    imp <- imp[, ncol(imp)-1]
                } else {
                    if (classRF) {
                        whichCol <- match(class, colnames(x$importance))
                    } else {
                        whichCol <- 1
                    }
                    imp <- imp[, whichCol]
                }
                names(imp) <- rownames(x$importance)
            }
        }
    }
    imp
}
