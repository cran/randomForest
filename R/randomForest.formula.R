"randomForest.formula" <-
function(formula, data = NULL, subset, ...){	
### formula interface for randomForest.
### code gratefully stolen from svm.formula (package e1071).
###
    call <- match.call()
    if (!inherits(formula, "formula"))
        stop("method is only for formula objects")
    m <- match.call(expand = FALSE)
    if (is.matrix(eval(m$data, sys.frame(sys.parent()))))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.frame(sys.parent()))
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 0
    x <- model.matrix(Terms, m)
    y <- model.extract(m, response)
    ret <- randomForest.default(x, y, ...)
    ret[["call"]] <- call
    ret$terms <- Terms
    class(ret) <- c("randomForest.formula", "randomForest")
    return(ret)
}
