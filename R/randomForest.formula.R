"randomForest.formula" <-
    function(x, data = NULL, ..., subset, na.action = na.fail) {	
### formula interface for randomForest.
### code gratefully stolen from svm.formula (package e1071).
###
    if (!inherits(x, "formula"))
        stop("method is only for formula objects")
    call <- match.call()
    m <- match.call(expand = FALSE)
    names(m)[2] <- "formula"
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 0
    y <- model.response(m)
    if(!is.null(y)) m <- m[, -1, drop=FALSE]
    for (i in seq(along=ncol(m))) {
        if(is.ordered(m[[i]])) m[[i]] <- as.numeric(m[[i]])
    }
    ret <- randomForest.default(m, y, ...)
    cl <- match.call()
    cl[[1]] <- as.name("randomForest")
    ret$call <- cl
    ret$terms <- Terms
    if (!is.null(attr(m, "na.action"))) 
        ret$na.action <- attr(m, "na.action")
    class(ret) <- c("randomForest.formula", "randomForest")
    return(ret)
}
