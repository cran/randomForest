combine <- function(...) {
  pad0 <- function(x, len) c(x, rep(0, len-length(x)))
  padm0 <- function(x, len) rbind(x, matrix(0, nrow=len-nrow(x), ncol=ncol(x)))
  rflist <- list(...)
  are.forest <- sapply(rflist, function(x) inherits(x, "randomForest")) 
  if(any(!are.forest)) stop("Argument must be a list of randomForest objects")
  ## Use the first component as a template
  rf <- rflist[[1]]
  classRF <- rf$type == "classification"
  trees <- sapply(rflist, function(x) x$ntree)
  ntree <- sum(trees)
  rf$ntree <- ntree
  nforest <- length(rflist)
  have.test <- !any(sapply(rflist, function(x) is.null(x$test)))
  
  ## Combine the forest component, if any
  have.forest <- sapply(rflist, function(x) !is.null(x$forest))
  if(all(have.forest)) {
    nrnodes <- max(sapply(rflist, function(x) x$forest$nrnodes))
    rf$forest$nrnodes <- nrnodes
    rf$forest$ndbigtree <-
      unlist(sapply(rflist, function(x) x$forest$ndbigtree))
    rf$forest$nodestatus <-
      do.call("cbind", lapply(rflist, function(x) padm0(x$forest$nodestatus,
                                                        nrnodes)))
    rf$forest$bestvar <-
      do.call("cbind", lapply(rflist, function(x) padm0(x$forest$bestvar,
                                                        nrnodes)))
    rf$forest$xbestsplit <-
      do.call("cbind", lapply(rflist, function(x) padm0(x$forest$xbestsplit,
                                                        nrnodes)))
    if(classRF) {
      rf$forest$nodeclass <-
        do.call("cbind", lapply(rflist, function(x) padm0(x$forest$nodeclass,
                                                         nrnodes)))
    } else {
      rf$forest$avnode <-
        do.call("cbind", lapply(rflist, function(x) padm0(x$forest$avnode,
                                                          nrnodes)))
    }
    tree.dim <- dim(rf$forest$treemap)
    rf$forest$treemap <-
      array(unlist(lapply(rflist, function(x) apply(x$forest$treemap, 2:3,
                                                    pad0, nrnodes))),
            c(nrnodes, 2, ntree))
    rf$forest$ntree <- ntree
  } else {
    rf$forest <- NULL
  }

  if(classRF) {
    ## Combine the votes matrix: if raw counts, just add.  Otherwise average.  
    rf$votes <- 0
    if(any(rf$votes > 1)) {
      for(i in 1:nforest) rf$votes <- rf$votes + rflist[[i]]$votes
    } else {
      warning("randomForest created with norm.votes=TRUE.  Result is averaged (weighted by number of trees), and thus not the most accurate.")
      for(i in 1:nforest)
        rf$votes <- rf$votes + rflist[[i]]$votes * rflist[[i]]$ntree
      rf$votes <- rf$votes / ntree    
    }
    rf$predicted <- factor(colnames(rf$votes)[max.col(rf$votes)],
                           levels=levels(rf$predicted))
    if(have.test) {
      rf$test$votes <- 0
      if(any(rf$test$votes > 1)) {
        for(i in 1:nforest)
          rf$test$votes <- rf$test$votes + rflist[[i]]$test$votes
      } else {
        for(i in 1:nforest)
          rf$test$votes <- rf$test$votes +
            rflist[[i]]$test$votes * rflist[[i]]$ntree
      }
      rf$test$predicted <-
        factor(colnames(rf$test$votes)[max.col(rf$test$votes)],
               levels=levels(rf$test$predicted))
    }
  } else {
    rf$predicted <- 0
    for(i in 1:nforest) rf$predicted <- rf$predicted +
      rflist[[i]]$predicted * rflist[[i]]$ntree
    rf$predicted <- rf$predicted / ntree
    if(have.test) {
      rf$test$predicted <- 0
      for(i in 1:nforest) rf$test$predicted <- rf$test$predicted +
        rflist[[i]]$test$predicted * rflist[[i]]$ntree
      rf$test$predicted <- rf$test$predicted / ntree
    }
  }

  ## If variable importance is in all of them, compute the average
  ## (weighted by the number of trees in each forest)
  have.imp <- !any(sapply(rflist, function(x) is.null(x$importance)))
  if(have.imp) {
    rf$importance <- 0
    for(i in 1:nforest)
      rf$importance <- rf$importance + rflist[[i]]$importance * rflist[[i]]$ntree
    rf$importance <- rf$importance / ntree
  }

  ## If proximity is in all of them, compute the average
  ## (weighted by the number of trees in each forest)
  have.prox <- !any(sapply(rflist, function(x) is.null(x$proximity)))
  if(have.prox) {
    rf$proximity <- 0
    for(i in 1:nforest)
      rf$proximity <- rf$proximity + rflist[[i]]$proximity * rflist[[i]]$ntree
    rf$proximity <- rf$proximity / ntree
  }

  ## Set confusion matrix and error rates to NULL
  if(classRF) {
    rf$confusion <- NULL
    rf$err.rate <- NULL
    if(have.test) {
      rf$test$confusion <- NULL
      rf$err.rate <- NULL
    }
  } else {
    rf$mse <- rf$rsq <- NULL
    if(have.test) rf$test$mse <- rf$test$rsq <- NULL
  }
  
  rf
}
