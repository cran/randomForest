"predict.randomForest" <-
  function (object, newdata, type = "response", norm.votes = TRUE,
            predict.all=FALSE, proximity = FALSE, nodes=FALSE, ...) 
{
  if (!inherits(object, "randomForest")) 
    stop("object not of class randomForest")
  if (is.null(object$forest)) stop("No forest component in the object")
  out.type <- charmatch(tolower(type),
                        c("response", "prob", "vote", "class"))
  if (is.na(out.type)) 
    stop("type must be one of 'response', 'prob', 'vote'")
  if (out.type == 4) out.type <- 1
  if (out.type != 1 && object$type == "regression")
    error("'prob' or 'vote' not meaningful for regression")
  if (out.type == 2) 
    norm.votes <- TRUE
  if (missing(newdata)) {
    if (object$type == "regression") return(object$predicted)
    if (proximity & is.null(object$proximity))
      warning("cannot return proximity without new data if random forest object does not already have proximity")
    if (out.type == 1) {
      if (proximity) {
        return(list(pred = object$predicted,
                    proximity = object$proximity))
      } else return(object$predicted)
    }
    if (norm.votes) { 
      t1 <- t(apply(object$votes, 1, function(x) { x/sum(x) }))
      if(proximity) return(list(pred = t1, proximity = object$proximity))
      else return(t1)
    } else {
      if(proximity) return(list(pred = object$votes, proximity = object$proximity))
      else return(object$votes)
    }
  }

  if (object$type == "unsupervised") 
    stop("Can't predict unsupervised forest.")

  if (inherits(object, "randomForest.formula")) {
    newdata <- as.data.frame(newdata)
    rn <- row.names(newdata)
    Terms <- delete.response(object$terms)
    x <- model.frame(Terms, newdata, na.action = na.omit)
    keep <- match(row.names(x), rn)
  } else {
    if (is.null(dim(newdata))) 
      dim(newdata) <- c(1, length(newdata))
    x <- newdata
    if (nrow(x) == 0) 
      stop("newdata has 0 rows")
    if (any(is.na(x))) 
      stop("missing values in newdata")
    keep <- 1:nrow(x)
    rn <- rownames(x)
  }
  if (is.data.frame(x)) {
    for(i in seq(along=ncol(x))) {
      if(is.ordered(x[[i]])) x[[i]] <- as.numeric(x[[i]])
    }
    cat.new <- sapply(x, function(x) if (is.factor(x) && !is.ordered(x)) 
                      length(levels(x))
    else 1)
    if (!all(object$forest$ncat == cat.new)) 
      stop("Type of predictors in new data do not match that of the training data.")
  }
  vname <- if (is.null(dim(object$importance))) {
      names(object$importance)
  } else {
      rownames(object$importance)
  }
  if (any(colnames(x) != vname)) stop("names of predictor variables do not match")
  mdim <- ncol(x)
  ntest <- nrow(x)
  ntree <- object$forest$ntree
  maxcat <- object$forest$maxcat
  nclass <- object$forest$nclass
  nrnodes <- object$forest$nrnodes
  show.error <- 0
  x <- t(data.matrix(x))

  if (predict.all) {
    if (object$type == "regression") {
      treepred <- double(ntest * ntree)
    } else {
      treepred <- integer(ntest * ntree)
    }
  } else {
    treepred <- numeric(ntest)
  }
  
  if (proximity) {
    proxmatrix <- matrix(0, ntest, ntest)
  } else {
    proxmatrix <- numeric(1)
  }

  nodexts <- if (nodes) integer(ntest*ntree) else integer(ntest)
  
  if(object$type == "regression") {
    keepIndex <- "ypred"
    if (predict.all) keepIndex <- c(keepIndex, "treepred")
    if (proximity) keepIndex <- c(keepIndex, "proximity")
    ans <- .C("runrforest",
              as.double(x),
              ypred = double(ntest),
              as.integer(mdim),
              as.integer(ntest),
              as.integer(ntree),
              as.integer(object$forest$ndbigtree),
              as.integer(aperm(object$forest$treemap, c(2, 1, 3))),
              as.integer(object$forest$nodestatus),
              as.integer(object$forest$nrnodes),
              as.double(object$forest$xbestsplit),
              as.double(object$forest$nodepred),
              as.integer(object$forest$bestvar),
              as.integer(object$forest$ncat),
              as.integer(predict.all),
              treepred = as.double(treepred),
              as.integer(proximity),
              proximity = as.double(proxmatrix),
              DUP=FALSE,
              PACKAGE = "randomForest")[keepIndex]
    ## Apply bias correction if needed.
    if (!is.null(object$coefs)) {
      yhat <- object$coefs[1] + object$coefs[2] * ans$ypred
    } else {
      yhat <- ans$ypred
    }
    if (predict.all) {
      treepred <- matrix(ans$treepred, length(keep),
                         dimnames=list(rn[keep], NULL))
    }
    if (!proximity) {
      res <- if (predict.all)
        list(aggregate=yhat, individual=treepred) else yhat
    } else {
      res <- list(predicted = yhat, proximity = structure(ans$proximity,
                                     dim=c(ntest, ntest),
                                     dimnames=list(rn, rn)))
    }
  } else {
    countts <- matrix(0, ntest, nclass)
    t1 <- .C("runforest",
             mdim = as.integer(mdim),
             ntest = as.integer(ntest), 
             nclass = as.integer(object$forest$nclass),
             maxcat = as.integer(maxcat), 
             nrnodes = as.integer(nrnodes),
             jbt = as.integer(ntree),
             xts = as.double(x),
             xbestsplit = as.double(object$forest$xbestsplit), 
             pid = as.double(object$forest$pid),
             cutoff = as.double(object$forest$cutoff),
             countts = as.double(countts),
             treemap = as.integer(aperm(object$forest$treemap, 
               c(2, 1, 3))),
             nodestatus = as.integer(object$forest$nodestatus), 
             cat = as.integer(object$forest$ncat),
             cbestsplit = as.integer(numeric(maxcat * nrnodes)),
             nodepred = as.integer(object$forest$nodepred), 
             treepred = as.integer(treepred),
             jet = as.integer(numeric(ntest)), 
             bestvar = as.integer(object$forest$bestvar),
             nodexts = nodexts,
             ndbigtree = as.integer(object$forest$ndbigtree), 
             predict.all = as.integer(predict.all),
             prox = as.integer(proximity),
             proxmatrix = as.double(proxmatrix),
             nodes = as.integer(nodes),
             DUP=TRUE,
             PACKAGE = "randomForest")
    if (out.type > 1) {
      out.class.votes <- t(matrix(t1$countts, nr = nclass, nc = ntest))
      if (norm.votes) 
        out.class.votes <-
          sweep(out.class.votes, 1, rowSums(out.class.votes), "/")
      z <- matrix(NA, ntest, nclass, dimnames = list(rn, levels(object$predicted)))
      z[keep, ] <- out.class.votes
      res <- z
    } else {
      out.class <- factor(rep(NA, length(rn)),
                          levels=1:length(object$classes),
                          labels=object$classes)
      out.class[keep] <- object$classes[t1$jet]
      names(out.class[keep]) <- rn[keep]
      res <- out.class
    }
    if (predict.all) {
      treepred <- matrix(object$classes[t1$treepred],
                         nrow=length(keep), dimnames=list(rn[keep], NULL))
      res <- list(aggregate=res, individual=treepred)
    }
    if(proximity)
      res <- list(predicted = res, proximity = structure(t1$proxmatrix,
                                dim = c(ntest, ntest),
                                dimnames = list(rn[keep], rn[keep])))
    if (nodes) attr(res, "nodes") <- matrix(t1$nodexts, ntest, ntree,
                                            dimnames=list(rn[keep], 1:ntree))
  }
  res
}
