"predict.randomForest" <-
  function (object, newdata, type = "response", norm.votes = TRUE, 
            proximity = FALSE, ...) 
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
    if(proximity & is.null(object$proximity))
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
  mdim <- ncol(x)
  ntest <- nrow(x)
  ntree <- object$forest$ntree
  maxcat <- object$forest$maxcat
  nclass <- object$forest$nclass
  nrnodes <- object$forest$nrnodes
  show.error <- 0
  x <- t(data.matrix(x))
  if(proximity) {
    proxmatrix <- matrix(0, ntest, ntest)
  } else {
    proxmatrix <- numeric(1)
  }
  if(object$type == "regression") {
    ans <- .C("runrforest",
              as.double(x),
              ypred = as.double(numeric(ntest)),
              as.integer(mdim),
              as.integer(ntest),
              as.integer(ntree),
              as.integer(object$forest$ndbigtree),
              as.integer(object$forest$treemap),
              as.integer(object$forest$nodestatus),
              as.integer(object$forest$nrnodes),
              as.double(object$forest$xbestsplit),
              as.double(object$forest$avnode),
              as.integer(object$forest$bestvar),
              as.integer(object$forest$ncat),
              as.integer(proximity),
              proximity = as.double(proxmatrix),
              DUP=FALSE,
              PACKAGE = "randomForest")[c(2, 15)]
    if (!is.null(object$coefs)) {
      yhat <- ans$ypred + (object$coefs[1] + object$coefs[2] * ans$ypred)
    } else {
      yhat <- ans$ypred
    }
    if (!proximity) {
      res <- yhat
    } else {
      res <- list(pred = yhat, proximity = structure(ans$proximity,
                                     dim=c(ntest, ntest),
                                     dimnames=list(rn, rn)))
    }
  } else {
    t1 <- .Fortran("runforest",
                   mdim = as.integer(mdim),
                   ntest = as.integer(ntest), 
                   nclass = as.integer(object$forest$nclass),
                   maxcat = as.integer(maxcat), 
                   nrnodes = as.integer(nrnodes),
                   labelsts = as.integer(as.numeric(show.error)), 
                   jbt = as.integer(ntree),
                   clts = as.integer(numeric(ntest)), 
                   xts = as.double(x),
                   xbestsplit = as.double(object$forest$xbestsplit), 
                   pid = as.double(object$forest$pid),
                   cutoff = as.double(object$forest$cutoff),
                   countts = as.double(numeric(nclass * ntest)),
                   treemap = as.integer(aperm(object$forest$treemap, 
                     c(2, 1, 3))),
                   nodestatus = as.integer(object$forest$nodestatus), 
                   cat = as.integer(object$forest$ncat),
                   cbestsplit = as.integer(numeric(maxcat * nrnodes)),
                   nodeclass = as.integer(object$forest$nodeclass), 
                   jts = as.integer(numeric(ntest)),
                   jet = as.integer(numeric(ntest)), 
                   bestvar = as.integer(object$forest$bestvar),
                   nodexts = as.integer(numeric(ntest)), 
                   ndbigtree = as.integer(object$forest$ndbigtree), 
                   prox = as.integer(as.numeric(proximity)),
                   proxmatrix = as.double(proxmatrix),
                   DUP=FALSE,
                   PACKAGE = "randomForest")
    out.class.votes <- t(matrix(t1$countts, nr = nclass, nc = ntest))
    if (norm.votes) 
      out.class.votes <- t(apply(out.class.votes, 1, function(x) x/sum(x)))
    z <- matrix(NA, ntest, nclass, dimnames = list(rn, levels(object$predicted)))
    z[keep, ] <- out.class.votes
    res <- z
    if (out.type == 1) {
      out.class <- max.col(out.class.votes)
      out.class <- if (ncol(z) > 1) 
        levels(object$predicted)[max.col(z)]
      else levels(object$predicted)[1 + (z > 0.5)]
      out.class <- factor(out.class, levels = levels(object$predicted))
      res <- out.class
    }
    if(proximity)
      res <- list(pred = res, proximity = structure(t1$proxmatrix,
                                dim = c(ntest, ntest),
                                dimnames = list(rn, rn)))
  }
  res
}
