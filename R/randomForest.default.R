"randomForest.default" <-
function(x, y=NULL,  xtest=NULL, ytest=NULL, addclass=0, ntree=500,
         mtry=ifelse(!is.null(y) && !is.factor(y), max(floor(ncol(x)/3), 1),
           floor(sqrt(ncol(x)))), classwt=NULL, cutoff, sampsize,
         nodesize= ifelse(!is.null(y) && !is.factor(y), 5, 1),
         importance=FALSE, proximity=FALSE, oob.prox=TRUE, outscale=FALSE,
         norm.votes=TRUE, do.trace=FALSE, keep.forest=is.null(xtest),
         corr.bias=FALSE, ...)
{
  classRF <- is.null(y) || is.factor(y)
  if (!classRF && length(unique(y)) <= 5) {
    warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
  }
  n <- nrow(x)
  p <- ncol(x)
  if(n == 0)
    stop("data (x) has 0 rows")
  x.row.names <- rownames(x)
  if(is.null(colnames(x))) {
    x.col.names <- 1:ncol(x)
  } else {
    x.col.names <- colnames(x)
  }
  
  ## bypass R's lazy evaluation:
  keep.forest <- keep.forest

  if(!is.null(xtest)) {
    if(ncol(x) != ncol(xtest))
      stop("x and xtest must have same number of columns") 
    ntest <- nrow(xtest)
    xts.row.names <- rownames(xtest)
  }

  if(mtry > p) {
    mtry <- p
    warning("mtry can not be larger than number of predictors.  Reset to equal to number of predictors")
  }
  if(!is.null(y)) {
    if(length(y) != n) stop("length of response must be the same as predictors")
    addclass <- 0
  } else {
    if(addclass == 0) addclass <- 1
    y <- factor(c(rep(1, n), rep(2, n)))
    x <- rbind(x, x)
    keep.forest <- FALSE
  }
  if(!is.element(addclass,0:2))
    stop("addclass can only take on values 0, 1, or 2")
  
  if(any(is.na(x))) stop("NA not permitted in predictors")
  if(!is.null(xtest) && any(is.na(xtest))) stop("NA not permitted in xtest")
  if(any(is.na(y))) stop("NA not permitted in response")
  if(!is.null(ytest) && any(is.na(ytest)))
    stop("NA not permitted in ytest")

  if(is.data.frame(x)) {
    ncat <- sapply(x, function(x) if(is.factor(x) && !is.ordered(x))
                   length(levels(x)) else 1)
    x <- data.matrix(x)
    if(!is.null(xtest)) {
      if(!is.data.frame(xtest))
        stop("xtest must be data frame if x is")
      ncatts <- sapply(xtest, function(x) if(is.factor(x) && !is.ordered(x))
                       length(levels(x)) else 1)
      if(!all(ncat == ncatts))
        stop("columns of x and xtest must be the same type")
      xtest <- data.matrix(xtest)
    }
  } else {
    ncat <- rep(1, p)
  }
  maxcat <- max(ncat)
  if (maxcat > 32)
    stop("Can not handle categorical predictors with more than 32 categories.")
  x <- t(x)
  if(!is.null(xtest)) xtest <- t(xtest)
  
  if(classRF) {
    nclass <- length(levels(y))
    if(!all(levels(y) == levels(ytest)))
      stop("y and ytest must have the same levels")
    if (missing(cutoff)) {
      cutoff <- rep(1 / nclass, nclass)
    } else {
      if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 0) ||
         length(cutoff) != nclass) {
        stop("Incorrect cutoff specified.")
      }
      if (!is.null(names(cutoff))) {
        if (!all(names(cutoff) %in% levels(y))) {
          stop("Wrong name(s) for cutoff")
        }
        cutoff <- cutoff[levels(y)]
      }
    }
    if (!is.null(classwt)) {
      if (length(classwt) != nclass)
        stop("length of classwt not equal to number of classes")
      ## If classwt has names, match to class labels.
      if (!is.null(names(classwt))) {
        if (!all(names(cutoff) %in% levels(y))) {
          stop("Wrong name(s) for cutoff")
        }
        classwt <- classwt[levels(y)]
      }
      if (any(classwt <= 0)) stop("classwt must be positive")
      ipi <- 1
    } else {
      classwt <- rep(1, nclass)
      ipi <- 0
    }
  } else addclass <- 0  
  
  if(outscale) {
    outscale <- 1
    proximity <- TRUE
    outlier <- rep(0, n)
  } else {
    outlier <- 0
    outscale <- 0
  }

  if(proximity) {
    prox <- matrix(0, n, n)
    proxts <- if (!is.null(xtest)) matrix(0, ntest, ntest + n) else 0
  } else {
    prox <- proxts <- 0
  }
  
  if (importance) {
    if (classRF) {
      impout <- matrix(0, p, nclass + 2)
    } else {
      impout <- matrix(0, p, 2)
    }
  } else {
      impout <- numeric(p)
  }

  nsample <- if (addclass == 0) n else 2*n
  if (classRF) {
    nboot <- if (missing(sampsize)) nsample else sum(sampsize)
    nrnodes <- 2 * trunc(nboot / nodesize) + 1
  } else {
    ## For regression trees, need to do this to get maximal trees.
    nrnodes <- 2*nsample + 1
  }

  testdat <- !is.null(xtest)
  if(testdat) {
    if(is.null(ytest)) {
      ytest <- labelts <- 0
    } else {
      labelts <- TRUE
    }
  } else {
    ntest <- xtest <- ytest <- 1
    labelts <- FALSE
  }

  if(keep.forest) nt <- ntree else nt <- 1

  if(labelts) error.test <- numeric(ntree) else error.test <- 0
  
  if(classRF) {
    if (!missing(sampsize)) {
      nsum <- sum(sampsize)
      if (any(sampsize <= 0) || nsum == 0) stop("Bad sampsize specification")
      ## If sampsize has names, match to class labels.
      if (!is.null(names(sampsize))) {
        sampsize <- sampsize[levels(y)]
      }
    } else {
      sampsize <- rep(0, nclass)
    }
    rfout <- .C("rf",
                x = as.double(x),
                xdim = as.integer(c(p, n)),
                y = as.integer(y),
                nclass = as.integer(nclass),
                ncat = as.integer(ncat), 
                maxcat = as.integer(maxcat),
                sampsize = as.integer(sampsize),
                Options = as.integer(c(addclass,
                                       importance,
                                       proximity,
                                       oob.prox,
                                       outscale,
                                       do.trace,
                                       keep.forest)),
                ntree=as.integer(ntree),
                mtry=as.integer(mtry),
                ipi=as.integer(ipi),
                classwt=as.double(classwt),
                cutoff = as.double(cutoff),
                nodesize=as.integer(nodesize),
                outlier=as.double(outlier),
                outcl=as.integer(numeric(nsample)),
                counttr=as.integer(numeric(nclass * nsample)),
                prox=as.double(prox),
                impout=as.double(impout),
                nrnodes = as.integer(nrnodes),
                ndbigtree = as.integer(numeric(ntree)),
                nodestatus = as.integer(numeric(nt * nrnodes)),
                bestvar = as.integer(numeric(nt * nrnodes)),
                treemap = as.integer(numeric(nt * 2 * nrnodes)),
                nodepred = as.integer(numeric(nt * nrnodes)),
                xbestsplit = as.double(numeric(nt * nrnodes)),
                pid = as.double(numeric(max(2, nclass))),
                errtr = as.double(numeric(ntree)),
                testdat = as.integer(testdat),
                xts = as.double(xtest),
                clts = as.integer(ytest),
                nts = as.integer(ntest),
                countts = as.double(numeric(nclass * ntest)),
                outclts = as.integer(numeric(ntest)),
                labelts = as.integer(labelts),
                proxts = as.double(proxts),
                errts = as.double(error.test),
                DUP=FALSE,
                PACKAGE="randomForest")
    if(addclass == 0) {
      if(keep.forest) {
        ## deal with the random forest outputs
        max.nodes <- max(rfout$ndbigtree)
        treemap <- array(rfout$treemap, dim = c(2, nrnodes, ntree))
        treemap <- aperm(treemap, c(2,1,3))[1:max.nodes,,]
      }
      ## Turn the predicted class into a factor like y.
      out.class <- levels(y)[rfout$outcl]
      names(out.class) <- x.row.names
      con <- table(observed = y, predicted = out.class)
      con <- cbind(con, class.error = 1 - diag(con)/rowSums(con))
    }
    out.votes <- t(matrix(rfout$counttr, nclass, nsample))[1:n, ]
    if(norm.votes) 
      out.votes <- t(apply(out.votes, 1, function(x) x/sum(x)))
    dimnames(out.votes) <- list(x.row.names, levels(y))
    if(testdat) {
      out.class.ts <- levels(y)[rfout$outclts]
      names(out.class.ts) <- xts.row.names
      out.votes.ts <- t(matrix(rfout$countts, nclass, ntest))
      dimnames(out.votes.ts) <- list(xts.row.names, levels(y))
      if (norm.votes)
        out.votes.ts <- t(apply(out.votes.ts, 1, function(x) x/sum(x)))
      if(labelts) {
        testcon <- table(observed = ytest, predicted = out.class.ts)
        testcon <- cbind(testcon,
                         class.error = 1 - diag(testcon)/rowSums(testcon))
      }
    }
    out <- list(call = match.call(),
                type = ifelse(addclass == 0, "classification",
                  "unsupervised"),
                predicted = if(addclass == 0) out.class else NULL,
                err.rate = if(addclass == 0) rfout$errtr else NULL, 
                confusion = if(addclass == 0) con else NULL,
                votes = out.votes,
                classes = levels(y),
                importance = if(importance) 
                  matrix(rfout$impout, p, nclass+2,
                         dimnames = list(x.col.names,
                           c(levels(y), "MeanDecreaseAccuracy",
                             "MeanDecreaseGini")))
                else structure(rfout$impout, names=x.col.names),
                proximity = if(proximity) matrix(rfout$prox, n, n,
                  dimnames = list(x.row.names, x.row.names)) else NULL,
                outlier = if(outscale) rfout$outlier else NULL,
                ntree = ntree,
                mtry = mtry,
                forest = if(addclass > 0 || !keep.forest) NULL else {
                  list(ndbigtree = rfout$ndbigtree, 
                       nodestatus = matrix(rfout$nodestatus,
                         nc = ntree)[1:max.nodes,],
                       bestvar = matrix(rfout$bestvar, nc = ntree)[1:max.nodes,],
                       treemap = treemap,
                       nodepred = matrix(rfout$nodepred,
                         nc = ntree)[1:max.nodes,],
                       xbestsplit = matrix(rfout$xbestsplit,
                         nc = ntree)[1:max.nodes,],
                       pid = rfout$pid, cutoff = cutoff, ncat = ncat, maxcat = maxcat, 
                       nrnodes = max.nodes, ntree = ntree, nclass = nclass)
                },
                test = if(!testdat) NULL else list(
                  predicted = out.class.ts,
                  err.rate=if(labelts) rfout$errts else NULL,
                  confusion=if(labelts) testcon else NULL,
                  votes=out.votes.ts,
                  proximity = if(proximity) matrix(rfout$proxts, nrow=ntest,
                    dimnames = list(xts.row.names, c(xts.row.names,
                      x.row.names))) else NULL))
  } else {
    ypred <- numeric(n)
    rfout <- .C("regrf",
                as.double(x),
                as.double(y),
                as.integer(n),
                as.integer(p),
                as.integer(nodesize),
                as.integer(nrnodes),
                as.integer(ntree),
                as.integer(mtry),
                as.integer(importance),
                as.integer(ncat),
                as.integer(do.trace),
                as.integer(proximity),
                as.integer(oob.prox),
                as.integer(corr.bias),
                ypred = as.double(ypred),
                impout = as.double(impout),
                prox = as.double(prox),
                ndbigtree = as.integer(numeric(ntree)),
                nodestatus = as.integer(numeric(nt * nrnodes)),
                treemap = as.integer(numeric(nt * 2 * nrnodes)),
                nodepred = as.double(numeric(nt * nrnodes)),
                bestvar = as.integer(numeric(nt * nrnodes)),
                xbestsplit = as.double(numeric(nt * nrnodes)),
                mse = as.double(numeric(ntree)),
                rsq = as.double(numeric(1)),
                keepf = as.integer(keep.forest),
                testdat = as.integer(testdat),
                xts = as.double(xtest),
                ntest = as.integer(ntest),
                yts = as.double(ytest),
                labelts = as.integer(labelts),
                ytestpred = as.double(numeric(ntest)),
                proxts = as.double(proxts),
                msets = as.double(error.test),
                coef = as.double(numeric(2)),
                DUP=FALSE,
                PACKAGE="randomForest")[c(15:25, 32:35)]
    ## Format the forest component, if present.
    if(keep.forest) {
      rfout$nodestatus <- matrix(rfout$nodestatus, ncol = ntree)
      rfout$bestvar <- matrix(rfout$bestvar, ncol = ntree)
      rfout$nodepred <- matrix(rfout$nodepred, ncol = ntree)
      rfout$xbestsplit <- matrix(rfout$xbestsplit, ncol = ntree)
      rfout$treemap <- aperm(array(rfout$treemap, dim = c(2, nrnodes, ntree)),
                             c(2, 1, 3))
    }

    out <- list(call = match.call(),
                type = "regression",
                predicted = structure(rfout$ypred, names=x.row.names),
                mse = rfout$mse,
                rsq = 1 - rfout$mse / (var(y)*(n-1)/n),
                importance = if(importance) structure(matrix(rfout[[2]],p,2),
                  dimnames=list(x.col.names, c("%IncMSE","IncNodePurity"))) else structure(rfout[[2]], names=x.col.names),
                proximity = if(proximity) structure(rfout$prox,
                  dim = c(n, n), dimnames = list(x.row.names, x.row.names)) else NULL,
                ntree = ntree,
                mtry = mtry,
                forest = if(keep.forest) c(rfout[4:9], list(ncat = ncat),
                  list(nrnodes=nrnodes), list(ntree=ntree)) else NULL,
                coefs = if (corr.bias) rfout$coef else NULL,
                test = if(testdat) {
                  list(predicted = structure(rfout$ytestpred,
                         names=xts.row.names),
                  mse = if(labelts) rfout$msets else NULL,
                  rsq = if(labelts) 1 - rfout$msets / (var(ytest)*(n-1)/n) else NULL,
                  proximity = if(proximity) 
                       matrix(rfout$proxts / ntree, nrow = ntest,
                              dimnames = list(xts.row.names, c(xts.row.names,
                                x.row.names))) else NULL)
                  } else NULL)
  }
  class(out) <- "randomForest"
  return(out)
}
