"randomForest.default" <-
    function(x, y=NULL,  xtest=NULL, ytest=NULL, addclass=0, ntree=500,
             mtry=ifelse(!is.null(y) && !is.factor(y),
             max(floor(ncol(x)/3), 1), 
             floor(sqrt(ncol(x)))), replace=TRUE, classwt=NULL, cutoff,
             sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x)),
             nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1, 
             importance=FALSE, localImp=FALSE,
             proximity=FALSE, oob.prox=proximity,
             outscale=FALSE, norm.votes=TRUE, do.trace=FALSE,
             keep.forest=is.null(xtest), corr.bias=FALSE, ...)
{
    classRF <- is.null(y) || is.factor(y)
    if (!classRF && length(unique(y)) <= 5) {
        warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
    }
    if (classRF && !is.null(y) && addclass==0 && length(unique(y)) < 2) {
        stop("Need at least two classes to do classification.")
    }
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
    
    ## overcome R's lazy evaluation:
    keep.forest <- keep.forest
    
    testdat <- !is.null(xtest)
    if (testdat) {
        if (ncol(x) != ncol(xtest))
            stop("x and xtest must have same number of columns") 
        ntest <- nrow(xtest)
        xts.row.names <- rownames(xtest)
    }
    
    if(mtry > p) {
        mtry <- p
        warning("mtry can not be larger than number of predictors.  Reset to equal to number of predictors")
    }
    if (!is.null(y)) {
        if (length(y) != n) stop("length of response must be the same as predictors")
        addclass <- 0
    } else {
        if (addclass == 0) addclass <- 1
        y <- factor(c(rep(1, n), rep(2, n)))
        x <- rbind(x, x)
        keep.forest <- FALSE
    }
  if (!is.element(addclass,0:2))
      stop("addclass can only take on values 0, 1, or 2")
  
  if (any(is.na(x))) stop("NA not permitted in predictors")
  if (testdat && any(is.na(xtest))) stop("NA not permitted in xtest")
  if (any(is.na(y))) stop("NA not permitted in response")
  if (!is.null(ytest) && any(is.na(ytest))) stop("NA not permitted in ytest")

    if (is.data.frame(x)) {
        ncat <- sapply(x, function(x) if(is.factor(x) && !is.ordered(x))
                       length(levels(x)) else 1)
        x <- data.matrix(x)
        if(testdat) {
            if(!is.data.frame(xtest))
                stop("xtest must be data frame if x is")
            ncatts <- sapply(xtest, function(x) if(is.factor(x) &&
                                                   !is.ordered(x))
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
    
    if (classRF) {
        nclass <- length(levels(y))
        if (!is.null(ytest)) {
            if (!is.factor(ytest)) stop("ytest must be a factor")
            if (!all(levels(y) == levels(ytest)))
                stop("y and ytest must have the same levels")
        }
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
        prox <- matrix(0.0, n, n)
        proxts <- if (testdat) matrix(ntest, ntest + n) else double(1)
    } else {
        prox <- proxts <- double(1)
    }

    if (localImp) {
        importance <- TRUE
        impmat <- matrix(0, p, n)
    } else impmat <- double(1)
    
    if (importance) {
        if (classRF) {
            impout <- matrix(0.0, p, nclass + 2)
            impSD <- matrix(0.0, p, nclass + 1)
        } else {
            impout <- matrix(0.0, p, 2)
            impSD <- double(p)
            names(impSD) <- x.col.names
        }
    } else {
        impout <- double(p)
        impSD <- double(1)
    }
    
    nsample <- if (addclass == 0) n else 2*n
    Stratify <- length(sampsize) > 1
    if ((!Stratify) && sampsize > nrow(x)) stop("sampsize too large")
    if (Stratify && (!classRF)) stop("sampsize should be of length one")
    if (classRF) {
        if (Stratify) {
            nsum <- sum(sampsize)
            if (length(sampsize) > nlevels(y))
                stop("sampsize has too many elements.")
            if (any(sampsize <= 0) || nsum == 0)
                stop("Bad sampsize specification")
            ## If sampsize has names, match to class labels.
            if (!is.null(names(sampsize))) {
                sampsize <- sampsize[levels(y)]
            }
            if (any(sampsize > table(y)))
              stop("sampsize can not be larger than class frequency")
        } else {
            nsum <- sampsize
        }
        nrnodes <- 2 * trunc(nsum / nodesize) + 1
    } else {
        ## For regression trees, need to do this to get maximal trees.
        nrnodes <- 2 * trunc(sampsize/max(1, nodesize - 3)) + 1
    }
    

    x <- t(x)
    storage.mode(x) <- "double"
    if (testdat) {
        xtest <- t(xtest)
        storage.mode(xtest) <- "double"
        if (is.null(ytest)) {
            ytest <- labelts <- 0
        } else {
            labelts <- TRUE
        }
    } else {
        xtest <- double(1)
        ytest <- double(1)
        ntest <- 1
        labelts <- FALSE
    }
    nt <- if (keep.forest) ntree else 1

    if (classRF) {
        error.test <- if (labelts) double((nclass+1) * ntree) else double(1)
        rfout <- .C("classRF",
                    x = x,
                    xdim = as.integer(c(p, n)),
                    y = as.integer(y),
                    nclass = as.integer(nclass),
                    ncat = as.integer(ncat), 
                    maxcat = as.integer(maxcat),
                    sampsize = as.integer(sampsize),
                    Options = as.integer(c(addclass,
                    importance,
                    localImp,
                    proximity,
                    oob.prox,
                    outscale,
                    do.trace,
                    keep.forest,
                    replace,
                    Stratify)),
                    ntree = as.integer(ntree),
                    mtry = as.integer(mtry),
                    ipi = as.integer(ipi),
                    classwt = as.double(classwt),
                    cutoff = as.double(cutoff),
                    nodesize = as.integer(nodesize),
                    outlier = as.double(outlier),
                    outcl = integer(nsample),
                    counttr = integer(nclass * nsample),
                    prox = prox,
                    impout = impout,
                    impSD = impSD,
                    impmat = impmat,
                    nrnodes = as.integer(nrnodes),
                    ndbigtree = integer(ntree),
                    nodestatus = integer(nt * nrnodes),
                    bestvar = integer(nt * nrnodes),
                    treemap = integer(nt * 2 * nrnodes),
                    nodepred = integer(nt * nrnodes),
                    xbestsplit = double(nt * nrnodes),
                    pid = double(max(2, nclass)),
                    errtr = double((nclass+1) * ntree),
                    testdat = as.integer(testdat),
                    xts = as.double(xtest),
                    clts = as.integer(ytest),
                    nts = as.integer(ntest),
                    countts = double(nclass * ntest),
                    outclts = as.integer(numeric(ntest)),
                    labelts = as.integer(labelts),
                    proxts = proxts,
                    errts = error.test,
                    DUP=FALSE,
                    PACKAGE="randomForest")[-1]
        if (addclass == 0) {
            if (keep.forest) {
        ## deal with the random forest outputs
                max.nodes <- max(rfout$ndbigtree)
                treemap <- array(rfout$treemap, dim = c(2, nrnodes, ntree))
                treemap <- aperm(treemap, c(2,1,3))[1:max.nodes,,]
            }
            ## Turn the predicted class into a factor like y.
            out.class <- factor(rfout$outcl, levels=1:nclass,
                                label=levels(y))
            names(out.class) <- x.row.names
            con <- table(observed = y,
                         predicted = out.class)[levels(y), levels(y)]
            con <- cbind(con, class.error = 1 - diag(con)/rowSums(con))
        }
        out.votes <- t(matrix(rfout$counttr, nclass, nsample))[1:n, ]
        oob.times <- rowSums(out.votes)
        if(norm.votes) 
            out.votes <- t(apply(out.votes, 1, function(x) x/sum(x)))
        dimnames(out.votes) <- list(x.row.names, levels(y))
        if(testdat) {
            out.class.ts <- factor(rfout$outclts, levels=1:nclass,
                                   label=levels(y))
            names(out.class.ts) <- xts.row.names
            out.votes.ts <- t(matrix(rfout$countts, nclass, ntest))
            dimnames(out.votes.ts) <- list(xts.row.names, levels(y))
            if (norm.votes)
                out.votes.ts <- t(apply(out.votes.ts, 1,
                                        function(x) x/sum(x)))
            if (labelts) {
                testcon <- table(observed = ytest,
                                 predicted = out.class.ts)[levels(y), levels(y)]
                testcon <- cbind(testcon,
                                 class.error = 1 - diag(testcon)/rowSums(testcon))
            }
        }
        out <- list(call = match.call(),
                    type = ifelse(addclass == 0, "classification",
                    "unsupervised"),
                    predicted = if (addclass == 0) out.class else NULL,
                    err.rate = if (addclass == 0) t(matrix(rfout$errtr,
                                   nclass+1,
                                   ntree, dimnames=list(c("OOB", levels(y)),
                                          NULL))) else NULL, 
                    confusion = if(addclass == 0) con else NULL,
                    votes = out.votes,
                    oob.times = oob.times,
                    classes = levels(y),
                    importance = if (importance) 
                    matrix(rfout$impout, p, nclass+2,
                           dimnames = list(x.col.names,
                           c(levels(y), "MeanDecreaseAccuracy",
                             "MeanDecreaseGini")))
                    else structure(rfout$impout, names=x.col.names),
                    importanceSD = if (importance)
                    matrix(rfout$impSD, p, nclass + 1,
                           dimnames = list(x.col.names,
                           c(levels(y), "MeanDecreaseAccuracy")))
                    else NULL,
                    localImportance = if (localImp)
                    matrix(rfout$impmat, p, n,
                           dimnames = list(x.col.names,x.row.names)) else NULL,
                    proximity = if (proximity) matrix(rfout$prox, n, n,
                    dimnames = list(x.row.names, x.row.names)) else NULL,
                    outlier = if (outscale) rfout$outlier else NULL,
                    ntree = ntree,
                    mtry = mtry,
                    forest = if (addclass > 0 || !keep.forest) NULL else {
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
                    err.rate = if (labelts) t(matrix(rfout$errts, nclass+1,
                    ntree,
                dimnames=list(c("Test", levels(y)), NULL))) else NULL,
                    confusion = if (labelts) testcon else NULL,
                    votes = out.votes.ts,
                    proximity = if(proximity) matrix(rfout$proxts, nrow=ntest,
                    dimnames = list(xts.row.names, c(xts.row.names,
                    x.row.names))) else NULL))
    } else {
        rfout <- .C("regRF",
                    x,
                    as.double(y),
                    as.integer(n),
                    as.integer(p),
                    as.integer(sampsize),
                    as.integer(nodesize),
                    as.integer(nrnodes),
                    as.integer(ntree),
                    as.integer(mtry),
                    as.integer(c(importance, localImp)),
                    as.integer(ncat),
                    as.integer(do.trace),
                    as.integer(proximity),
                    as.integer(oob.prox),
                    as.integer(corr.bias),
                    ypred = double(n),
                    impout = impout,
                    impmat = impmat,
                    impSD = impSD,
                    prox = prox,
                    ndbigtree = integer(ntree),
                    nodestatus = integer(nt * nrnodes),
                    treemap = integer(nt * 2 * nrnodes),
                    nodepred = double(nt * nrnodes),
                    bestvar = integer(nt * nrnodes),
                    xbestsplit = double(nt * nrnodes),
                    mse = double(ntree),
                    keepf = as.integer(keep.forest),
                    replace = as.integer(replace),
                    testdat = as.integer(testdat),
                    xts = xtest,
                    ntest = as.integer(ntest),
                    yts = as.double(ytest),
                    labelts = as.integer(labelts),
                    ytestpred = double(ntest),
                    proxts = proxts,
                    msets = double(if (labelts) ntree else 1),
                    coef = double(2),
                    oob.times = integer(n),
                    DUP=FALSE,
                    PACKAGE="randomForest")[c(16:27, 35:39)]
        ## Format the forest component, if present.
        if (keep.forest) {
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
                    rsq = 1 - rfout$mse / (var(y) * (n-1) / n),
                    oob.times = rfout$oob.times,
                    importance = if (importance) matrix(rfout$impout, p, 2,
                    dimnames=list(x.col.names,
                                  c("%IncMSE","IncNodePurity"))) else
                        structure(rfout$impout, names=x.col.names),
                    importanceSD=if (importance) rfout$impSD else NULL,
                    localImportance = if (importance)
                    matrix(rfout$impmat, p, n, dimnames=list(x.col.names,
                                               x.row.names)) else NULL,
                    proximity = if (proximity) matrix(rfout$prox, n, n,
                    dimnames = list(x.row.names, x.row.names)) else NULL,
                    ntree = ntree,
                    mtry = mtry,
                    forest = if (keep.forest)
                    c(rfout[c("ndbigtree", "nodestatus", "treemap",
                              "nodepred", "bestvar", "xbestsplit")],
                      list(ncat = ncat), list(nrnodes=nrnodes),
                      list(ntree=ntree)) else NULL,
                    coefs = if (corr.bias) rfout$coef else NULL,
                    test = if(testdat) {
                        list(predicted = structure(rfout$ytestpred,
                             names=xts.row.names),
                             mse = if(labelts) rfout$msets else NULL,
                             rsq = if(labelts) 1 - rfout$msets /
                                        (var(ytest)*(n-1)/n) else NULL,
                             proximity = if (proximity) 
                             matrix(rfout$proxts / ntree, nrow = ntest,
                                    dimnames = list(xts.row.names,
                                    c(xts.row.names,
                                    x.row.names))) else NULL)
                    } else NULL)
    }
    class(out) <- "randomForest"
    return(out)
}
