"randomForest.default" <-
function(x, y=NULL, addclass=0, ntree=100,
                                 mtry=ceiling(sqrt(ncol(x))),
                                 classwt=NULL, nodesize=1, importance=FALSE,
                                 proximity=FALSE, outscale=FALSE,
                                 norm.votes=TRUE, do.trace=FALSE, ...)
{
  n <- nrow(x)
  p <- ncol(x)
  if(n == 0)
    stop("data (x) has 0 rows")
  x.row.names <- row.names(x)
  if(mtry > p) {
    mtry <- p
    warning("mtry can not be larger than number of predictors.  Reset to equal to number of predictors")
  }
  if(!is.null(y)) {
    if(length(y) != n) stop("length of response must be the same as predictors")
    if(!is.factor(y)) stop("response must be factor")
    addclass <- 0
  } else {
    if(addclass == 0) addclass <- 1
    y <- factor(c(rep(1, n), rep(2, n)))
    x <- rbind(x, x)
  }
  if(!is.element(addclass,0:2))
    stop("addclass can only take on values 0, 1, or 2")
  
  if(any(is.na(x))) stop("NA not permitted in predictors")
  if(any(is.na(y))) stop("NA not permitted in response")
  nclass <- length(levels(y))
  if(is.data.frame(x)) {
    ncat <- sapply(x, function(x) if(is.factor(x)) length(levels(x)) else 1)
    x <- data.matrix(x)
  } else {
    ncat <- rep(1, p)
  }
  maxcat <- max(ncat)
  x <- t(x)
  
  if(!is.null(classwt)) {
    if(length(classwt) != nclass)
      stop("length of classwt not equal to number of classes")
    if(any(classwt<=0)) stop("classwt must be positive")
    ipi <- 1
  } else {
    classwt <- rep(1, nclass)
    ipi <- 0
  }

  if(importance) impout <- matrix(0, p, 4) else impout <- 0
  
  if(outscale) {
    outscale <- 1
    proximity <- TRUE
    outlier <- rep(0, n)
  } else {
    outlier <- 0
    outscale <- 0
  }
  
  if(proximity) prox <- matrix(0, n, n) else prox <- 0

  nsample <- if(addclass == 0) n else 2*n
  nrnodes <- 2 * (nsample/nodesize) + 1
  outcl <- rep(0, nsample)

  rfout <- .C("rf",
              x=as.double(x),
              p=as.integer(p),
              n=as.integer(n),
              y=as.integer(y),
              nclass=as.integer(nclass),
              ncat=as.integer(ncat), 
              maxcat=as.integer(maxcat),
              addclass=as.integer(addclass),
              ntree=as.integer(ntree),
              mtry=as.integer(mtry),
              ipi=as.integer(ipi),
              classwt=as.double(classwt),
              nodesize=as.integer(nodesize),
              importance=as.integer(importance),
              proximity=as.integer(proximity),
              outscale=as.integer(outscale),
              outlier=as.double(outlier),
              outcl=as.integer(numeric(nsample)),
              counttr=as.integer(numeric(nclass * nsample)),
              prox=as.double(prox),
              impout=as.double(impout), 
              trace=as.integer(as.numeric(do.trace)),
              ndbigtree = as.integer(numeric(ntree)),
              nodestatus = as.integer(numeric(ntree * nrnodes)),
              bestvar = as.integer(numeric(ntree * nrnodes)),
              treemap = as.integer(numeric(ntree * 2 * nrnodes)),
              nodeclass = as.integer(numeric(ntree * nrnodes)),
              xbestsplit = as.double(numeric(ntree * nrnodes)),
              pid = as.double(numeric(max(2, nclass))))

  if(addclass == 0) {
    ## deal with the random forest outputs
    max.nodes <- max(rfout$ndbigtree)
    treemap <- array(rfout$treemap, dim = c(2, nrnodes, ntree))
    treemap <- aperm(treemap, c(2,1,3))[1:max.nodes,,]
    ## Turn the predicted class into a factor like y.
  }
  out.class <- factor(rfout$outcl, levels = 1:length(levels(y)), 
                      labels = levels(y))
  out.votes <- t(matrix(rfout$counttr, nclass, nsample))[1:n, ]
  dimnames(out.votes) <- list(x.row.names, levels(out.class))
  if(norm.votes) out.votes <- t(apply(out.votes, 1, function(x) x/sum(x)))
  
  out <- list(call = match.call(),
              supervised = addclass == 0,
              predicted = if(addclass == 0) out.class else NULL,
              err.rate = mean(y != out.class),
              confusion = if(addclass == 0) table(observed=y, predicted=out.class) else NULL,
              votes = out.votes,
              importance = if(importance) matrix(rfout$impout,p,4) else NULL,
              proximity = if(proximity) matrix(rfout$prox, n, n) else NULL,
              outlier = if(outscale) rfout$outlier else NULL,
              ntree = ntree,
              mtry = mtry,
              forest = if(addclass > 0) NULL else {
                list(ndbigtree = rfout$ndbigtree, 
                     nodestatus = matrix(rfout$nodestatus,
                       nc = ntree)[1:max.nodes,],
                     bestvar = matrix(rfout$bestvar, nc = ntree)[1:max.nodes,],
                     treemap = treemap,
                     nodeclass = matrix(rfout$nodeclass,
                       nc = ntree)[1:max.nodes,],
                     xbestsplit = matrix(rfout$xbestsplit,
                       nc = ntree)[1:max.nodes,],
                     pid = rfout$pid, ncat = ncat, maxcat = maxcat, 
                     nrnodes = max.nodes, ntree = ntree, nclass = nclass)
                })
  class(out) <- "randomForest"
  return(out)
}
