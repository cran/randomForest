"predict.randomForest" <-
function (object, newdata, type = "class", norm.votes = TRUE, ...) 
{
  if (!inherits(object, "randomForest")) 
    stop("object not of class randomForest")
  out.type <- charmatch(tolower(type), c("class", "prob", "vote"))
  if(is.na(out.type))
	stop("type must be one of 'class', 'prob', 'vote'")
  if(out.type == 2)  norm.votes <- TRUE
  if (missing(newdata)) {
	if(out.type == 1) return(object$predicted)
	if(norm.votes) return(t(apply(object$votes, 1, function(x){x/sum(x)})))
	  else return(object$votes)
  }
  if (!object$supervised) 
    stop("Can't predict unsupervised forest.")
  if (inherits(object, "randomForest.formula")) {
    newdata <- as.data.frame(newdata)
    rn <- row.names(newdata)
    Terms <- delete.response(object$terms)
    m <- model.frame(Terms, newdata, na.action = na.omit)
    keep <- match(row.names(m), rn)
    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0) 
      x <- x[, -xint, drop = FALSE]
  } else {
    if (is.null(dim(newdata))) 
      dim(newdata) <- c(1, length(newdata))
    x <- as.matrix(newdata)
    if(nrow(x) == 0) stop("newdata has 0 rows")
    if (any(is.na(x))) 
      stop("missing values in newdata")
    keep <- 1:nrow(x)
    rn <- rownames(x)
  }
  if (is.data.frame(x)) {
    cat.new <- sapply(x, function(x) if (is.factor(x)) 
                      length(levels(x))
    else 1)
    if (!all.equal(object$ncat, cat.new)) 
      stop("Type of predictors in new data do not match that of the training data.")
  }
  mdim <- ncol(x)
  ntest <- nrow(x)
  ntree <- object$forest$ntree
  maxcat <- object$forest$maxcat
  nclass <- object$forest$nclass
  nrnodes <- object$forest$nrnodes
  show.error <- 0
  x <- t(x)
  t1 <- .Fortran("runforest", mdim = as.integer(mdim), ntest = as.integer(ntest), 
                 nclass = as.integer(object$forest$nclass), maxcat = as.integer(maxcat), 
                 nrnodes = as.integer(nrnodes), labelsts = as.integer(as.numeric(show.error)), 
                 jbt = as.integer(ntree), clts = as.integer(numeric(ntest)), 
                 xts = as.double(x), xbestsplit = as.double(object$forest$xbestsplit), 
                 pid = as.double(object$forest$pid), countts = as.double(numeric(nclass * 
                                                       ntest)), treemap = as.integer(aperm(object$forest$treemap, 
                                                                  c(2, 1, 3))), nodestatus = as.integer(object$forest$nodestatus), 
                 cat = as.integer(object$forest$ncat), cbestsplit = as.integer(numeric(maxcat * 
                                                         nrnodes)), nodeclass = as.integer(object$forest$nodeclass), 
                 jts = as.integer(numeric(ntest)), jet = as.integer(numeric(ntest)), 
                 bestvar = as.integer(object$forest$bestvar), nodexts = as.integer(numeric(ntest)), 
                 ndbigtree = as.integer(object$forest$ndbigtree))
  out.class.votes <- t(matrix(t1$countts, nr = nclass, nc = ntest))
  if (norm.votes) 
    out.class.votes <- t(apply(out.class.votes, 1, function(x) x/sum(x)))
    z <- matrix(NA, length(rn), nclass, dimnames = list(rn, 
                                          levels(object$predicted)))
    z[keep, ] <- out.class.votes
   res <- z
  if(out.type == 1){
    out.class <- max.col(out.class.votes)
    out.class <- if (ncol(z) > 1) 
      levels(object$predicted)[max.col(z)]
    else levels(object$predicted)[1 + (z > 0.5)]
    out.class <- factor(out.class, levels = levels(object$predicted))
   res <- out.class
  }
  res
}
