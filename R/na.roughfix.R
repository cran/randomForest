na.roughfix <- function(object, ...)
  UseMethod("na.roughfix")

na.roughfix.data.frame <- function(object, ...) {
##  n <- length(object)
##  vars <- seq(length = n)
  isfac <- sapply(object, is.factor)
  isnum <- sapply(object, is.numeric)
  hasNA <- which(sapply(object, function(x) any(is.na(x))))
  if (any(!(isfac | isnum)))
      stop("na.roughfix only works for numeric or factor")
  for (j in hasNA) {
    if (isfac[j]) {
      freq <- table(object[[j]])
      xmode <- names(freq)[max.col(rbind(freq))]
      object[[j]][is.na(object[[j]])] <- xmode
    } else {
      xmed <- median(object[[j]], na.rm=TRUE)
      object[[j]][is.na(object[[j]])] <- xmed
    }
  }
  object
}

na.roughfix.default <- function(object, ...) {
  if (!is.atomic(object)) 
    return(object)
  d <- dim(object)
  if (length(d) > 2) 
    stop("can't handle objects with more than two dimensions")
  if (all(!is.na(object)))
    return(object)
  if (!is.numeric(object))
    stop("roughfix can only deal with numeric data.")
  if (length(d) == 2) {
      hasNA <- which(apply(object, 2, function(x) any(is.na(x))))
      for (j in hasNA) 
          object[is.na(object[, j]), j] <- median(object[, j], na.rm=TRUE)
  } else {
      object[is.na(object)] <- median(object, na.rm=TRUE)
  }
  object
}
