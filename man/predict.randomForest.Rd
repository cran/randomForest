\name{predict.randomForest}
\alias{predict.randomForest}
\title{predict method for random forest objects}
\description{
  Prediction of test data using random forest.
}
\usage{
\method{predict}{randomForest}(object, newdata, type="response",
  norm.votes=TRUE, proximity=FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{randomForest}, as that
    created by the function \code{randomForest}.}
  \item{newdata}{a data frame or matrix containing new data.}
  \item{type}{one of \code{response}, \code{prob}. or \code{votes},
  indicating the type of output: predicted values, matrix of class
  probabilities, or matrix of vote counts.  \code{class} is allowed, but
  automatically converted to "response", for backward compatibility.}
  \item{norm.votes}{Should the vote counts be normalized (i.e.,
    expressed as fractions)?  Ignored if \code{object$type} is
    \code{regression}.}
  \item{proximity}{Should proximity measures be computed?  An error is
    issued if \code{object$type} is \code{regression}.}
  \item{...}{not used currently.}
}

\value{
  If \code{object$type} is \code{regression}, a vector of predicted
  values is returned.

  If \code{object$type} is \code{classification}, the vector returned
  depends on the argument \code{type}:
  \item{response}{predicted classes (the classes with majority vote).}
  \item{prob}{matrix of class probabilities (one column for each class
  and one row for each input).}
  \item{votes}{matrix of vote counts (one column for each class
  and one row for each new input); either in raw counts or in fractions
  (if \code{norm.votes=TRUE}).}

If \code{proximity=TRUE}, the returned object is a list with two
components: \code{pred} is the prediction (as described above) and
\code{proximity} is the proximitry matrix.  An error is issued if
\code{object$type} is \code{regression}.
}
\references{
  Breiman, L., \emph{Random Forests},
  \url{http://oz.berkeley.edu/users/breiman/randomforest2001.pdf}.
}
\author{ Andy Liaw \email{andy\_liaw@merck.com} and Matthew Wiener
  \email{matthew\_wiener@merck.com}, based on original Fortran code by
  Leo Breiman and Adele Cutler.}

\seealso{\code{\link{randomForest}}}

\examples{
data(iris)
set.seed(111)
n.iris <- nrow(iris)
iris.test <- sample(n.iris, round(n.iris/3), replace=FALSE)
iris.rf <- randomForest(Species ~ ., data=iris[-iris.test,])
iris.pred <- predict(iris.rf, iris[iris.test,])
table(observed = iris[iris.test,"Species"], predicted = iris.pred)
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{regression}
