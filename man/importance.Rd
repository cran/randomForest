\name{importance}
\alias{importance}
\alias{importance.default}
\alias{importance.randomForest}
\title{Extract variable importance measure}
\description{
  This is the extractor function for variable importance measures as
  produced by \code{\link{randomForest}}.
}
\usage{
\method{importance}{randomForest}(x, type=2, class, ...)
}
\arguments{
  \item{x}{an object of class \code{\link{randomForest}}}.
  \item{type}{either 1 or 2, specifying the type of importance measure
    (1=mean decrease in accuracy, 2=mean decrease in node impurity).}
  \item{class}{for classification problem, which class-specific measure
    to return.}
  \item{...}{not used.}
}
\value{
  A (named) vector of importance measure, one for each predictor variable.
}
\details{
  See the documentation for \code{\link{randomForest}} for explanation
  of how the importance measures are computed.
}
%\references{
%}
\seealso{
  \code{\link{randomForest}}, \code{\link{var.imp.plot}}
}
set.seed(4543)
data(mtcars)
mtcars.rf <- randomForest(mpg ~ ., data=mtcars, ntree=1000, keep=FALSE,
                          importance=TRUE)
importance(mtcars.rf)
importance(mtcars.rf, type=1)			  
%\author{}
\keyword{regression}
\keyword{classif}
\keyword{tree}

