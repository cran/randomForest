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
\method{importance}{randomForest}(x, type=NULL, class=NULL, scale=TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{\link{randomForest}}}.
  \item{type}{either 1 or 2, specifying the type of importance measure
    (1=mean decrease in accuracy, 2=mean decrease in node impurity).}
  \item{class}{for classification problem, which class-specific measure
    to return.}
  \item{scale}{For permutation based measures, should the measures be
    divided their ``standard errors''?}
  \item{...}{not used.}
}
\value{
  If \code{class} and \code{type} are \code{NULL}, A matrix of
  \code{p} rows and \code{nclass + 2} columns (where \code{p} is the
  number of variables in the data and \code{nclass} is the number of
  classes) for classification problem, or \code{p} rows and \code{2}
  columns for regression.  In the classification case, the first
  \code{nclass} columns are the class-specific importance measures
  (based on permutation of out-of-bag data).  The \code{nclass+1}st
  column is the overall importance, and the last column is the overall
  measure based on the decrease in node purity (or `Gini-based'
  measure).

  If either \code{type} or \code{class} is given, it returns a
  (named) vector of importance measure, one for each predictor variable.
}
\details{
  See the documentation for \code{\link{randomForest}} for explanation
  of how the importance measures are computed.
}
\seealso{
  \code{\link{randomForest}}, \code{\link{varImpPlot}}
}
\example{
set.seed(4543)
data(mtcars)
mtcars.rf <- randomForest(mpg ~ ., data=mtcars, ntree=1000, keep=FALSE,
                          importance=TRUE)
importance(mtcars.rf)
importance(mtcars.rf, type=1)
}
\keyword{regression}
\keyword{classif}
\keyword{tree}

