\name{varImpPlot}
\alias{varImpPlot}
\title{Variable Importance Plot}
\description{
  Dotchart of variable importance as measured by a Random Forest
}
\usage{
varImpPlot(x, sort=TRUE, n.var=min(30, if(is.null(dim(x$importance)))
           length(x$importance) else nrow(x$importance)),
           class = NULL, scale=TRUE, xlab="Importance", ylab="",
           main=deparse(substitute(x)), ...) 
}
\arguments{
  \item{x}{An object of class \code{randomForest}.}
  \item{sort}{Should the variables be sorted in decreasing order of
    importance?}
  \item{n.var}{How many variables to show? (Ignored if
    \code{sort=FALSE}.)}
  \item{class}{For classification data, an integer or string indicating
    the class for which variable importance is seeked.}
  \item{scale}{For permutation-based measures, should the measures be
    divided by their ``standard errors''?}
  \item{xlab}{label for the x-axis.}
  \item{ylab}{label for the y-axis.}
  \item{main}{plot title.}
  \item{...}{Other graphical parameters.}
}
\value{
  Invisibly, the importance of the variables.
}
\seealso{
  \code{\link{randomForest}}
}
\examples{
set.seed(4543)
data(mtcars)
mtcars.rf <- randomForest(mpg ~ ., data=mtcars, ntree=1000, keep.forest=FALSE,
                          importance=TRUE)
varImpPlot(mtcars.rf)
}
\author{Andy Liaw \email{andy\_liaw@merck.com}}
\keyword{regression}
\keyword{classif}
\keyword{tree}
