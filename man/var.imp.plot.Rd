\name{var.imp.plot}
\alias{var.imp.plot}
\title{Variable Importance Plot}
\description{
  Dotchart of variable importance as measured by a Random Forest
}
\usage{
var.imp.plot(x, sort=TRUE, 
    n.var=min(50, if(is.null(dim(x$importance)))
            length(x$importance) else nrow(x$importance)),
    class=NULL,
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
mtcars.rf <- randomForest(mpg ~ ., data=mtcars, ntree=1000, keep=FALSE,
                          importance=TRUE)
var.imp.plot(mtcars.rf)
}
\author{Andy Liaw \email{andy\_liaw@merck.com}}
\keyword{regression}
\keyword{classif}
\keyword{tree}
