\name{var.imp.plot}
\alias{var.imp.plot}
\title{Variable Importance Plot}
\description{
  Barchart of variable importance as measured by a Random Forest
}
\usage{
var.imp.plot(x, sort=TRUE, label=TRUE,
    n.var=min(50, if(is.null(dim(x$importance)))
            length(x$importance) else nrow(x$importance)),
    main=deparse(substitute(x)), ...)
}
\arguments{
  \item{x}{An object of class \code{randomForest}.}
  \item{sort}{Should the variables be sorted in decreasing order of
    importance?}
  \item{label}{Should the variable names be added to the importance?}
  \item{n.var}{How many variables to show? (Ignored if \code{sort=FALSE}.)}
  \item{main}{plot title.}
  \item{...}{Other graphical parameters.}
}
\value{
  Invisibly, the importance of the variables.
}
\note{
If \code{x$type} is \code{classification}, all four measures are plotted
in a 2 by 2 array.
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
