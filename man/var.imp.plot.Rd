\name{var.imp.plot}
\alias{var.imp.plot}
\title{Variable Importance Plot}
\description{
  Barchart of variable importance as measured by a Random Forest
}
\usage{
var.imp.plot(x, sort=TRUE, label=TRUE, main=deparse(substitute(x)), ...)
}
\arguments{
  \item{x}{An object of class \code{randomForest}.}
  \item{sort}{Should the variables be sorted in decreasing order of
    importance?}
  \item{label}{Should the variable names be added to the importance?}
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
data(mtcars)
var.imp.plot(randomForest(mpg ~ ., data=mtcars, ntree=1000, keep=FALSE,
                          importance=TRUE))
}
\author{Andy Liaw}
\keyword{regression}
\keyword{classif}
\keyword{tree}
