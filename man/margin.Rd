\name{margin}
\alias{margin}
\alias{plot.margin}
\title{Margins of randomForest Classifier}
\description{
  Compute or plot the margin of predictions from a randomForest classifier.
}
\usage{
margin(rf, observed)
\method{plot}{margin}(x, sort=TRUE, ...)
}
\arguments{
  \item{rf}{an object of class \code{\link{randomForest}}, whose
    \code{type} is not \code{regression}.}
  \item{observed}{the response variable (a factor) used to train
    \code{rf}.}
  \item{x}{an object of class \code{margin}, as returned by the function
    \code{margin}.}
  \item{sort}{Should the data be sorted by their class labels?}
  \item{...}{other graphical parameters to be passed to \code{plot.default}.}
}
\value{
  For \code{margin}, the \emph{margin} of observations from the
  \code{\link{randomForest}} classifier.  The margin of a data point is defined as the
  proportion of votes for the correct class minus maximum proportion of
  votes for the other classes.  Thus under majority votes, positive
  margin means correct classification, and vice versa.
}
\seealso{
  \code{\link{randomForest}}
}
\examples{
set.seed(1)
data(iris)
iris.rf <- randomForest(Species ~ ., iris, keep.forest=FALSE)
plot(margin(iris.rf, iris$Species))
}
\author{Robert Gentlemen, with slight modifications by Andy Liaw}
\keyword{classif}
