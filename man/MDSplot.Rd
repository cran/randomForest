\name{MDSplot}
\alias{MDSplot}
\title{Multi-dimensional Scaling Plot of Proximity matrix from randomForest}
\description{
  Plot the scaling coordinates of the proximity matrix from randomForest.
}
\usage{
MDSplot(rf, fac, k=2, ...)
}
\arguments{
  \item{rf}{an object of class \code{\link{randomForest}} that contains
    the \code{proximity} component.}
  \item{fac}{a factor that was used as response to train \code{rf}.}
  \item{k}{number of dimensions for the scaling coordinates.}
  \item{...}{other graphical parameters.}
}
\value{
  The output of \code{\link[mva]{cmdscale}} on 1 - \code{rf$proximity} is
  returned invisibly.
}
\Note{
  If \code{k > 2}, \code{\link{pairs}} is used to produce the
  scatterplot matrix of the coordinates.
}
\seealso{
  \code{\link{randomForest}}
}
\examples{
set.seed(1)
data(iris)
iris.rf <- randomForest(Species ~ ., iris, proximity=TRUE,
                        keep.forest=FALSE)
MDSplot(iris.rf, iris$Species)
}
\author{Robert Gentleman, with slight modifications by Andy Liaw}
\keyword{classif}
\keyword{tree}
