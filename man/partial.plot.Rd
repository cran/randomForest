\name{partial.plot}
\alias{partial.plot}
\alias{partial.plot.default}
\alias{partial.plot.randomForest}

\title{Partial dependence plot}

\description{
Partial dependence plot gives a graphical depiction of the marginal
effect of a variable on the response, currently only implemented for
random forests.
}

\usage{
\method{partial.plot}{randomForest}(x, pred.data, x.var, add=FALSE, n.pt=min(nrow(pred.data), 51),
             rug=TRUE, ...)
}

\arguments{
  \item{x}{an object of class \code{randomForest}, which contains a
    \code{forest} component.}
  \item{pred.data}{a data frame used for contructing the plot, usually
    the training data used to contruct the random forest.}
  \item{x.var}{(unquoted) name of the variable for which partial
    dependence is to be examined.}
  \item{add}{whether to add to existing plot (\code{TRUE}) or create a
    new plot (\code{FALSE}).}
  \item{n.pt}{if \code{x.var} is continuous, the number of points on the
    grid for evaluating partial dependence.}
  \item{rug}{whether to draw hash marks at the bottom of the plot
    indicating the deciles of \code{x.var}.}
  \item{...}{other graphical parameters to be passed on to \code{plot}
    or \code{lines}.}
}

\value{
A list with two components: \code{x} and \code{y}, which are the values
used in the plot.
}

\note{
  Currently the function is only implemented for regression.
  The \code{randomForest} object must contain the \code{forest}
  component; i.e., created with \code{randomForest(...,
    keep.forest=TRUE)}.

  This function is likely to be slow for large data sets.
}
\references{
Freidman, J. (2001). Greedy function approximation: the gradient
boosting machine, \emph{Ann. of Stat.}}

\seealso{\code{\link{randomForest}}}

\author{Andy Liaw \email{andy\_liaw@merck.com}}

\examples{
data(airquality)
airquality <- na.omit(airquality)
ozone.rf <- randomForest(Ozone ~ ., airquality)
partial.plot(ozone.rf, airquality, Temp)
}

\keyword{regression}
