\name{partial.plot}
\alias{partial.plot}
\alias{partial.plot.default}
\alias{partial.plot.randomForest}

\title{Partial dependence plot}

\description{
Partial dependence plot gives a graphical depiction of the marginal
effect of a variable on the class probability (classification) or
response (regression).
}

\usage{
\method{partial.plot}{randomForest}(x, pred.data, x.var, which.class,
    add=FALSE, n.pt=min(length(unique(pred.data[,deparse(substitute(x.var)])),
       51), rug=TRUE, ...)
}

\arguments{
  \item{x}{an object of class \code{randomForest}, which contains a
    \code{forest} component.}
  \item{pred.data}{a data frame used for contructing the plot, usually
    the training data used to contruct the random forest.}
  \item{x.var}{(unquoted) name of the variable for which partial
    dependence is to be examined.}
  \item{which.class}{For classification data, the class to focus on
    (default the first class).}
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

\details{
  The function being plotted is defined as:
  \deqn{
    \tilde{f}(x) = \frac{1}{n} \sum_{i=1}^n f(x, x_{iC}),
  }
  where \eqn{x} is the variable for which partial dependence is sought,
  and \eqn{x_{iC}} is the other variables in the data.  The summand is
  the predicted regression function for regression, and predicted
  probabilities (i.e., fraction of votes) for \code{which.class} for
  classification.  For classification, Friedman (2001) proposed the
  logit, but as the estimated probabilities could be 0, for which the
  logit is undefined, we chose to use just the probabilities for the
  class of interest instead.
}
\note{
  The \code{randomForest} object must contain the \code{forest}
  component; i.e., created with \code{randomForest(...,
    keep.forest=TRUE)}.

  This function runs quite slow for large data sets.
}
\references{
Friedman, J. (2001). Greedy function approximation: the gradient
boosting machine, \emph{Ann. of Stat.}}

\seealso{\code{\link{randomForest}}}

\author{Andy Liaw \email{andy\_liaw@merck.com}}

\examples{
data(airquality)
airquality <- na.omit(airquality)
set.seed(131)
ozone.rf <- randomForest(Ozone ~ ., airquality)
partial.plot(ozone.rf, airquality, Temp)

data(iris)
set.seed(543)
iris.rf <- randomForest(Species~., iris)
partial.plot(iris.rf, iris, Petal.Width, "versicolor")
}
\keyword{classif}
\keyword{regression}
\keyword{tree}

