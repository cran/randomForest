\name{randomForest}
\alias{randomForest}
\alias{randomForest.formula}
\alias{randomForest.default}
\alias{print.randomForest}

\title{Breiman's Random Forest classifier}
\description{
  \code{randomForest} implements Breiman's random forest algorithm (based on
  Breiman and Cutler's original Fortran code) for classification.  It
  can also be used in unsupervised mode for locating outliers or
  assessing proximities among data points.
}
\usage{
\method{randomForest}{formula}(formula, data=NULL, subset, ...)
\method{randomForest}{default}(x, y=NULL, addclass=0, ntree=100,
mtry=ceiling(sqrt(ncol(x))), classwt=NULL, nodesize=1, importance=FALSE,
proximity=FALSE, outscale=FALSE, norm.votes=TRUE, do.trace=FALSE, ...)
\method{print}{randomForest}(x, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which
    \code{randomForest} is called from.}
  \item{subset}{an index vector indicating which rows should be used.}
  \item{x}{a data frame or a matrix of predictors (for the
    \code{print} method, an \code{randomForest} object).}
  \item{y}{A response vector, must be a factor, with the levels
    specifying the class membership.  If omitted, \code{randomForest}
    will run in unsupervised mode with \code{addclass=1} (unless
    explicitly set otherwise).}
  \item{addclass}{\code{=0} (default) do not add a synthetic class to
    the data. \code{=1} label the input data as class 1 and add a
    synthetic class by randomly sampling from the product of empirical
    marginal distributions of the input.  \code{=2} is similar to
    \code{=1}, but the synthetic data are sampled from the uniform
    hyperrectangle that contain the input. }
  \item{ntree}{Number of trees to grow.  This should not be set to too
    small a number, to ensure that every input row gets predicted at
    least a few times. }
  \item{mtry}{Number of variables randomly sampled as candidates at each
    split. }
  \item{classwt}{Priors of the classes.  Need not add up to one. }
  \item{nodesize}{Minimum size of terminal nodes.  Setting this number
    larger causes smaller trees to be grown (and thus take less time).}
  \item{importance}{Should importance of predictors be assessed? }
  \item{proximity}{Should proximity measure among the rows be calculated? }
  \item{outscale}{Should outlyingness of rows be assessed? }
  \item{norm.votes}{If \code{TRUE} (default), the final result of votes
    are expressed as fractions.  If \code{FALSE}, raw vote counts are
    returned (useful for combining results from different runs).}
  \item{do.trace}{If set to \code{TRUE}, give a more verbose output as
    \code{randomForest} is run. }
  \item{...}{optional parameters to be passed to the low level function
    \code{randomForest.default}.}
}
%\details{
%}
\value{
  An object of class \code{randomForest}, which is a list with the
  following components:
  \item{call}{the original call to \code{randomForest}}
  \item{supervised}{\code{TRUE} if input data have class labels,
    \code{FALSE} otherwise.}
  \item{predicted}{the predicted class of the input data.}
  \item{err.rate}{final error rate of the prediction on the input data.}
  \item{confusion}{the confusion matrix of the prediction.}
  \item{votes}{a matrix with one row for each input data point and one
    column for each class, giving the fraction or number of `votes' from
    the random forest.}
  \item{importance}{a matrix with four columns, each a different measure
    of importance of the predictors (\code{NULL} if
    \code{importance=FALSE} when \code{randomForest} is called, this
    component is set to \code{NULL}).}
  \item{proximity}{if \code{proximity=TRUE} when \code{randomForest} is
    called, a matrix of proximity measures among the input (based on the
    frequency that pairs of data points are in the same terminal
    nodes).}
  \item{outlier}{if \code{outscale=TRUE} when \code{randomForest} is
    called, a vector indicating how outlying the data points are (based
    on the proximity measures).}
  \item{ntree}{number of trees grown.}
  \item{mtry}{number of predictors sampled for spliting at each node.}
  \item{forest}{a list that contains the entire forest; \code{NULL} if
    \code{randomForest} is run in unsupervised mode.}  
}
\references{
  Breiman, L., \emph{Random Forests}, 
  \url{http://oz.berkeley.edu/users/breiman/randomforest2001.pdf}.
}
\author{Andy Liaw \email{andy\_liaw@merck.com} and Mathew Wiener
  \email{mathew\_wiener@merck.com}, based on original Fortran code by
  Leo Breiman and Adele Cutler.}

\seealso{\code{\link{predict.randomForest}}}

\examples{
data(iris)
iris.rf <- randomForest(Species ~ ., data=iris)
print(iris.rf)
}
\keyword{classif}% at least one, from doc/KEYWORDS
