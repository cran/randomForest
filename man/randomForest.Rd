\name{randomForest}
\alias{randomForest}
\alias{randomForest.formula}
\alias{randomForest.default}
\alias{print.randomForest}

\title{Classification and Regression with Random Forest}
\description{
  \code{randomForest} implements Breiman's random forest algorithm (based on
  Breiman and Cutler's original Fortran code) for classification and
  regression.  It can also be used in unsupervised mode for locating
  outliers or assessing proximities among data points.
}
\usage{
\method{randomForest}{formula}(x, data=NULL, ..., subset, na.action=na.fail)
\method{randomForest}{default}(x, y = NULL, xtest = NULL, ytest = NULL, addclass = 0, 
    ntree = 500, mtry = ifelse(!is.null(y) && !is.factor(y), 
        max(floor(ncol(x)/3), 1), floor(sqrt(ncol(x)))), classwt = NULL, cutoff,
    nodesize = ifelse(!is.null(y) && !is.factor(y), 5, 1), importance = FALSE, 
    proximity = FALSE, outscale = FALSE, norm.votes = TRUE, do.trace = FALSE, 
    keep.forest = is.null(xtest), corr.bias=FALSE, ...) 
NULL

\method{print}{randomForest}(x, ...)
}
\arguments{
  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which
    \code{randomForest} is called from.}
  \item{subset}{an index vector indicating which rows should be used.
    (NOTE: If given, this argument must be named.)}
  \item{na.action}{A function to specify the action to be taken if NAs
    are found.  (NOTE: If given, this argument must be named.)}
  \item{x}{a data frame or a matrix of predictors, or a formula
    describing the model to be fitted (for the
    \code{print} method, an \code{randomForest} object).}
  \item{y}{A response vector.  If a factor, classification is assumed,
    otherwise regression is assumed.  If omitted, \code{randomForest}
    will run in unsupervised mode with \code{addclass=1} (unless
    explicitly set otherwise).}
  \item{xtest}{a data frame or matrix (like \code{x}) containing
    predictors for the test set.}
  \item{ytest}{response for the test set.}
  \item{addclass}{\code{=0} (default) do not add a synthetic class to
    the data. \code{=1} label the input data as class 1 and add a
    synthetic class by randomly sampling from the product of empirical
    marginal distributions of the input.  \code{=2} is similar to
    \code{=1}, but the synthetic data are sampled from the uniform
    hyperrectangle that contain the input.  Ignored for regression.}
  \item{ntree}{Number of trees to grow.  This should not be set to too
    small a number, to ensure that every input row gets predicted at
    least a few times. }
  \item{mtry}{Number of variables randomly sampled as candidates at each
    split.  Note that the default values are different for
    classification and regression}
  \item{classwt}{Priors of the classes.  Need not add up to one.
    Ignored for regression.}
  \item{cutoff}{(Classification only)  A vector of length equal to
    number of classes.  The `winning' class for an observation is the
    one with the maximum (positive) difference between proportion of
    votes and cutoff.  Default is 1/k where k is the number of classes
    (i.e., majority vote wins).}
  \item{nodesize}{Minimum size of terminal nodes.  Setting this number
    larger causes smaller trees to be grown (and thus take less time).
    Note that the default values are different for classification
    and regression.}
  \item{importance}{Should importance of predictors be assessed? }
  \item{proximity}{Should proximity measure among the rows be
    calculated?}
  \item{outscale}{Should outlyingness of rows be assessed?  Ignored for
    regression.}
  \item{norm.votes}{If \code{TRUE} (default), the final result of votes
    are expressed as fractions.  If \code{FALSE}, raw vote counts are
    returned (useful for combining results from different runs).
    Ignored for regression.}
  \item{do.trace}{If set to \code{TRUE}, give a more verbose output as
    \code{randomForest} is run.  If set to some integer, then running
    output is printed for every \code{do.trace} trees.}
  \item{keep.forest}{If set to \code{FALSE}, the forest will not be
    retained in the output object.  If \code{xtest} is given, defaults
    to \code{FALSE}.}
  \item{corr.bias}{perform bias correction for regression?  Note:
    Experimental.  Use at your own risk.}
  \item{...}{optional parameters to be passed to the low level function
    \code{randomForest.default}.}
}

\value{
  An object of class \code{randomForest}, which is a list with the
  following components:

  \item{call}{the original call to \code{randomForest}}
  \item{type}{one of \code{regression}, \code{classification}, or
    {unsupervised}.}
  \item{predicted}{the predicted values of the input data based on
    out-of-bag samples.}
  \item{importance}{a matrix with two columns, each a different measure
    of importance of the predictors.  If
    \code{importance=FALSE}, the last measure is still returned.  This
    measure corresponds to average decrease in node impurity caused by
    splitting on the variables, and is always computed.}
  \item{ntree}{number of trees grown.}
  \item{mtry}{number of predictors sampled for spliting at each node.}
  \item{forest}{(a list that contains the entire forest; \code{NULL} if 
    \code{randomForest} is run in unsupervised mode or if
    \code{keep.forest=FALSE}.}
  \item{err.rate}{(classification only) vector error rates of the
    prediction on the input data, the i-th element being the error rate
    for all trees up to the i-th.} 
  \item{confusion}{(classification only) the confusion matrix of the
    prediction.}
  \item{votes}{(classification only) a matrix with one row for each
    input data point and one column for each class, giving the fraction
    or number of `votes' from the random forest.}
  \item{proximity}{if \code{proximity=TRUE} when
    \code{randomForest} is called, a matrix of proximity measures among
    the input (based on the frequency that pairs of data points are in
    the same terminal nodes).}
  \item{outlier}{(classification only) if \code{outscale=TRUE} when
    \code{randomForest} is called, a vector indicating how outlying the
    data points are (based on the proximity measures).}

  \item{mse}{(regression only) vector of mean square errors: sum of squared
    residuals divided by \code{n}.}
  \item{rsq}{(regression only) ``pseudo R-squared'': 1 - \code{mse} /
    Var(y).} 
  \item{test}{if test set is given (through the \code{xtest} or additionally
  \code{ytest} arguments), this component is a list which contains the
  corresponding \code{predicted}, \code{err.rate}, \code{confusion},
  \code{votes} (for classification) or \code{predicted}, \code{mse} and
  \code{rsq} (for regression) for the test set.  If
  \code{proximity=TRUE}, there is also a component, \code{proximity},
  which contains the proximity among the test set as well as proximity
  between test and training data.} 
}

\note{
  The \code{forest} structure is slightly different between
  classification and regression.

  If \code{xtest} is given, prediction of the test set is done ``in
  place'' as the trees are grown.  If \code{ytest} is also given, and
  \code{do.trace} is set to some positive integer, then for every
  \code{do.trace} trees, the test set error is printed.  Results for the
  test set is returned in the \code{test} component of the resulting
  \code{randomForest} object.

  Here are the definitions of the variable importance measures.  For
  classification, define the margin of a case as the proportion of votes
  for the correct class minus the maximum proportion of votes for the
  other classes.  The first measure is the percent decrease in average
  margin when the out-of-bag portion of the data are permuted one
  variable at a time.  The second measure is the total decrease in Gini
  index in nodes split by the variable, averaged over all trees.  For
  regression, the first measure is the percent increase in mean squared
  residuals when the out-of-bag portion of the data are permuted one
  variable at a time.  The second measure is the total decrease in sum
  of squares in nodes split by the variable, averaged over all trees.
}

\references{
  Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1),
  5-32.

  Breiman, L (2002), ``Manual On Setting Up, Using, And Understanding
  Random Forests V3.1'', \url{http://oz.berkeley.edu/users/breiman/ Using_random_forests_V3.1.pdf}.
}
\author{Andy Liaw \email{andy\_liaw@merck.com} and Matthew Wiener
  \email{matthew\_wiener@merck.com}, based on original Fortran code by
  Leo Breiman and Adele Cutler.}

\seealso{\code{\link{predict.randomForest}}, \code{\link{var.imp.plot}}}

\examples{
## Classification:
data(iris)
set.seed(71)
iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE,
                        proximity=TRUE)
print(iris.rf)
## Look at variable importance:
print(round(iris.rf$importance, 2))
## Do MDS on 1 - proximity:
require(mva)
iris.mds <- cmdscale(1 - iris.rf$proximity, eig=TRUE)
op <- par(pty="s")
pairs(cbind(iris[,1:4], iris.mds$points), cex=0.6, gap=0, 
      col=c("red", "green", "blue")[codes(iris$Species)],
      main="Iris Data: Predictors and MDS of Proximity Based on RandomForest")
par(op)
print(iris.mds$GOF)

## The `unsupervised' case:
set.seed(17)
iris.urf <- randomForest(iris[, -5], proximity=TRUE, outscale=TRUE)
## Look for Outliers:
plot(iris.urf$out, type="h", ylab="",
     main="Measure of Outlyingness for Iris Data")

## Regression:
data(airquality)
set.seed(131)
ozone.rf <- randomForest(Ozone ~ ., data=airquality, mtry=3, importance=TRUE)
print(ozone.rf)
## Show "importance" of variables: higher value mean more important:
print(round(ozone.rf$importance, 2))
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{regression}
\keyword{tree}
