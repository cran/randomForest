\name{getTree}
\alias{getTree}
\title{Extract a single tree from a forest.}
\description{
  This function extract the structure of a tree from a
  \code{randomForest} object. 
}
\usage{
getTree(rfobj, k=1)
}
\arguments{
  \item{rfobj}{a \code{\link{randomForest}} object.}
  \item{k}{which tree to extract?}
}
\value{
  A matrix with six columns and number of rows equal to total number of
  nodes in the tree.  The six columns are:
  \item{left daughter}{the row where the left daughter node is; 0 if the
  node is terminal}
  \item{right daughter}{the row where the right daughter node is; 0 if
    the node is terminal}
  \item{split var}{which variable was used to split the node; 0 if the
    node is terminal}
  \item{split point}{where the best split is; see Details for
    categorical predictor}
  \item{status}{is the node terminal (-1) or not (1)}
  \item{prediction}{the prediction for the node; 0 if the node is not
    terminal}
}
\details{
  For categorical predictors, the splitting point is represented by an
  integer, whose binary expansion gives the identities of the categories
  that goes to left or right.  For example, if a predictor has three
  categories, and the split point is 5.  The binary expansion of 5 is
  (1, 0, 1) (because \eqn{5 = 1*2^0 + 0*2^1 + 1*2^2}), so cases with
  categories 1 or 3 in this predictor get sent to the left, and the rest
  to the right.
}
%\references{
%}
\seealso{
  \code{\link{randomForest}}
}
\examples{
data(iris)
getTree(randomForest(Species ~ ., iris, ntree=10), 3)
}
\author{Andy Liaw \email{andy\_liaw@merck.com}}
\keyword{tree}
