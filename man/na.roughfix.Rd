\name{na.roughfix}
\alias{na.roughfix}
\alias{na.roughfix.default}
\alias{na.roughfix.data.frame}
\title{Rough Imputation of Missing Values}
\description{
  Impute Missing Values by median/mode.
}
\usage{
na.roughfix(object, ...)
}
\arguments{
  \item{object}{a data frame or numeric matrix.}
  \item{\dots}{further arguments special methods could require.}
}
\value{
  A completed data matrix or data frame.  For numeric variables,
  \code{NA}s are replaced with column medians.  For factor variables,
  \code{NA}s are replaced with the most frequent levels (breaking ties
  at random).  If \code{object} contains no \code{NA}s, it is returned
  unaltered. 
}
\notes{
  This is used as a starting point for imputing missing values by random
  forest.
}
%\references{
%}
\seealso{
\code{\link{rfImpute}}, \code{\link{randomForest}}.
}
\examples{
}
\author{Andy Liaw}
\keyword{NA}
