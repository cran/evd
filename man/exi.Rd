\name{exi}
\alias{exi}
\title{Estimates of the Extremal Index}
\description{
  Estimates of the extremal index. 
}
\usage{
exi(data, u, r = 1, ulow = rep(-Inf, ncol(u)), rlow =
    rep(1, length(r)), dimnames = list(NULL, NULL), drop = TRUE)
}
\arguments{
  \item{data}{A numeric vector, which may contain missing values.}
  \item{u}{A numeric vector of thresholds, unless time varying
    thresholds are used, in which case \code{u} should be a matrix,
    typically with \code{length(data)} rows (or else the usual
    recycling rules are applied).}
  \item{r}{A numeric vector of postive integers denoting the
    clustering interval lengths.}
  \item{ulow}{A numeric vector of lower thresholds with length
    \code{length(u)} (if \code{u} is a vector) or \code{ncol(u)}
    (if \code{u} is a matrix), unless time varying lower
    thresholds are used, in which case \code{ulow} should be a
    matrix with \code{length(u)} or \code{ncol(u)} columns, and
    typically with \code{length(data)} rows (or else the
    usual recycling rules are applied). By default there are no
    lower thresholds (or equivalently, the lower thresholds are
    \code{-Inf}).}
  \item{rlow}{A numeric vector of postive integers, of length
    \code{length(r)}, denoting the lower clustering interval
    lengths. The lower clustering interval length is only
    relevant if it is less than the clustering interval length
    and if there exists a lower threshold (greater than
    \code{-Inf}).}
  \item{dimnames}{The \code{dimnames} of the result.}
  \item{drop}{Logical; return a vector if either \code{u} or
    \code{r} has one element/column?}
}
\details{
  The extremal index is estimated using the inverse of the
  average cluster size, using the clusters of exceedences
  derived from \code{\link{clusters}}.
  If the threshold is larger than (or equal to)
  \code{max(data)}, so that no clusters are derived, then
  the corresponding estimate is \code{NaN}.
}
\value{
  A matrix with \code{length(u)} (if \code{u} is a vector) or
  \code{ncol(u)} (if \code{u} is a matrix) rows and
  \code{length(r)} columns, such that the \code{ij}th element
  is the estimate of the extremal index using the threshold
  \code{u[i]} or \code{u[,i]}, the clustering interval length
  \code{r[j]}, the lower threshold \code{ulow[i]} or
  \code{ulow[,i]}, and the lower clustering interval length
  \code{rlow[j]}.
  If \code{drop} is \code{TRUE} (the default), a numeric vector is
  returned in preference to a matrix with one row/column.
}

\seealso{\code{\link{clusters}}}

\examples{
data(portpirie)
exi(portpirie, 4.2, 3, 3.8)
us <- seq(3.9, 4.2, len = 10)
exi(portpirie, us, 3)
exi(portpirie, us, 1:3)
tvu <- c(rep(4.2, 20), rep(4.1, 25), rep(4.2, 20))
exi(portpirie, as.matrix(tvu), 1:3)

us <- seq(3.9, 4.2, len = 100)
eis <- exi(portpirie, us, 1:3)
matplot(us, eis, type = "l", xlab = "Threshold", ylab = "Ext. Index")
}
\keyword{manip}

