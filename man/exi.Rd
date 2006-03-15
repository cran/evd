\name{exi}
\alias{exi}
\title{Estimates of the Extremal Index}
\description{
  Estimates of the extremal index. 
}
\synopsis{
exi(data, u, r = 1, ulow = -Inf, rlow = 1)
}
\usage{
exi(data, u, r = 1, ulow = -Inf)
}
\arguments{
  \item{data}{A numeric vector, which may contain missing values.}
  \item{u}{A single value giving the threshold, unless a time varying
    threshold is used, in which case \code{u} should be a vector of
    thresholds, typically with the same length as \code{data} (or else
    the usual recycling rules are applied).}
  \item{r}{Either a postive integer denoting the clustering interval
    length, or zero, in which case the intervals estimator of
    Ferro and Segers (2003) is used and following arguments are
    ignored. By default the interval length is one.}
  \item{ulow}{A single value giving the lower threshold, unless a time
    varying lower threshold is used, in which case \code{ulow} should
    be a vector of lower thresholds, typically with the same length as
    \code{data} (or else the usual recycling rules are applied).
    By default there is no lower threshold (or equivalently, the
    lower threshold is \code{-Inf}).}
}
\details{
  If \code{r} is a positive integer the extremal index is estimated
  using the inverse of the average cluster size, using the clusters
  of exceedences derived from \code{\link{clusters}}. If \code{r} is
  zero, an estimate based on inter-exceedance times is used (Ferro
  and Segers, 2003).
  
  If there are no exceedances of the threshold, the estimate is
  \code{NaN}. If there is only one exceedance, the estimate is
  one.
}
\value{
  A single value estimating the extremal index.
}
\references{
  Ferro, C. A. T. and Segers, J. (2003)
  Inference for clusters of extreme values.
  \emph{JRSS B}, \bold{65}, 545--556.
}

\seealso{\code{\link{clusters}}, \code{\link{exiplot}}}

\examples{
data(portpirie)
exi(portpirie, 4.2, r = 3, ulow = 3.8)
tvu <- c(rep(4.2, 20), rep(4.1, 25), rep(4.2, 20))
exi(portpirie, tvu, r = 1)
exi(portpirie, tvu, r = 0)
}
\keyword{manip}

