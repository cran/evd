\name{exiplot}
\alias{exiplot}
\title{Plot Estimates of the Extremal Index}
\description{
  Plots estimates of the extremal index. 
}
\usage{
exiplot(data, tlim, r = 1, ulow = -Inf, rlow = 1, add = FALSE, 
    nt = 100, lty = 1, xlab = "Threshold", ylab = "Ext. Index",
    ylim = c(0,1), \dots)
}
\arguments{
  \item{data}{A numeric vector, which may contain missing values.}
  \item{tlim}{A numeric vector of length two, giving the limits for
    the (time invariant) thresholds at which the estimates are
    evaluated.}
  \item{r, ulow, rlow}{The estimation method. See \code{\link{exi}}.}
  \item{add}{Add to an existing plot?}
  \item{nt}{The number of thresholds at which the estimates are
    evaluated.}
  \item{lty}{Line type.}     
  \item{xlab, ylab}{x and y axis labels.}
  \item{ylim}{y axis limits.} 
  \item{\dots}{Other arguments passed to \code{plot} or
    \code{lines}.}
}
\details{
  The estimates are calculated using the function \code{\link{exi}}.
}
\value{
  A list with components \code{x} and \code{y} is invisibly returned.
  The first component contains the thresholds, the second contains the
  estimates.
}

\seealso{\code{\link{clusters}}, \code{\link{exi}}}

\examples{
sdat <- mar(100, psi = 0.5)
tlim <- quantile(sdat, probs = c(0.4,0.9))
exiplot(sdat, tlim)
exiplot(sdat, tlim, r = 4, add = TRUE, lty = 2)
exiplot(sdat, tlim, r = 0, add = TRUE, lty = 4)
}
\keyword{hplot}


