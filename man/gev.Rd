\name{gev}
\alias{dgev}
\alias{pgev}
\alias{qgev}
\alias{rgev}
\title{The Generalized Extreme Value Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the generalized extreme value (GEV) distribution
  with location, scale and shape parameters. 
}
\usage{
dgev(x, loc=0, scale=1, shape=0, log = FALSE) 
pgev(q, loc=0, scale=1, shape=0, lower.tail = TRUE) 
qgev(p, loc=0, scale=1, shape=0, lower.tail = TRUE)
rgev(n, loc=0, scale=1, shape=0)
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations.}
  \item{loc, scale, shape}{Location, scale and shape parameters; the 
    \code{shape} argument cannot be a vector (must have length one).}
  \item{log}{Logical; if TRUE, the log density is returned.}
  \item{lower.tail}{Logical; if TRUE (default), probabilities are
    P[X <= x], otherwise, P[X > x]}
}
\details{
  The GEV distribution function with parameters
  \eqn{\code{loc} = a}, \eqn{\code{scale} = b} and
  \eqn{\code{shape} = s} is
  \deqn{G(z) = \exp\left[-\{1+s(z-a)/b\}^{-1/s}\right]}{
    G(x) = exp[-{1+s(z-a)/b}^(-1/s)]}
  for \eqn{1+s(z-a)/b > 0}, where \eqn{b > 0}.
  If \eqn{s = 0} the distribution is defined by continuity.
  The parametric form of the GEV encompasses that of the Gumbel,
  Frechet and reversed Weibull distributions, which are obtained
  for \eqn{s = 0}, \eqn{s > 0} and \eqn{s < 0} respectively.
  It was first introduced by Jenkinson (1955).
}
\value{
  \code{dgev} gives the density, \code{pgev} gives the
  distribution function, \code{qgev} gives the quantile function,
  and \code{rgev} generates random deviates.
}
\references{ 
  Jenkinson, A. F. (1955)
  The frequency distribution of the annual maximum (or minimum) of
  meteorological elements.
  \emph{Quart. J. R. Met. Soc.}, \bold{81}, 158--171.
}	

\seealso{\code{\link{rfrechet}}, \code{\link{rgumbel}}, \code{\link{rrweibull}}}

\examples{
dgev(2:4, 1, 0.5, 0.8)
pgev(2:4, 1, 0.5, 0.8)
qgev(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
rgev(6, 1, 0.5, 0.8)
p <- (1:9)/10
pgev(qgev(p, 1, 2, 0.8), 1, 2, 0.8)
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
}
\keyword{distribution}

