\name{frechet}
\alias{dfrechet}
\alias{pfrechet}
\alias{qfrechet}
\alias{rfrechet}
\title{The Frechet Distribution}
\description{
  Density function, distribution function, quantile function
  and random generation for the Frechet distribution with
  location, scale and shape parameters. 
}
\usage{
dfrechet(x, loc=0, scale=1, shape=1, log = FALSE) 
pfrechet(q, loc=0, scale=1, shape=1, lower.tail = TRUE) 
qfrechet(p, loc=0, scale=1, shape=1, lower.tail = TRUE)
rfrechet(n, loc=0, scale=1, shape=1)
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations.}
  \item{loc, scale, shape}{Location, scale and shape parameters (can be
    given as vectors).}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
  \item{lower.tail}{Logical; if \code{TRUE} (default), probabilities
    are P[X <= x], otherwise, P[X > x]}
}
\details{
  The Frechet distribution function with parameters
  \eqn{\code{loc} = a}, \eqn{\code{scale} = b} and
  \eqn{\code{shape} = s} is
  \deqn{G(z) = \exp\left\{-\left(\frac{z-a}{b}\right)^{-s}
  \right\}}{G(x) = exp{-[(z-a)/b]^(-s)}}
  for \eqn{z > a} and zero otherwise, where \eqn{b > 0} and
  \eqn{s > 0}.
}
\value{
  \code{dfrechet} gives the density function, \code{pfrechet} gives
  the distribution function, \code{qfrechet} gives the quantile
  function, and \code{rfrechet} generates random deviates.
}

\seealso{\code{\link{rgev}}, \code{\link{rgumbel}}, \code{\link{rrweibull}}}

\examples{
dfrechet(2:4, 1, 0.5, 0.8)
pfrechet(2:4, 1, 0.5, 0.8)
qfrechet(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
rfrechet(6, 1, 0.5, 0.8)
p <- (1:9)/10
pfrechet(qfrechet(p, 1, 2, 0.8), 1, 2, 0.8)
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
}
\keyword{distribution}

