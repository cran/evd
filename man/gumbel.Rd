\name{gumbel}
\alias{dgumbel}
\alias{pgumbel}
\alias{qgumbel}
\alias{rgumbel}
\title{The Gumbel Distribution}
\description{
  Density function, distribution function, quantile function and
  random generation for the Gumbel distribution with location and
  scale parameters. 
}
\usage{
dgumbel(x, loc=0, scale=1, log = FALSE) 
pgumbel(q, loc=0, scale=1, lower.tail = TRUE) 
qgumbel(p, loc=0, scale=1, lower.tail = TRUE)
rgumbel(n, loc=0, scale=1)
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations.}
  \item{loc, scale}{Location and scale parameters (can be given
    as vectors).}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
  \item{lower.tail}{Logical; if \code{TRUE} (default), probabilities
    are P[X <= x], otherwise, P[X > x]}
}
\details{
  The Gumbel distribution function with parameters
  \eqn{\code{loc} = a} and \eqn{\code{scale} = b} is
  \deqn{G(z) = \exp\left\{-\exp\left[-\left(\frac{z-a}{b}\right)
    \right]\right\}}{G(x) = exp{-exp[-(z-a)/b]}}
  for all real \eqn{z}, where \eqn{b > 0}.
}
\value{
  \code{dgumbel} gives the density function, \code{pgumbel} gives
  the distribution function, \code{qgumbel} gives the quantile
  function, and \code{rgumbel} generates random deviates.
}
\seealso{\code{\link{rfrechet}}, \code{\link{rgev}}, \code{\link{rrweibull}}}

\examples{
dgumbel(-1:2, -1, 0.5)
pgumbel(-1:2, -1, 0.5)
qgumbel(seq(0.9, 0.6, -0.1), 2, 0.5)
rgumbel(6, -1, 0.5)
p <- (1:9)/10
pgumbel(qgumbel(p, -1, 2), -1, 2)
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
}
\keyword{distribution}


