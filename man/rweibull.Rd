\name{rweibull}
\alias{drweibull}
\alias{prweibull}
\alias{qrweibull}
\alias{rrweibull}
\title{The Reversed Weibull Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the reversed Weibull distribution with location,
  scale and shape parameters. 
}
\usage{
drweibull(x, loc=0, scale=1, shape=1, log = FALSE) 
prweibull(q, loc=0, scale=1, shape=1, lower.tail = TRUE) 
qrweibull(p, loc=0, scale=1, shape=1, lower.tail = TRUE)
rrweibull(n, loc=0, scale=1, shape=1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{loc, scale, shape}{location, scale and shape parameters.}
  \item{log}{logical; if TRUE, the log density is returned.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    P[X <= x], otherwise, P[X > x]}
}
\details{
  The reversed Weibull distribution function with parameters
  \eqn{\code{loc} = a}, \eqn{\code{scale} = b} and
  \eqn{\code{shape} = s} is
  \deqn{G(z) = \exp\left\{-\left[-\left(\frac{z-a}{b}\right)
    \right]^s\right\}}{G(x) = exp{-[-(z-a)/b]^s}}
  for \eqn{z < a} and one otherwise, where \eqn{b > 0} and
  \eqn{s > 0}.
}
\note{
  Within extreme value theory the reversed Weibull distibution is
  usually referred to as the Weibull distribution.
  I make a distinction to avoid confusion with the distribution used
  in survival analysis, which is related to the distribution given
  above by a change in sign.
}
\value{
  \code{drweibull} gives the density, \code{prweibull} gives the
  distribution function, \code{qrweibull} gives the quantile function,
  and \code{rrweibull} generates random deviates.
}

\seealso{\code{\link{rfrechet}}, \code{\link{rgev}}, \code{\link{rgumbel}}}

\examples{
drweibull(-5:-3, -1, 0.5, 0.8)
prweibull(-5:-3, -1, 0.5, 0.8)
qrweibull(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
rrweibull(6, -1, 0.5, 0.8)
p <- (1:9)/10
prweibull(qrweibull(p, -1, 2, 0.8), -1, 2, 0.8)
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
}
\keyword{distribution}




