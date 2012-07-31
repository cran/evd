\name{gpd}
\alias{dgpd}
\alias{pgpd}
\alias{qgpd}
\alias{rgpd}
\title{The Generalized Pareto Distribution}
\description{
  Density function, distribution function, quantile function and
  random generation for the generalized Pareto distribution (GPD)
  with location, scale and shape parameters. 
}
\usage{
dgpd(x, loc=0, scale=1, shape=0, log = FALSE) 
pgpd(q, loc=0, scale=1, shape=0, lower.tail = TRUE) 
qgpd(p, loc=0, scale=1, shape=0, lower.tail = TRUE)
rgpd(n, loc=0, scale=1, shape=0)
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations.}
  \item{loc, scale, shape}{Location, scale and shape parameters; the 
    \code{shape} argument cannot be a vector (must have length one).}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
  \item{lower.tail}{Logical; if \code{TRUE} (default), probabilities
    are P[X <= x], otherwise, P[X > x]}
}
\details{
  The generalized Pareto distribution function (Pickands, 1975) with
  parameters \eqn{\code{loc} = a}, \eqn{\code{scale} = b} and
  \eqn{\code{shape} = s} is
  \deqn{G(z) = 1 - \{1+s(z-a)/b\}^{-1/s}}{
    G(z) = 1 - {1+s(z-a)/b}^(-1/s)}
  for \eqn{1+s(z-a)/b > 0} and \eqn{z > a}, where \eqn{b > 0}.
  If \eqn{s = 0} the distribution is defined by continuity.
}
\value{
  \code{dgpd} gives the density function, \code{pgpd} gives the
  distribution function, \code{qgpd} gives the quantile function,
  and \code{rgpd} generates random deviates.
}
\references{ 
  Pickands, J. (1975)
  Statistical inference using extreme order statistics.
  \emph{Annals of Statistics}, \bold{3}, 119--131.
}	

\seealso{\code{\link{fpot}}, \code{\link{rgev}}}

\examples{
dgpd(2:4, 1, 0.5, 0.8)
pgpd(2:4, 1, 0.5, 0.8)
qgpd(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
rgpd(6, 1, 0.5, 0.8)
p <- (1:9)/10
pgpd(qgpd(p, 1, 2, 0.8), 1, 2, 0.8)
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
}
\keyword{distribution}

