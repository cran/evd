\name{gumbelx}
\alias{dgumbelx}
\alias{pgumbelx}
\alias{qgumbelx}
\alias{rgumbelx}
\title{Maxima of Two Gumbel Distributions}
\description{
  Density function, distribution function, quantile function and
  random generation for the maxima of two Gumbel distributions,
  each with different location and scale parameters. 
}
\usage{
dgumbelx(x, loc1=0, scale1=1, loc2=0, scale2=1, log = FALSE) 
pgumbelx(q, loc1=0, scale1=1, loc2=0, scale2=1, lower.tail = TRUE) 
qgumbelx(p, interval, loc1=0, scale1=1, loc2=0, scale2=1, lower.tail = TRUE, \dots)
rgumbelx(n, loc1=0, scale1=1, loc2=0, scale2=1)
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations.}
  \item{interval}{A length two vector containing the end-points of the interval 
    to be searched for the quantiles, passed to the uniroot function.}
  \item{loc1, scale1, loc2, scale2}{Location and scale parameters of the two 
    Gumbel distributions. The second location parameter must be greater than or
    equal to the first location parameter.}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
  \item{lower.tail}{Logical; if \code{TRUE} (default), probabilities
    are P[X <= x], otherwise, P[X > x]}
  \item{\dots}{Other arguments passed to uniroot.}
}
\value{
  \code{dgumbelx} gives the density function, \code{pgumbelx} gives the
  distribution function, \code{qgumbelx} gives the quantile function,
  and \code{rgumbelx} generates random deviates.
}	

\seealso{\code{\link{fgev}}, \code{\link{rfrechet}},
  \code{\link{rgumbel}}, \code{\link{rrweibull}}, \code{\link{uniroot}}}

\examples{
dgumbelx(2:4, 0, 1.1, 1, 0.5)
pgumbelx(2:4, 0, 1.1, 1, 0.5)
qgumbelx(seq(0.9, 0.6, -0.1), interval = c(0,10), 0, 1.2, 2, 0.5)
rgumbelx(6, 0, 1.1, 1, 0.5)
p <- (1:9)/10
pgumbelx(qgumbelx(p, interval = c(0,10), 0, 0.5, 1, 2), 0, 0.5, 1, 2)
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
}
\keyword{distribution}

