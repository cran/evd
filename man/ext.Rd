\name{ext}
\alias{dext}
\alias{pext}
\alias{qext}
\alias{rext}
\title{Distributions of Maxima and Minima}
\description{
  Density, distribution function, quantile function and random
  generation for the maximum/minimum of a given number of
  independent variables from a specified distribution. 
}
\usage{
dext(x, densfun, distnfun, \dots, distn, mlen = 1, largest = TRUE,
     log = FALSE)
pext(q, distnfun, \dots, distn, mlen = 1, largest = TRUE, lower.tail = TRUE) 
qext(p, quantfun, \dots, distn, mlen = 1, largest = TRUE, lower.tail = TRUE) 
rext(n, quantfun, \dots, distn, mlen = 1, largest = TRUE) 
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations.}
  \item{densfun, distnfun, quantfun}{Density, distribution and
    quantile function of the distribution.}
  \item{\dots}{Parameters of the distribution.}
  \item{distn}{A character string, optionally specified as an
    alternative to \code{densfun}, \code{distnfun} and \code{quantfun}
    such that the density, distribution and quantile functions are
    formed upon the addition of the prefixes \code{d}, \code{p} and
    \code{q} respectively.}
  \item{mlen}{The number of independent variables over which the
    maximum/minimum is taken.}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
  \item{largest}{Logical; if \code{TRUE} (default) use maxima,
    otherwise minima.}
  \item{lower.tail}{Logical; if \code{TRUE} (default) probabilities
    are P[X <= x], otherwise  P[X > x].}
}
\details{
  The density function must have a \code{log} argument (a simple
  wrapper can always be constructed to achieve this).
  If any parameter arguments are missing the default values
  (where they exist) in the density, distribution and quantile
  functions are used.
}
\value{
  \code{dext} gives the density, \code{pext} gives the
  distribution function and \code{qext} gives the quantile function
  of the maximum/minimum of \code{mlen} independent variables from
  a specified distibution. \code{rext} generates random deviates.
}
\note{
  The simulation speed is invariant to \code{mlen} since \code{rext}
  does not simulate directly from the specified distribution.
}

\seealso{\code{\link{rgev}},\code{\link{rorder}}}

\examples{
dext(2:4, dnorm, pnorm, mean = 0.5, sd = 1.2, mlen = 5)
dext(2:4, distn = "norm", mean = 0.5, sd = 1.2, mlen = 5)
dext(2:4, distn = "exp", mlen = 2)
pext(2:4, distn = "exp", rate = 1.2, mlen = 2)
qext(seq(0.9, 0.6, -0.1), distn = "exp", rate = 1.2, mlen = 2)
rext(5, qgamma, shape = 1, mlen = 10)
p <- (1:9)/10
pexp(qext(p, distn = "exp", rate = 1.2, mlen = 1), rate = 1.2)
## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
}
\keyword{distribution}


