\name{bvalog}
\alias{dbvalog}
\alias{pbvalog}
\alias{rbvalog}
\title{The Bivariate Asymmetric Logistic Distribution}
\description{
  Density, distribution function and random generation for the
  bivariate asymmetric logistic distribution. 
}
\usage{
dbvalog(x, dep, asy, mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvalog(q, dep, asy, mar1 = c(0, 1, 0), mar2 = mar1) 
rbvalog(n, dep, asy, mar1 = c(0, 1, 0), mar2 = mar1) 
}
\arguments{
  \item{x, q}{a vector of length two or a matrix with two columns,
    in which case the density/distribution is evaluated across
    the rows.}
  \item{n}{number of observations.}
  \item{dep}{dependence parameter.}
  \item{asy}{a vector containing the two asymmetry parameters.}
  \item{mar1, mar2}{vectors of length three containing marginal
    parameters.}
  \item{log}{logical; if \code{TRUE}, the log density is returned.}
}
\details{
  The bivariate asymmetric logistic distribution function with
  parameters \eqn{\code{dep} = r} and
  \eqn{\code{asy} = (t_1,t_2)}{\code{asy} = (t1,t2)} is
  \deqn{G(z_1,z_2) = \exp\left\{-(1-t_1)y_1-(1-t_2)y_2-
    [(t_1y_1)^{1/r}+(t_2y_2)^{1/r}]^r\right\}}{
    G(z1,z2) = exp{-(1-t1)y1-(1-t2)y2-[(t1y1)^(1/r)+(t2y2)^(1/r)]^r}}
  where \eqn{0 < r \leq 1}{0 < r <= 1},
  \eqn{0 \leq t_1,t_2 \leq 1}{0 <= t1,t2 <= 1}, and   
  \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,2}, where the marginal
  parameters are given by
  \eqn{\code{mari} = (a_i,b_i,s_i)}{\code{mari} = (ai,bi,si)},
  \eqn{b_i > 0}{bi > 0}.
  If \eqn{s_i = 0}{si = 0} then \eqn{y_i}{yi} is defined by
  continuity.
  The univariate marginal distributions are generalized extreme
  value.
}
\value{
  \code{dbvalog} gives the density, \code{pbvalog} gives the
  distribution function and \code{rbvalog} generates random deviates.
}
\note{
  \code{rbvalog} uses the bivariate version of Algorithm 1.2 in
  Stephenson(2002).
}
\references{ 
  Stephenson, A. G. (2002)
  Simulating multivariate extreme value distributions of logistic type.
  To be published - available on request.
}

\seealso{\code{\link{abvlog}}, \code{\link{rbvlog}},
  \code{\link{rgev}}, \code{\link{rmvalog}}}

\examples{
dbvalog(matrix(rep(0:4,2),ncol=2), .7, c(0.5,1))
pbvalog(matrix(rep(0:4,2),ncol=2), .7, c(0.5,1))  
rbvalog(10, .7, c(0.5,1))
}
\keyword{distribution}

