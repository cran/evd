\name{bvlog}
\alias{dbvlog}
\alias{pbvlog}
\alias{rbvlog}
\title{The Bivariate Logistic Distribution}
\description{
  Density, distribution function and random generation for the
  bivariate logistic distribution. 
}
\usage{
dbvlog(x, dep, mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvlog(q, dep, mar1 = c(0, 1, 0), mar2 = mar1) 
rbvlog(n, dep, mar1 = c(0, 1, 0), mar2 = mar1) 
}
\arguments{
  \item{x, q}{a vector of length two or a matrix with two columns,
    in which case the density/distribution is evaluated across
    the rows.}
  \item{n}{number of observations.}
  \item{dep}{dependence parameter.}
  \item{mar1, mar2}{vectors of length three containing marginal
    parameters.}
  \item{log}{logical; if \code{TRUE}, the log density is returned.}
}
\details{
  The bivariate logistic distribution function with
  parameter \eqn{\code{dep} = r} is
  \deqn{G(z_1,z_2) = \exp\left[-(y_1^{1/r}+y_2^{1/r})^r\right]}{
    G(z1,z2) = exp{-[y1^(1/r)+y2^(1/r)]^r}}
  where \eqn{0 < r \leq 1}{0 < r <= 1} and   
  \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,2}, where the marginal
  parameters are given by
  \eqn{\code{mari} = (a_i,b_i,s_i)}{\code{mari} = (ai,bi,si)},
  \eqn{b_i > 0}{bi > 0}.
  If \eqn{s_i = 0}{si = 0} then \eqn{y_i}{yi} is defined by
  continuity.
  This is a special case of the bivariate asymmetric logistic
  distribution.
  The univariate marginal distributions are generalized extreme
  value.
}
\value{
  \code{dbvlog} gives the density, \code{pbvlog} gives the
  distribution function and \code{rbvlog} generates random deviates.
}
\note{
  \code{rbvlog} uses the bivariate version of Algorithm 1.1 in
  Stephenson(2002).
}
\references{ 
  Stephenson, A. G. (2002)
  Simulating multivariate extreme value distributions of logistic type.
  To be published - available on request.
}

\seealso{\code{\link{abvlog}}, \code{\link{rbvalog}},
  \code{\link{rgev}}, \code{\link{rmvlog}}}

\examples{
dbvlog(matrix(rep(0:4,2),ncol=2), .7)
pbvlog(matrix(rep(0:4,2),ncol=2), .7)  
rbvlog(10, .7)
}
\keyword{distribution}

