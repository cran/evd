\name{mvlog}
\alias{pmvlog}
\alias{rmvlog}
\title{The Multivariate Logistic Model}
\description{
  Distribution function and random generation for the
  multivariate logistic model. 
}
\usage{
pmvlog(q, dep, d = 2, mar = c(0, 1, 0)) 
rmvlog(n, dep, d = 2, mar = c(0, 1, 0))
}
\arguments{
  \item{q}{A vector of length \code{d} or a matrix with \code{d}
    columns, in which case the distribution is evaluated across
    the rows.}
  \item{n}{Number of observations.}
  \item{dep}{Dependence parameter.}
  \item{d}{Dimension.}
  \item{mar}{A vector of length three containing marginal parameters
    for every univariate margin, or a matrix with three columns where
    each column represents a vector of values to be passed to the
    corresponding marginal parameter on every univariate margin.}
}
\details{
  Let \eqn{z = (z_1,z_2,\ldots,z_d)}{z = (z1,z2,\ldots,zd)}. 
  The \code{d} dimensional multivariate logistic distribution
  function (Gumbel, 1960) with parameter \eqn{\code{dep} = r} is
  \deqn{G(z) = \exp\left[-(y_1^{1/r}+\ldots+y_d^{1/r})^r\right]}{
    G(z) = exp{-[y1^(1/r)+\ldots+yd^(1/r)]^r}}
  where \eqn{0 < r \leq 1}{0 < r <= 1} and   
  \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,\ldots,d}{i = 1,\ldots,d}.
  Different parameters on each margin are not implemented, so
  \eqn{\code{mar} = (a_i,b_i,s_i)}{\code{mar} = (ai,bi,si)}
  for every \eqn{i}.
  If \eqn{s_i = 0}{si = 0} then \eqn{y_i}{yi} is defined by continuity.
  This is a special case of the multivariate asymmetric logistic
  model.
  The univariate marginal distributions are generalized extreme value.
}
\value{
  \code{pmvlog} gives the distribution function and \code{rmvlog}
  generates random deviates.
}
\note{
  The multivariate asymmetric logistic model
  (see \code{\link{rmvalog}}) provides a more
  flexible extension to the bivariate logistic model.
  Multivariate extensions of other bivariate models are more complex.
  A multivariate extension of the Husler-Reiss model exists, involving
  a multidimensional integral and one parameter for each bivariate
  margin.
  Multivariate extensions for the negative logistic model can be derived
  but are considerably more complex and appear to be less flexible. The
  obvious multivariate extension, which parallels that for the logistic
  model, is often called the multivariate negative logistic model
  (e.g Kotz \emph{et al}, 2000), but it is not a valid distribution
  function and should not be used.

  \code{rmvlog} uses Algorithm 2.1 in Stephenson(2002).
}
\references{
  Gumbel, E. J. (1960)
  Distributions des valeurs extremes en plusieurs dimensions.
  \emph{Publ. Inst. Statist. Univ. Paris}, \bold{9}, 171--173.
  
  Kotz, S. and Balakrishnan, N. and Johnson, N. L. (2000)
  \emph{Continuous Multivariate Distributions}, vol. 1.
  New York: John Wiley & Sons, 2nd edn.
  
  Stephenson, A. G. (2002)
  Simulating multivariate extreme value distributions of logistic type.
  To be published - available on request.
}

\seealso{\code{\link{rbvlog}}, \code{\link{rmvalog}},
  \code{\link{rgev}}}

\examples{
pmvlog(matrix(rep(0:4,5), ncol=5), .7, d = 5)
pmvlog(rep(4,5), .7, d = 5)
rmvlog(10, .7, d = 5)
}
\keyword{distribution}

