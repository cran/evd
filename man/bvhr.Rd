\name{bvhr}
\alias{dbvhr}
\alias{pbvhr}
\alias{rbvhr}
\title{The Husler-Reiss Distribution}
\description{
  Density, distribution function and random generation for the
  Husler-Reiss distribution. 
}
\usage{
dbvhr(x, dep, mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvhr(q, dep, mar1 = c(0, 1, 0), mar2 = mar1) 
rbvhr(n, dep, mar1 = c(0, 1, 0), mar2 = mar1) 
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
  The Husler-Reiss distribution function with parameter
  \eqn{\code{dep} = r} is
  \deqn{G(z_1,z_2) = \exp\left(-y_1\Phi\{r^{-1}+{\textstyle\frac{1}{2}}
    r[\log(y_1/y_2)]\} - y_2\Phi\{r^{-1}+{\textstyle\frac{1}{2}}r
    [\log(y_2/y_1)]\}\right)}{
    G(z1,z2) = exp(-y1 Phi{r^{-1}+r[log(y1/y2)]/2} -
    y2 Phi{r^{-1}+r[log(y2/y1)]/2}}
  where \eqn{\Phi(\cdot)}{Phi()} is the standard normal distribution
  function, \eqn{r > 0} and   
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
  \code{dbvhr} gives the density, \code{pbvhr} gives the
  distribution function and \code{rbvhr} generates random deviates.
}
\note{
  \code{rbvhr} uses a root-finding algorithm for each realization,
  and is therefore relatively slow.
}

\seealso{\code{\link{abvhr}}, \code{\link{rgev}}}

\examples{
dbvhr(matrix(rep(0:4,2),ncol=2), 1.7)
pbvhr(matrix(rep(0:4,2),ncol=2), 1.7)  
rbvhr(10, 1.7)
}
\keyword{distribution}

