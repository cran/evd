\name{bvaneglog}
\alias{dbvaneglog}
\alias{pbvaneglog}
\alias{rbvaneglog}
\title{The Bivariate Asymmetric Negative Logistic Distribution}
\description{
  Density, distribution function and random generation for the
  bivariate asymmetric negative logistic distribution. 
}
\usage{
dbvaneglog(x, dep, asy, mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvaneglog(q, dep, asy, mar1 = c(0, 1, 0), mar2 = mar1) 
rbvaneglog(n, dep, asy, mar1 = c(0, 1, 0), mar2 = mar1) 
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
  The bivariate asymmetric negative logistic distribution function
  with parameters parameters \eqn{\code{dep} = r} and
  \eqn{\code{asy} = (t_1,t_2)}{\code{asy} = (t1,t2)} is
  \deqn{G(z_1,z_2) = \exp\left\{-y_1-y_2+
    [(t_1y_1)^{-r}+(t_2y_2)^{-r}]^{-1/r}\right\}}{
    G(z1,z2) = exp{-y1-y2+[(t1y1)^(-r)+(t2y2)^(-r)]^(-1/r)}}
  where \eqn{r > 0}, \eqn{0 \leq t_1,t_2 \leq 1}{0 <= t1,t2 <= 1}, and   
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
  \code{dbvaneglog} gives the density, \code{pbvaneglog} gives the
  distribution function and \code{rbvaneglog} generates random deviates.
}
\note{
  \code{rbvaneglog} uses a root-finding algorithm for each realization,
  and is therefore relatively slow.
}

\seealso{\code{\link{abvneglog}}, \code{\link{rbvneglog}},
  \code{\link{rgev}}}

\examples{
dbvaneglog(matrix(rep(0:4,2),ncol=2), 1.2, c(0.5,1))
pbvaneglog(matrix(rep(0:4,2),ncol=2), 1.2, c(0.5,1))  
rbvaneglog(10, 1.2, c(0.5,1))
}
\keyword{distribution}

