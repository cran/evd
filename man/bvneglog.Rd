\name{bvneglog}
\alias{dbvneglog}
\alias{pbvneglog}
\alias{rbvneglog}
\title{The Bivariate Negative Logistic Distribution}
\description{
  Density, distribution function and random generation for the
  bivariate negative logistic distribution. 
}
\usage{
dbvneglog(x, dep, mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvneglog(q, dep, mar1 = c(0, 1, 0), mar2 = mar1) 
rbvneglog(n, dep, mar1 = c(0, 1, 0), mar2 = mar1) 
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
  The bivariate negative logistic distribution function
  with parameter \eqn{\code{dep} = r} is
  \deqn{G(z_1,z_2) = \exp\left\{-y_1-y_2+
    [y_1^{-r}+y_2^{-r}]^{-1/r}\right\}}{
    G(z1,z2) = exp{-y1-y2+[y1^(-r)+y2^(-r)]^(-1/r)}}
  where \eqn{r > 0} and   
  \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,2}, where the marginal
  parameters are given by
  \eqn{\code{mari} = (a_i,b_i,s_i)}{\code{mari} = (ai,bi,si)},
  \eqn{b_i > 0}{bi > 0}.
  If \eqn{s_i = 0}{si = 0} then \eqn{y_i}{yi} is defined by
  continuity.
  This is a special case of the bivariate asymmetric negative
  logistic distribution.
  The univariate marginal distributions are generalized extreme
  value.
}
\value{
  \code{dbvneglog} gives the density, \code{pbvneglog} gives the
  distribution function and \code{rbvneglog} generates random deviates.
}
\note{
  \code{rbvneglog} uses a root-finding algorithm for each realization,
  and is therefore relatively slow.
}

\seealso{\code{\link{abvneglog}}, \code{\link{rbvaneglog}},
  \code{\link{rgev}}}

\examples{
dbvneglog(matrix(rep(0:4,2),ncol=2), 1.2)
pbvneglog(matrix(rep(0:4,2),ncol=2), 1.2)  
rbvneglog(10, 1.2)
}
\keyword{distribution}

