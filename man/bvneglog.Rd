\name{bvneglog}
\alias{dbvneglog}
\alias{pbvneglog}
\alias{rbvneglog}
\title{The Bivariate Negative Logistic Model}
\description{
  Density, distribution function and random generation for the
  bivariate negative logistic model. 
}
\usage{
dbvneglog(x, dep, mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvneglog(q, dep, mar1 = c(0, 1, 0), mar2 = mar1) 
rbvneglog(n, dep, mar1 = c(0, 1, 0), mar2 = mar1) 
}
\arguments{
  \item{x, q}{A vector of length two or a matrix with two columns,
    in which case the density/distribution is evaluated across
    the rows.}
  \item{n}{Number of observations.}
  \item{dep}{Dependence parameter.}
  \item{mar1, mar2}{Vectors of length three containing marginal
    parameters, or matrices with three columns where each
    column represents a vector of values to be passed to the
    corresponding marginal parameter.}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
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
  logistic model.
  The univariate marginal distributions are generalized extreme
  value.

  Independence is obtained in the limit as \eqn{r} approaches zero.
  Complete dependence is obtained as \eqn{r} tends to infinity.
  The earliest reference to this model appears to be
  Galambos (1975, Section 4).
}
\value{
  \code{dbvneglog} gives the density, \code{pbvneglog} gives the
  distribution function and \code{rbvneglog} generates random deviates.
}
\references{ 
  Galambos, J. (1975)
  Order statistics of samples from multivariate distributions.
  \emph{J. Amer. Statist. Assoc.}, \bold{70}, 674--680.
}

\seealso{\code{\link{abvneglog}}, \code{\link{rbvaneglog}},
  \code{\link{rgev}}}

\examples{
dbvneglog(matrix(rep(0:4,2),ncol=2), 1.2)
pbvneglog(matrix(rep(0:4,2),ncol=2), 1.2)  
rbvneglog(10, 1.2)
}
\keyword{distribution}

