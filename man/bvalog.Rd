\name{bvalog}
\alias{dbvalog}
\alias{pbvalog}
\alias{rbvalog}
\title{The Bivariate Asymmetric Logistic Model}
\description{
  Density, distribution function and random generation for the
  bivariate asymmetric logistic model. 
}
\usage{
dbvalog(x, dep, asy = c(1, 1), mar1 = c(0, 1, 0), mar2 = mar1,
    log = FALSE) 
pbvalog(q, dep, asy = c(1, 1), mar1 = c(0, 1, 0), mar2 = mar1) 
rbvalog(n, dep, asy = c(1, 1), mar1 = c(0, 1, 0), mar2 = mar1) 
}
\arguments{
  \item{x, q}{A vector of length two or a matrix with two columns,
    in which case the density/distribution is evaluated across
    the rows.}
  \item{n}{Number of observations.}
  \item{dep}{Dependence parameter.}
  \item{asy}{A vector containing the two asymmetry parameters.}
  \item{mar1, mar2}{Vectors of length three containing marginal
    parameters, or matrices with three columns where each
    column represents a vector of values to be passed to the
    corresponding marginal parameter.}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
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

  When \eqn{t_1 = t_2 = 1}{t1 = t2 = 1} the asymmetric logistic
  model is equivalent to the logistic model.
  Independence is obtained when either \eqn{r = 1},
  \eqn{t_1 = 0}{t1 = 0} or \eqn{t_2 = 0}{t2 = 0}.
  Complete dependence is obtained in the limit when
  \eqn{t_1 = t_2 = 1}{t1 = t2 = 1} and \eqn{r}
  approaches zero.
  Different limits occur when \eqn{t_1}{t1} and \eqn{t_2}{t2}
  are fixed and \eqn{r} approaches zero.
  The model was first introduced by Tawn (1988). 
}
\value{
  \code{dbvalog} gives the density, \code{pbvalog} gives the
  distribution function and \code{rbvalog} generates random deviates.
}
\note{
  \code{rbvalog} uses the bivariate version of Algorithm 1.2 in
  Stephenson (2002).
}
\references{
  Stephenson, A. G. (2002)
  Simulating multivariate extreme value distributions of logistic type.
  To be published - available on request.

  Tawn, J. A. (1988)
  Bivariate extreme value theory: models and estimation.
  \emph{Biometrika}, \bold{75}, 397--415.  
}

\seealso{\code{\link{abvlog}}, \code{\link{rbvlog}},
  \code{\link{rgev}}, \code{\link{rmvalog}}}

\examples{
dbvalog(matrix(rep(0:4,2),ncol=2), .7, c(0.5,1))
pbvalog(matrix(rep(0:4,2),ncol=2), .7, c(0.5,1))  
rbvalog(10, .7, c(0.5,1))
}
\keyword{distribution}

