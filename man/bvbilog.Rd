\name{bvbilog}
\alias{dbvbilog}
\alias{pbvbilog}
\alias{rbvbilog}
\title{The Bivariate Bilogistic Model}
\description{
  Density, distribution function and random generation for the
  bivariate bilogistic model. 
}
\usage{
dbvbilog(x, alpha, beta, mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvbilog(q, alpha, beta, mar1 = c(0, 1, 0), mar2 = mar1) 
rbvbilog(n, alpha, beta, mar1 = c(0, 1, 0), mar2 = mar1)  
}
\arguments{
  \item{x, q}{A vector of length two or a matrix with two columns,
    in which case the density/distribution is evaluated across
    the rows.}
  \item{n}{Number of observations.}
  \item{alpha}{Alpha parameter.}
  \item{beta}{Beta parameter.}
  \item{mar1, mar2}{Vectors of length three containing marginal
    parameters, or matrices with three columns where each
    column represents a vector of values to be passed to the
    corresponding marginal parameter.}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
}
\details{
  The bilogistic distribution function with
  parameters \eqn{\code{alpha} = \alpha}{\code{alpha} = alpha}
  and \eqn{\code{beta} = \beta}{\code{beta} = beta} is
  \deqn{G(z_1,z_2) = \exp\left\{-y_1 q^{1-\alpha} -
    y_2 (1-q)^{1-\beta}\right\}}{
    G(z1,z2) = exp{- y1 q^(1-alpha) - y2 (1-q)^(1-beta)}}
  where
  \eqn{q = q(y_1,y_2;\alpha,\beta)}{q = q(y1,y2;alpha,beta)}
  is the root of the equation
  \deqn{(1-\alpha) y_1 (1-q)^\beta - (1-\beta) y_2 q^\alpha = 0,}{
    (1-alpha) y1 (1-q)^beta - (1-beta) y2 q^alpha = 0,}
  \eqn{0 < \alpha,\beta < 1}{0 < alpha,beta < 1} and   
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

  When \eqn{\alpha = \beta}{alpha = beta} the bilogistic model
  is equivalent to the logistic model with dependence parameter
  \eqn{\code{dep} = \alpha = \beta}{\code{dep} = alpha = beta}.
  Complete dependence is obtained in the limit as
  \eqn{\alpha = \beta}{alpha = beta} approaches zero.
  Independence is obtained as
  \eqn{\alpha = \beta}{alpha = beta} approaches one, and when
  one of \eqn{\alpha,\beta}{alpha,beta} is fixed and the other
  approaches one.
  Different limits occur when one of
  \eqn{\alpha,\beta}{alpha,beta} is fixed and the other
  approaches zero.
  A bilogistic model is fitted in Smith (1990), where it appears
  to have been first introduced. 
}
\value{
  \code{dbvbilog} gives the density, \code{pbvbilog} gives the
  distribution function and \code{rbvbilog} generates random deviates.
}
\note{
  \code{rbvbilog} uses a root finding algorithm to evaluate \eqn{q}
  within the root finding algorithm used to simulate from the
  conditional distribution.
  It is therefore relatively slow (about 2.7 seconds per 1000 random
  vectors on a 450MHz PIII, 512Mb RAM).

  The bilogistic and negative bilogistic models can be represented under
  a single model, using the integral of the maximum of two beta
  distributions (Joe, 1997).
}
\references{ 
  Joe, H. (1997)
  \emph{Multivariate Models and Dependence Concepts},
  London: Chapman & Hall.

  Smith, R. L. (1990)
  Extreme value theory. In
  \emph{Handbook of Applicable Mathematics} (ed. W. Ledermann),
  vol. 7. Chichester: John Wiley, pp. 437--471. 
}

\seealso{\code{\link{abvbilog}}, \code{\link{rbvnegbilog}},
  \code{\link{rgev}}}

\examples{
dbvbilog(matrix(rep(0:4,2),ncol=2), .7, 0.52)
pbvbilog(matrix(rep(0:4,2),ncol=2), .7, 0.52)  
rbvbilog(10, .7, 0.52)
}
\keyword{distribution}

