\name{bvnegbilog}
\alias{dbvnegbilog}
\alias{pbvnegbilog}
\alias{rbvnegbilog}
\title{The Bivariate Negative Bilogistic Model}
\description{
  Density, distribution function and random generation for the
  bivariate negative bilogistic model. 
}
\usage{
dbvnegbilog(x, alpha, beta, mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvnegbilog(q, alpha, beta, mar1 = c(0, 1, 0), mar2 = mar1) 
rbvnegbilog(n, alpha, beta, mar1 = c(0, 1, 0), mar2 = mar1)  
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
  The negative bilogistic distribution function with
  parameters \eqn{\code{alpha} = \alpha}{\code{alpha} = alpha}
  and \eqn{\code{beta} = \beta}{\code{beta} = beta} is
  \deqn{G(z_1,z_2) = \exp\left\{- y_1 - y_2 + y_1 q^{1+\alpha} +
    y_2 (1-q)^{1+\beta}\right\}}{
    G(z1,z2) = exp{- y1 - y2 + y1 q^(1+alpha) + y2 (1-q)^(1+beta)}}
  where
  \eqn{q = q(y_1,y_2;\alpha,\beta)}{q = q(y1,y2;alpha,beta)}
  is the root of the equation
  \deqn{(1+\alpha) y_1 q^\alpha - (1+\beta) y_2 (1-q)^\beta = 0,}{
    (1+alpha) y1 q^alpha - (1+beta) y2 (1-q)^beta = 0,}
  \eqn{\alpha > 0}{alpha > 0}, \eqn{\beta > 0}{beta > 0} and   
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

  When \eqn{\alpha = \beta}{alpha = beta} the negative bilogistic
  model is equivalent to the negative logistic model with dependence
  parameter
  \eqn{\code{dep} = 1/\alpha = 1/\beta}{
    \code{dep} = 1/alpha = 1/beta}.
  Complete dependence is obtained in the limit as
  \eqn{\alpha = \beta}{alpha = beta} approaches zero.
  Independence is obtained as
  \eqn{\alpha = \beta}{alpha = beta} tends to infinity, and when
  one of \eqn{\alpha,\beta}{alpha,beta} is fixed and the other
  tends to infinity.
  Different limits occur when one of
  \eqn{\alpha,\beta}{alpha,beta} is fixed and the other
  approaches zero.
  The model was first introduced by Coles and Tawn (1994).
}
\value{
  \code{dbvnegbilog} gives the density, \code{pbvnegbilog} gives the
  distribution function and \code{rbvnegbilog} generates random deviates.
}
\note{
  \code{rbvnegbilog} uses a root finding algorithm to evaluate \eqn{q}
  within the root finding algorithm used to simulate from the
  conditional distribution.
  It is therefore relatively slow (about 2.8 seconds per 1000 random
  vectors on a 450MHz PIII, 512Mb RAM).

  The bilogistic and negative bilogistic models can be represented under
  a single model, using the integral of the maximum of two beta
  distributions (Joe, 1997).
}
\references{
  Coles, S. G. and Tawn, J. A. (1994)
  Statistical methods for multivariate extremes: an application to
  structural design (with discussion).
  \emph{Appl. Statist.}, \bold{43}, 1--48.
  
  Joe, H. (1997)
  \emph{Multivariate Models and Dependence Concepts},
  London: Chapman & Hall.
}

\seealso{\code{\link{abvnegbilog}}, \code{\link{rbvbilog}},
  \code{\link{rgev}}}

\examples{
dbvnegbilog(matrix(rep(0:4,2),ncol=2), .7, 1.52)
pbvnegbilog(matrix(rep(0:4,2),ncol=2), .7, 1.52)  
rbvnegbilog(10, .7, 1.52)
}
\keyword{distribution}

