\name{bvct}
\alias{dbvct}
\alias{pbvct}
\alias{rbvct}
\title{The Bivariate Coles-Tawn Model}
\description{
  Density, distribution function and random generation for the
  bivariate Coles-Tawn model (also known as the bivariate
  Dirichelet model). 
}
\usage{
dbvct(x, alpha, beta, mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvct(q, alpha, beta, mar1 = c(0, 1, 0), mar2 = mar1) 
rbvct(n, alpha, beta, mar1 = c(0, 1, 0), mar2 = mar1)  
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
  The Coles-Tawn (Coles and Tawn, 1991) distribution function with
  parameters \eqn{\code{alpha} = \alpha > 0}{\code{alpha} = alpha > 0}
  and \eqn{\code{beta} = \beta > 0}{\code{beta} = beta > 0} is
  \deqn{G(z_1,z_2) =
    \exp\left\{-y_1 [1 - \mbox{Be}(q;\alpha+1,\beta)] -
    y_2 \mbox{Be}(q;\alpha,\beta+1) \right\}}{
    G(z1,z2) = exp{- y1 [1 - Be(q;alpha+1,beta)] - y2 Be(q;alpha,beta+1)}}
  where
  \eqn{q = \alpha y_2 / (\alpha y_2 + \beta y_1)}{
    q = alpha y2 / (alpha y2 + beta y1)} and
  \eqn{\mbox{Be}(q;\alpha,\beta)}{Be(q;alpha,beta)} is the beta
  distribution function evaluated at \eqn{q} with
  \eqn{\code{shape1} = \alpha}{\code{shape1} = alpha},
  \eqn{\code{shape2} = \beta}{\code{shape2} = beta}, and where   
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

  Complete dependence is obtained in the limit as
  \eqn{\alpha = \beta}{alpha = beta} tends to infinity.
  Independence is obtained as
  \eqn{\alpha = \beta}{alpha = beta} approaches zero, and when
  one of \eqn{\alpha,\beta}{alpha,beta} is fixed and the other
  approaches zero.
  Different limits occur when one of
  \eqn{\alpha,\beta}{alpha,beta} is fixed and the other
  tends to infinity.
}
\value{
  \code{dbvct} gives the density, \code{pbvct} gives the
  distribution function and \code{rbvct} generates random deviates.
}

\references{
  Coles, S. G. and Tawn, J. A. (1991)
  Modelling extreme multivariate events.
  \emph{J. Roy. Statist. Soc., B}, \bold{53}, 377--392.
}

\seealso{\code{\link{abvct}}, \code{\link{rgev}}}

\examples{
dbvct(matrix(rep(0:4,2),ncol=2), .7, 0.52)
pbvct(matrix(rep(0:4,2),ncol=2), .7, 0.52)  
rbvct(10, .7, 0.52)
}
\keyword{distribution}

