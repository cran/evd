\name{bvhr}
\alias{dbvhr}
\alias{pbvhr}
\alias{rbvhr}
\title{The Husler-Reiss Model}
\description{
  Density, distribution function and random generation for the
  Husler-Reiss model. 
}
\usage{
dbvhr(x, dep, mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvhr(q, dep, mar1 = c(0, 1, 0), mar2 = mar1) 
rbvhr(n, dep, mar1 = c(0, 1, 0), mar2 = mar1) 
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
  The Husler-Reiss distribution function (Husler and Reiss, 1989)
  with parameter \eqn{\code{dep} = r} is
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

  Independence is obtained in the limit as \eqn{r} approaches zero.
  Complete dependence is obtained as \eqn{r} tends to infinity.
}
\value{
  \code{dbvhr} gives the density, \code{pbvhr} gives the
  distribution function and \code{rbvhr} generates random deviates.
}
\references{ 
  Husler, J. and Reiss, R.-D. (1989)
  Maxima of normal random vectors: between independence 
  and complete dependence.
  \emph{Statist. Probab. Letters}, \bold{7}, 283--286.
}

\seealso{\code{\link{abvhr}}, \code{\link{rgev}}}

\examples{
dbvhr(matrix(rep(0:4,2),ncol=2), 1.7)
pbvhr(matrix(rep(0:4,2),ncol=2), 1.7)  
rbvhr(10, 1.7)
}
\keyword{distribution}

