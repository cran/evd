\name{mvalog}
\alias{pmvalog}
\alias{rmvalog}
\title{The Multivariate Asymmetric Logistic Distribution}
\description{
  Distribution function and random generation for the
  multivariate asymmetric logistic distribution. 
}
\usage{
pmvalog(q, dep, asy, d = 2, mar = c(1, 1, 0)) 
rmvalog(n, dep, asy, d = 2, mar = c(1, 1, 0))
}
\arguments{
  \item{q}{a vector of length \code{d} or a matrix with \code{d}
    columns, in which case the distribution is evaluated across
    the rows.}
  \item{n}{number of observations.}
  \item{dep}{a vector of length \eqn{2^d-d-1} containing the dependence
    parameters.}
  \item{asy}{a list with \eqn{2^d-1} vector elements containing
    asymmetry parameters for each separate component.}
  \item{d}{dimension.}
  \item{mar}{a vector of length three containing the marginal
    parameters for every univariate margin.}
}
\details{
  Let \eqn{z = (z_1,z_2,\ldots,z_d)}{z = (z1,z2,\ldots,zd)}.
  Let \eqn{B} be the set of all non-empty subsets of
  \eqn{\{1,\ldots,d\}}{{1,\ldots,d}}, let
  \eqn{B_1=\{b \in B:|b|=1\}}{B1={b in B:|b|=1}} and let
  \eqn{B_{(i)}=\{b \in B:i \in b\}}{B(i)={b in B:i in b}}.
  The \code{d} dimensional multivariate asymmetric logistic distribution
  function is
  \deqn{G(z)=\exp\left\{-\sum\nolimits_{b \in B} \left[\sum\nolimits_
    {i\in b}(t_{i,b}y_i)^{1/r_b}\right]^{r_b}\right\},}{G(z) =
    exp{-sum{b in B} [sum{i in b}(t{i,b}yi)^(1/r{b})]^r{b}},}
  where the dependence parameters \eqn{r_b\in(0,1]}{r{b} in (0,1]} for
  all \eqn{b\in B \setminus B_1}{b in B\B1}, and the asymmetry parameters
  \eqn{t_{i,b}\in[0,1]}{t{i,b} in [0,1]} for all \eqn{b\in B}{b in B}
  and \eqn{i\in b}{i in b}, and where 
  \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,\ldots,d}.
  If \eqn{s_i = 0}{si = 0} then \eqn{y_i}{yi} is defined by continuity.
  Different parameters on each margin are not implemented, so
  \eqn{\code{mar} = (a_i,b_i,s_i)}{\code{mar} = (ai,bi,si)}
  for every \eqn{i}.
  The constraints
  \eqn{\sum_{b \in B_{(i)}}t_{i,b}=1}{sum{b in B(i)} t{i,b}=1} for
  \eqn{i = 1,\ldots,d}
  ensure that the marginal distributions are generalized extreme value.
  More subtle constraints exist and can be found in the help manual.

  \code{dep} should be a vector of length \eqn{2^d-d-1} which contains
  \eqn{\{r_b:b\in B \setminus B_1\}}{{r{b}:b in B\B1}}, with
  the order defined by the natural set ordering on the index.
  In particular, for the trivariate model \code{dep} should be given
  as \eqn{(r_{12},r_{13},r_{23},r_{123})}{(r{12},r{13},r{23},r{123})}.
  \code{asy} should be a list with \eqn{2^d-1} elements.
  Each element is a vector which corresponds to a set
  \eqn{b\in B}{b in B}, containing \eqn{t_{i,b}}{t{i,b}} for
  every integer \eqn{i\in b}{i in b}.
  The elements should be given using the natural set ordering on the
  \eqn{b\in B}{b in B}, so that the first \eqn{d} elements are vectors
  of length one corresponding to the sets
  \eqn{\{1\},\ldots,\{d\}}{{1},\ldots,{d}}, and the last element is a
  a vector of length \eqn{d}, corresponding to the set
  \eqn{\{1,\ldots,d\}}{{1,\ldots,d}}.
  \code{asy} must be constructed to ensure that all constraints are
  satisfied or an error will occur.
  See the help manual for detailed examples.
  The univariate marginal distributions are generalized extreme value.
}
\value{
  \code{pmvalog} gives the distribution function and \code{rmvalog}
  generates random deviates.
}
\note{
  The interface for \code{rmvalog} and \code{pmvalog} is not
  particularly good. 
  (In particular, the set-up time is large in high dimensions even if a
  large number of the asymmetry parameters are zero.)
  If anybody uses \code{rmvalog} or \code{pmvalog} above about five
  dimensions it may inspire me to make improvements. Possibly.

  Multivariate extensions of other bivariate models are more complex.
  A multivariate extension of the Husler-Reiss model exists, involving
  a multidimensional integral and one parameter for each bivariate
  margin.
  Multivariate extensions for the negative logistic model can be derived
  but are considerably more complex and appear to be less flexible (the
  obvious multivariate extension, which parallels that for the logistic
  model, is often called the multivariate negative logistic distribution
  (e.g Kotz \emph{et al}, 2000) but it is not a valid distribution function and
  should not be used).

  \code{rmvalog} uses Algorithm 2.2 in Stephenson(2002).
}
\references{
  Kotz, S. and Balakrishnan, N. and Johnson, N. L. (2000)
  \emph{Continuous Multivariate Distributions}, vol. 1.
  New York: John Wiley & Sons, 2nd edn.
  
  Stephenson, A. G. (2002)
  Simulating multivariate extreme value distributions of logistic type.
  To be published - available on request.
}

\seealso{\code{\link{rbvalog}}, \code{\link{rmvlog}},
  \code{\link{rgev}}}

\examples{
pmvalog(c(2, 2, 2), dep = c(.6,.5,.8,.3),
  asy = list(.4, .1, .6, c(.3,.2), c(.1,.1), c(.4,.1), c(.2,.3,.2)),
  d = 3)
rmvalog(3, dep = c(.6,.5,.8,.3),
  asy = list(.4, .0, .6, c(.3,.2), c(.1,.1), c(.4,.1), c(.2,.4,.2)),
  d = 3)
set.seed(50)
rmvalog(3, dep = rep(.7,11),
  asy = list(0, 0, 0, 0, c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0),
  c(.2,.1,.2), c(.1,.1,.2), c(.3,.4,.1), c(.2,.2,.2), c(.4,.6,.2,.5)),
  d = 4)
set.seed(50)
rmvalog(3, dep = c(rep(1,6),rep(.7,5)),
  asy = list(0, 0, 0, 0, c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0),
  c(.2,.1,.2), c(.1,.1,.2), c(.3,.4,.1), c(.2,.2,.2), c(.4,.6,.2,.5)),
  d = 4)
}
\keyword{distribution}

