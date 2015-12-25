\name{mvevd}
\alias{pmvevd}
\alias{rmvevd}
\alias{dmvevd}
\title{Parametric Multivariate Extreme Value Distributions}
\description{
  Density function, distribution function and random generation
  for the multivariate logistic and multivariate asymmetric
  logistic models.
}
\usage{
pmvevd(q, dep, asy, model = c("log", "alog"), d = 2, mar = c(0,1,0),
    lower.tail = TRUE)
rmvevd(n, dep, asy, model = c("log", "alog"), d = 2, mar = c(0,1,0))
dmvevd(x, dep, asy, model = c("log", "alog"), d = 2, mar = c(0,1,0),
    log = FALSE) 
}
\arguments{
  \item{x, q}{A vector of length \code{d} or a matrix with \code{d}
    columns, in which case the density/distribution is evaluated
    across the rows.}
  \item{n}{Number of observations.}
  \item{dep}{The dependence parameter(s). For the logistic model,
    should be a single value. For the asymmetric logistic model,
    should be a vector of length \eqn{2^d-d-1}, or a single value,
    in which case the value is used for each of the \eqn{2^d-d-1}
    parameters (see \bold{Details}).}
  \item{asy}{The asymmetry parameters for the asymmetric logistic
    model. Should be a list with \eqn{2^d-1} vector elements
    containing the asymmetry parameters for each separate component
    (see \bold{Details}).}
  \item{model}{The specified model; a character string. Must be either
    \code{"log"} (the default) or \code{"alog"} (or any unique partial
    match), for the logistic and asymmetric logistic models respectively.}
  \item{d}{The dimension.}
  \item{mar}{A vector of length three containing marginal parameters
    for every univariate margin, or a matrix with three columns where
    each column represents a vector of values to be passed to the
    corresponding marginal parameter. It can also be a list with
    \code{d} elements, such that each element is either a vector of
    length three or a matrix with three columns, in which case the
    \eqn{i}th element represents the marginal parameters on the
    \eqn{i}th margin.}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
  \item{lower.tail}{Logical; if \code{TRUE} (default), the
    distribution function is returned; the survivor function
    is returned otherwise.}
}
\details{
  Define
  \deqn{y_i = y_i(z_i) = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = yi(zi) = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,\ldots,d}{i = 1,\ldots,d}, where the marginal
  parameters are given by
  \eqn{(a_i,b_i,s_i)}{(ai,bi,si)}, \eqn{b_i > 0}{bi > 0}.
  If \eqn{s_i = 0}{si = 0} then \eqn{y_i}{yi} is defined by
  continuity.
  Let \eqn{z = (z_1,z_2,\ldots,z_d)}{z = (z1,z2,\ldots,zd)}.
  In each of the multivariate distributions functions
  \eqn{G(z)} given below, the
  univariate margins are generalized extreme value, so that
  \eqn{G(z_i) = \exp(-y_i)}{G(zi) = \exp(-yi)} for
  \eqn{i = 1,\ldots,d}{i = 1,\ldots,d}.
  If \eqn{1+s_i(z_i-a_i)/b_i \leq 0}{1+si(zi-ai)/bi <= 0} for some
  \eqn{i = 1,\ldots,d}{i = 1,\ldots,d}, the value \eqn{z_i}{zi} is
  either greater than the upper end point (if \eqn{s_i < 0}{si < 0}),
  or less than the lower end point (if \eqn{s_i > 0}{si > 0}), of the
  \eqn{i}th univariate marginal distribution.
  
  \code{model = "log"} (Gumbel, 1960)
  
  The \code{d} dimensional multivariate logistic distribution
  function with parameter \eqn{\code{dep} = r} is
  \deqn{G(z) = \exp\left\{-\left(\sum\nolimits_{i = 1}^{d}
    y_i^{1/r}\right)^r\right\}}{
    G(z) = exp{-[sum_{i=1}^d yi^(1/r)]^r}}
  where \eqn{0 < r \leq 1}{0 < r <= 1}.
  This is a special case of the multivariate asymmetric logistic
  model.

  \code{model = "alog"} (Tawn, 1990)
  
  Let \eqn{B} be the set of all non-empty subsets of
  \eqn{\{1,\ldots,d\}}{{1,\ldots,d}}, let
  \eqn{B_1=\{b \in B:|b|=1\}}{B1={b in B:|b|=1}}, where \eqn{|b|}
  denotes the number of elements in the set \eqn{b}, and let
  \eqn{B_{(i)}=\{b \in B:i \in b\}}{B(i)={b in B:i in b}}.
  The \code{d} dimensional multivariate asymmetric logistic distribution
  function is
  \deqn{G(z)=\exp\left\{-\sum\nolimits_{b \in B} \left[\sum\nolimits_
    {i\in b}(t_{i,b}y_i)^{1/r_b}\right]^{r_b}\right\},}{G(z) =
    exp{-sum{b in B} [sum{i in b}(t{i,b}yi)^(1/r{b})]^r{b}},}
  where the dependence parameters \eqn{r_b\in(0,1]}{r{b} in (0,1]} for
  all \eqn{b\in B \setminus B_1}{b in B\B1}, and the asymmetry parameters
  \eqn{t_{i,b}\in[0,1]}{t{i,b} in [0,1]} for all \eqn{b\in B}{b in B}
  and \eqn{i\in b}{i in b}.
  The constraints
  \eqn{\sum_{b \in B_{(i)}}t_{i,b}=1}{sum{b in B(i)} t{i,b}=1} for
  \eqn{i = 1,\ldots,d}
  ensure that the marginal distributions are generalized extreme value.
  Further constraints arise from the possible redundancy of asymmetry
  parameters in the expansion of the distribution form.
  Let \eqn{b_{-i_0} = \{i \in b:i \neq i_0\}}{
    b_{-i0} = {i in b:i is not i_0}}.
  If \eqn{r_b = 1}{r{b} = 1} for some
  \eqn{b\in B \setminus B_1}{b in B\B1} then
  \eqn{t_{i,b} = 0}{t{i,b} = 0} for all \eqn{i\in b}{i in b}.
  Furthermore, if for some
  \eqn{b\in B \setminus B_1}{b in B\B1},
  \eqn{t_{i,b} = 0}{t{i,b} = 0} for all
  \eqn{i\in b_{-i_0}}{i in b_{-i0}}, then
  \eqn{t_{i_0,b} = 0}{t{i0,b} = 0}.

  \code{dep} should be a vector of length \eqn{2^d-d-1} which contains
  \eqn{\{r_b:b\in B \setminus B_1\}}{{r{b}:b in B\B1}}, with
  the order defined by the natural set ordering on the index.
  For example, for the trivariate model,
  \eqn{\code{dep} = (r_{12},r_{13},r_{23},r_{123})}{
    \code{dep} = (r{12},r{13},r{23},r{123})}.
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
}
\value{
  \code{pmvevd} gives the distribution function, \code{dmvevd} gives
  the density function and \code{rmvevd} generates random deviates, for
  the multivariate logistic or multivariate asymmetric logistic model.
}
\note{
  Multivariate extensions of other bivariate models are more complex.
  A multivariate extension of the Husler-Reiss model exists, involving
  a multidimensional integral and one parameter for each bivariate
  margin.
  Multivariate extensions for the negative logistic model can be derived
  but are considerably more complex and appear to be less flexible.
  The ``multivariate negative logistic model'' often presented in the
  literature (e.g. Kotz \emph{et al}, 2000) is not a valid distribution
  function and should not be used.

  The logistic and asymmetric logistic models respectively are
  simulated using Algorithms 2.1 and 2.2 in Stephenson(2003b).

  The density function of the logistic model is evaluated using
  the representation of Shi(1995).
  The density function of the asymmetric logistic model is evaluated
  using the representation given in Stephenson(2003a).
}
\references{
  Gumbel, E. J. (1960)
  Distributions des valeurs extremes en plusieurs dimensions.
  \emph{Publ. Inst. Statist. Univ. Paris}, \bold{9}, 171--173.
  
  Kotz, S. and Balakrishnan, N. and Johnson, N. L. (2000)
  \emph{Continuous Multivariate Distributions}, vol. 1.
  New York: John Wiley & Sons, 2nd edn.

  Shi, D. (1995)
  Fisher information for a multivariate extreme value distribution.
  \emph{Biometrika}, \bold{82}(3), 644--649.

  Stephenson, A. G. (2003a)
  \emph{Extreme Value Distributions and their Application}.
  Ph.D. Thesis, Lancaster University, Lancaster, UK.
  
  Stephenson, A. G. (2003b)
  Simulating multivariate extreme value distributions of logistic type.
  \emph{Extremes}, \bold{6}(1), 49--60.

  Tawn, J. A. (1990)
  Modelling multivariate extreme value distributions.
  \emph{Biometrika}, \bold{77}, 245--253.
}


\seealso{\code{\link{rbvevd}}, \code{\link{rgev}}}

\examples{
pmvevd(matrix(rep(0:4,5), ncol=5), dep = .7, model = "log", d = 5)
pmvevd(rep(4,5), dep = .7, model = "log", d = 5)
rmvevd(10, dep = .7, model = "log", d = 5)
dmvevd(rep(-1,20), dep = .7, model = "log", d = 20, log = TRUE)

asy <- list(.4, .1, .6, c(.3,.2), c(.1,.1), c(.4,.1), c(.2,.3,.2))
pmvevd(rep(2,3), dep = c(.6,.5,.8,.3), asy = asy, model = "alog", d = 3)
asy <- list(.4, .0, .6, c(.3,.2), c(.1,.1), c(.4,.1), c(.2,.4,.2))
rmvevd(10, dep = c(.6,.5,.8,.3), asy = asy, model = "alog", d = 3)
dmvevd(rep(0,3), dep = c(.6,.5,.8,.3), asy = asy, model = "alog", d = 3)

asy <- list(0, 0, 0, 0, c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0),
  c(.2,.1,.2), c(.1,.1,.2), c(.3,.4,.1), c(.2,.2,.2), c(.4,.6,.2,.5))
rmvevd(10, dep = .7, asy = asy, model = "alog", d = 4)
rmvevd(10, dep = c(rep(1,6), rep(.7,5)), asy = asy, model = "alog", d = 4)
}
\keyword{distribution}

