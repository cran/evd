\name{bvevd}
\alias{dbvevd}
\alias{pbvevd}
\alias{rbvevd}
\title{Parametric Bivariate Extreme Value Distributions}
\description{
  Density function, distribution function and random generation
  for nine parametric bivariate extreme value models. 
}
\synopsis{
dbvevd(x, dep, asy = c(1, 1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
    mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvevd(q, dep, asy = c(1, 1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
    mar1 = c(0, 1, 0), mar2 = mar1, lower.tail = TRUE) 
rbvevd(n, dep, asy = c(1, 1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
    mar1 = c(0, 1, 0), mar2 = mar1) 
}
\usage{
dbvevd(x, dep, asy = c(1, 1), alpha, beta, model = "log",
    mar1 = c(0, 1, 0), mar2 = mar1, log = FALSE) 
pbvevd(q, dep, asy = c(1, 1), alpha, beta, model = "log",
    mar1 = c(0, 1, 0), mar2 = mar1, lower.tail = TRUE) 
rbvevd(n, dep, asy = c(1, 1), alpha, beta, model = "log",
    mar1 = c(0, 1, 0), mar2 = mar1) 
}
\arguments{
  \item{x, q}{A vector of length two or a matrix with two columns,
    in which case the density/distribution is evaluated across
    the rows.}
  \item{n}{Number of observations.}
  \item{dep}{Dependence parameter for the logistic, asymmetric
    logistic, Husler-Reiss, negative logistic and asymmetric
    negative logistic models.}
  \item{asy}{A vector of length two, containing the two asymmetry
    parameters for the asymmetric logistic and asymmetric negative
    logistic models.}
  \item{alpha, beta}{Alpha and beta parameters for the bilogistic,
    negative bilogistic, Coles-Tawn and asymmetric mixed models.}
  \item{model}{The specified model; a character string. Must be
    either \code{"log"} (the default), \code{"alog"}, \code{"hr"},
    \code{"neglog"}, \code{"aneglog"}, \code{"bilog"},
    \code{"negbilog"}, \code{"ct"} or \code{"amix"} (or any unique
    partial match), for the logistic, asymmetric logistic,
    Husler-Reiss, negative logistic, asymmetric negative logistic,
    bilogistic, negative bilogistic, Coles-Tawn and asymmetric
    mixed models respectively. If parameter arguments are given
    that do not correspond to the specified model those arguments
    are ignored, with a warning.}
  \item{mar1, mar2}{Vectors of length three containing marginal
    parameters, or matrices with three columns where each
    column represents a vector of values to be passed to the
    corresponding marginal parameter.}
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
  \eqn{i = 1,2}, where the marginal parameters are given by
  \eqn{\code{mari} = (a_i,b_i,s_i)}{\code{mari} = (ai,bi,si)},
  \eqn{b_i > 0}{bi > 0}.
  If \eqn{s_i = 0}{si = 0} then \eqn{y_i}{yi} is defined by
  continuity.
  
  In each of the bivariate distributions functions
  \eqn{G(z_1,z_2)}{G(z1,z2)} given below, the univariate margins
  are generalized extreme value, so that
  \eqn{G(z_i) = \exp(-y_i)}{G(zi) = \exp(-yi)} for \eqn{i = 1,2}.
  If \eqn{1+s_i(z_i-a_i)/b_i \leq 0}{1+si(zi-ai)/bi <= 0} for some
  \eqn{i = 1,2}, the value \eqn{z_i}{zi} is either greater than the
  upper end point (if \eqn{s_i < 0}{si < 0}), or less than the lower
  end point (if \eqn{s_i > 0}{si > 0}), of the \eqn{i}th univariate
  marginal distribution.
  
  \code{model = "log"} (Gumbel, 1960)
  
  The bivariate logistic distribution function with
  parameter \eqn{\code{dep} = r} is
  \deqn{G(z_1,z_2) = \exp\left[-(y_1^{1/r}+y_2^{1/r})^r\right]}{
    G(z1,z2) = exp{-[y1^(1/r)+y2^(1/r)]^r}}
  where \eqn{0 < r \leq 1}{0 < r <= 1}.
  This is a special case of the bivariate asymmetric logistic
  model.
  Complete dependence is obtained in the limit as
  \eqn{r} approaches zero.
  Independence is obtained when \eqn{r = 1}.

  
  \code{model = "alog"} (Tawn, 1988)
  
  The bivariate asymmetric logistic distribution function with
  parameters \eqn{\code{dep} = r} and
  \eqn{\code{asy} = (t_1,t_2)}{\code{asy} = (t1,t2)} is
  \deqn{G(z_1,z_2) = \exp\left\{-(1-t_1)y_1-(1-t_2)y_2-
    [(t_1y_1)^{1/r}+(t_2y_2)^{1/r}]^r\right\}}{
    G(z1,z2) = exp{-(1-t1)y1-(1-t2)y2-[(t1y1)^(1/r)+(t2y2)^(1/r)]^r}}
  where \eqn{0 < r \leq 1}{0 < r <= 1} and
  \eqn{0 \leq t_1,t_2 \leq 1}{0 <= t1,t2 <= 1}.
  When \eqn{t_1 = t_2 = 1}{t1 = t2 = 1} the asymmetric logistic
  model is equivalent to the logistic model.
  Independence is obtained when either \eqn{r = 1},
  \eqn{t_1 = 0}{t1 = 0} or \eqn{t_2 = 0}{t2 = 0}.
  Complete dependence is obtained in the limit when
  \eqn{t_1 = t_2 = 1}{t1 = t2 = 1} and \eqn{r}
  approaches zero.
  Different limits occur when \eqn{t_1}{t1} and \eqn{t_2}{t2}
  are fixed and \eqn{r} approaches zero.


  \code{model = "hr"} (Husler and Reiss, 1989)
  
  The Husler-Reiss distribution function with parameter
  \eqn{\code{dep} = r} is
  \deqn{G(z_1,z_2) = \exp\left(-y_1\Phi\{r^{-1}+{\textstyle\frac{1}{2}}
    r[\log(y_1/y_2)]\} - y_2\Phi\{r^{-1}+{\textstyle\frac{1}{2}}r
    [\log(y_2/y_1)]\}\right)}{
    G(z1,z2) = exp(-y1 Phi{r^{-1}+r[log(y1/y2)]/2} -
    y2 Phi{r^{-1}+r[log(y2/y1)]/2}}
  where \eqn{\Phi(\cdot)}{Phi()} is the standard normal distribution
  function and \eqn{r > 0}.
  Independence is obtained in the limit as \eqn{r} approaches zero.
  Complete dependence is obtained as \eqn{r} tends to infinity.


  \code{model = "neglog"} (Galambos, 1975)

  The bivariate negative logistic distribution function
  with parameter \eqn{\code{dep} = r} is
  \deqn{G(z_1,z_2) = \exp\left\{-y_1-y_2+
    [y_1^{-r}+y_2^{-r}]^{-1/r}\right\}}{
    G(z1,z2) = exp{-y1-y2+[y1^(-r)+y2^(-r)]^(-1/r)}}
  where \eqn{r > 0}.
  This is a special case of the bivariate asymmetric negative
  logistic model.
  Independence is obtained in the limit as \eqn{r} approaches zero.
  Complete dependence is obtained as \eqn{r} tends to infinity.
  The earliest reference to this model appears to be
  Galambos (1975, Section 4).


  \code{model = "aneglog"} (Joe, 1990)
  
  The bivariate asymmetric negative logistic distribution function
  with parameters parameters \eqn{\code{dep} = r} and
  \eqn{\code{asy} = (t_1,t_2)}{\code{asy} = (t1,t2)} is
  \deqn{G(z_1,z_2) = \exp\left\{-y_1-y_2+
    [(t_1y_1)^{-r}+(t_2y_2)^{-r}]^{-1/r}\right\}}{
    G(z1,z2) = exp{-y1-y2+[(t1y1)^(-r)+(t2y2)^(-r)]^(-1/r)}}
  where \eqn{r > 0} and \eqn{0 < t_1,t_2 \leq 1}{0 < t1,t2 <= 1}.
  When \eqn{t_1 = t_2 = 1}{t1 = t2 = 1} the asymmetric negative
  logistic model is equivalent to the negative logistic model.
  Independence is obtained in the limit as either \eqn{r},
  \eqn{t_1}{t1} or \eqn{t_2}{t2} approaches zero.
  Complete dependence is obtained in the limit when
  \eqn{t_1 = t_2 = 1}{t1 = t2 = 1} and \eqn{r}
  tends to infinity.
  Different limits occur when \eqn{t_1}{t1} and \eqn{t_2}{t2}
  are fixed and \eqn{r} tends to infinity.
  The earliest reference to this model appears to be Joe (1990),
  who introduces a multivariate extreme value distribution which
  reduces to \eqn{G(z_1,z_2)}{G(z1,z2)} in the bivariate case.


  \code{model = "bilog"} (Smith, 1990)
  
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
  \eqn{0 < \alpha,\beta < 1}{0 < alpha,beta < 1}.
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


  \code{model = "negbilog"} (Coles and Tawn, 1994)

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
  \eqn{\alpha > 0}{alpha > 0} and \eqn{\beta > 0}{beta > 0}.
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


  \code{model = "ct"} (Coles and Tawn, 1991)
  
  The Coles-Tawn distribution function with
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
  \eqn{\code{shape1} = \alpha}{\code{shape1} = alpha} and
  \eqn{\code{shape2} = \beta}{\code{shape2} = beta}.
  Complete dependence is obtained in the limit as
  \eqn{\alpha = \beta}{alpha = beta} tends to infinity.
  Independence is obtained as
  \eqn{\alpha = \beta}{alpha = beta} approaches zero, and when
  one of \eqn{\alpha,\beta}{alpha,beta} is fixed and the other
  approaches zero.
  Different limits occur when one of
  \eqn{\alpha,\beta}{alpha,beta} is fixed and the other
  tends to infinity.

  \code{model = "amix"} (Tawn, 1988)

  The asymmetric mixed distribution function with
  parameters \eqn{\code{alpha} = \alpha}{\code{alpha} = alpha}
  and \eqn{\code{beta} = \beta}{\code{beta} = beta} has
  a dependence function with the following cubic polynomial
  form.
  \deqn{A(t) = 1 - (\alpha +\beta)t + \alpha t^2 + \beta t^3}{
    A(t) = 1 - (\alpha +\beta)t + \alpha t^2 + \beta t^3}
  where \eqn{\alpha}{alpha} and \eqn{\alpha + 3\beta}{alpha + 3beta}
  are non-negative, and where \eqn{\alpha + \beta}{alpha + beta}
  and \eqn{\alpha + 2\beta}{alpha + 2beta} are less than or equal
  to one.
  These constraints imply that beta lies in the interval [-0.5,0.5]
  and that alpha lies in the interval [0,1.5], though alpha can
  only be greater than one if beta is negative. The strength
  of dependence increases for increasing alpha (for fixed beta).
  Complete dependence cannot be obtained.
  Independence is obtained when both parameters are zero.
  For the definition of a dependence function, see
  \code{\link{abvevd}}.
}
\value{
  \code{dbvevd} gives the density function, \code{pbvevd} gives the
  distribution function and \code{rbvevd} generates random deviates,
  for one of nine parametric bivariate extreme value models.
}
\note{
  The logistic and asymmetric logistic models respectively are
  simulated using bivariate versions of Algorithms 1.1 and 1.2 in
  Stephenson(2003).
  All other models are simulated using a root finding algorithm
  to simulate from the conditional distributions.

  The simulation of the bilogistic and negative bilogistic models
  requires a root finding algorithm to evaluate \eqn{q}
  within the root finding algorithm used to simulate from the
  conditional distributions.
  The generation of bilogistic and negative bilogistic random
  deviates is therefore relatively slow (about 2.8 seconds per
  1000 random vectors on a 450MHz PIII, 512Mb RAM).

  The bilogistic and negative bilogistic models can be represented
  under a single model, using the integral of the maximum of two
  beta distributions (Joe, 1997).

  The Coles-Tawn model is called the Dirichelet model in Coles
  and Tawn (1991).
}
\references{
  Coles, S. G. and Tawn, J. A. (1991)
  Modelling extreme multivariate events.
  \emph{J. Roy. Statist. Soc., B}, \bold{53}, 377--392.
  
  Coles, S. G. and Tawn, J. A. (1994)
  Statistical methods for multivariate extremes: an application to
  structural design (with discussion).
  \emph{Appl. Statist.}, \bold{43}, 1--48.
  
  Galambos, J. (1975)
  Order statistics of samples from multivariate distributions.
  \emph{J. Amer. Statist. Assoc.}, \bold{70}, 674--680.
  
  Gumbel, E. J. (1960)
  Distributions des valeurs extremes en plusieurs dimensions.
  \emph{Publ. Inst. Statist. Univ. Paris}, \bold{9}, 171--173.

  Husler, J. and Reiss, R.-D. (1989)
  Maxima of normal random vectors: between independence 
  and complete dependence.
  \emph{Statist. Probab. Letters}, \bold{7}, 283--286.

  Joe, H. (1990)
  Families of min-stable multivariate exponential and multivariate
  extreme value distributions.
  \emph{Statist. Probab. Letters}, \bold{9}, 75--81.

  Joe, H. (1997)
  \emph{Multivariate Models and Dependence Concepts},
  London: Chapman & Hall.

  Smith, R. L. (1990)
  Extreme value theory. In
  \emph{Handbook of Applicable Mathematics} (ed. W. Ledermann),
  vol. 7. Chichester: John Wiley, pp. 437--471.
  
  Stephenson, A. G. (2003)
  Simulating multivariate extreme value distributions of logistic type.
  \emph{Extremes}, \bold{6}(1), 49--60.

  Tawn, J. A. (1988)
  Bivariate extreme value theory: models and estimation.
  \emph{Biometrika}, \bold{75}, 397--415.  
}

\seealso{\code{\link{abvevd}}, \code{\link{rgev}}, \code{\link{rmvevd}}}


\examples{
pbvevd(matrix(rep(0:4,2), ncol=2), dep = 0.7, model = "log")
pbvevd(c(2,2), dep = 0.7, asy = c(0.6,0.8), model = "alog")
pbvevd(c(1,1), dep = 1.7, model = "hr")

margins <- cbind(0, 1, seq(-0.5,0.5,0.1))
rbvevd(11, dep = 1.7, model = "hr", mar1 = margins)
rbvevd(10, dep = 1.2, model = "neglog", mar1 = c(10, 1, 1))
rbvevd(10, alpha = 0.7, beta = 0.52, model = "bilog")

dbvevd(c(0,0), dep = 1.2, asy = c(0.5,0.9), model = "aneglog")
dbvevd(c(0,0), alpha = 0.75, beta = 0.5, model = "ct", log = TRUE)
dbvevd(c(0,0), alpha = 0.7, beta = 1.52, model = "negbilog")
}
\keyword{distribution}




