\name{fbvall}
\alias{fbvall}
\title{Simultaneous Maximum-likelihood Fitting of Bivariate Extreme
  Value Distributions}
\description{
  Maximum-likelihood estimates for the bivariate extreme value
  distributions considered in this package, including linear modelling
  of the marginal location parameters. 
}
\usage{
fbvall(x, nsloc1 = NULL, nsloc2 = NULL, omit = NULL, boxcon = TRUE, 
    std.err = TRUE, orderby = c("AIC", "BIC", "SC"),
    control = list(maxit = 250)) 
}
\arguments{
  \item{x}{A matrix or data frame with two columns, which may contain
    missing values.}
  \item{nsloc1, nsloc2}{A data frame with the same number of rows as the
    length of \code{x}, to be passed to the individual fitting
    functions, for linear modelling of the location parameter on the
    first/second margin.
    The data frames are treated as covariate matrices (excluding the
    intercept).
    A numeric vector can be given as an alternative to a single column
    data frame.}
  \item{omit}{A character vector giving the models which are not to
    be fit, containing any of ``log'', ``alog'', ``hr'', ``neglog'',
    ``aneglog'', ``bilog'', ``negbilog'' or ``ct'' for the logistic,
    asymmetric logistic, Husler-Reiss, negative logistic, asymmetric
    negative logistic, bilogistic, negative bilogistic and Coles-Tawn
    models respectively.}
  \item{boxcon}{Logical; if \code{TRUE} (the default), the
    \code{L-BFGS-B} optimization method is used, which includes box
    constraints. If \code{FALSE} the \code{BFGS} method is used.
    The \code{BFGS} method is much faster, and should be preferred
    if the estimates are known to lie in the interior of the parameter
    space.}
  \item{std.err}{Logical; if \code{TRUE} (the default), the ``standard
    errors'' are returned. A logical vector can also be given. This
    should be the same length as the number of models being fitted.}
  \item{orderby}{The (column) order of the models in the returned list
    components. The models are ordered by best fit, as measured by
    Akaike's information criterion (the default), Bayesian information
    criterion or Schwarz's criterion.}
  \item{control}{To be passed to \code{\link{optim}}. Only options
    which are independent of the number of parameters within the
    optimization should be given. Some options are related to the
    optimization method used. This is defined by the \code{boxcon}
    argument (see above).}
}
\details{
  Maximization of the log-likelihood for every bivariate model is performed.
  The \code{std.err} component of the returned list is taken from the
  observed information, calculated by a numerical approximation.
  The ``standard errors'' must be interpreted with caution because
  the usual asymptotic properties of maximum likelihood estimators
  may not hold (Smith, 1985).

  Any bivariate extreme value distribution can be written as
  \deqn{G(z_1,z_2) = \exp\left[-(y_1+y_2)A\left(
    \frac{y_1}{y_1+y_2}\right)\right]}{
    G(z1,z2) = exp{-(y1+y2)A[y1/(y1+y2)]}}
  for some function \eqn{A(\cdot)}{A()} defined on \eqn{[0,1]}, where    
  \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,2}, with the marginal
  parameters given by \eqn{(a_i,b_i,s_i)}{(ai,bi,si)},
  \eqn{b_i > 0}{bi > 0}.

  \eqn{A(\cdot)}{A()} is called (by some authors) the dependence
  function.
  It follows that \eqn{A(0)=A(1)=1}, and that \eqn{A(\cdot)}{A()} is
  a convex function with \eqn{\max(x,1-x) \leq A(x)\leq 1}{max(x,1-x)
    <= A(x) <= 1} for all \eqn{0\leq x\leq1}{0 <= x <= 1}.
  \eqn{A(\cdot)}{A()} does not depend on the marginal parameters.

  The component \code{dep.summary} of the returned list contains
  three values summarizing the dependence structure for each fitted
  model, based on \eqn{A(\cdot)}{A()}.
  There are two measures of dependence.
  The first is given by \eqn{2(1-A(1/2))}.
  The second is the integral of \eqn{4(1 - A(x))}, taken over
  \eqn{0\leq x\leq1}{0 <= x <= 1}.
  These appear in the rows of \code{dep.summary} labelled by \code{dep}
  and \code{intdep} respectively.
  Both measures are zero at independence and one at complete dependence.
  
  The final row of \code{dep.summary}, labelled \code{intasy}, contains
  a measure of asymmetry given by the integral of
  \eqn{4(A(x) - A(1-x))/(3 - 2\sqrt2)}{
  4(A(x) - A(1-x))/(3 - 2 sqrt(2))},
  taken over \eqn{0 \leq x \leq 0.5}{0 <= x <= 0.5}.
  This integral lies in the closed interval [-1,1] (conjecture), with
  larger absolute values representing stronger asymmetry.
  As a rough guide, any value within the interval \eqn{(-0.25,0.25)}
  suggests that the dependence structure is close to symmetric.
  For the symmetric models \eqn{A(x) = A(1-x)} for all
  \eqn{0 \leq x \leq 0.5}{0 <= x <= 0.5}, so the integral will be zero.
}
\note{
  The measure \eqn{2(1-A(1/2))} is standard theory but the two integral
  measures of dependence and asymmetry are not.
  I have proposed a second measure of dependence because a comparison
  of \eqn{2(1-A(1/2))} across symmetric and strongly asymmetric models
  may not be very enlightening.
  I am open to suggestions of other simple one-number summaries of
  the dependence structure (where simple excludes bivariate integrals).
}
\value{
  A list with with components
  \item{estimate}{A matrix containing the maximum likelihood estimates.}
  \item{std.err}{A matrix containing the ``standard errors'' (if
    \eqn{\code{std.err} = \code{TRUE}}).}
  \item{deviance}{The deviance at the maximum likelihood estimates
    for each fitted model.}
  \item{criteria}{A matrix containing Akaike's information criterion,
    Bayesian information criterion and Schwarz's criterion for each
    fitted model.}
  \item{dep.summary}{A matrix containing three values summarizing the
    dependence structure for each fitted model (see \bold{Details}).}
}
\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}

\seealso{\code{\link{fbvlog}}, \code{\link{optim}}}

\examples{
\dontrun{data(sealevel)}
\dontrun{fbvall(sealevel)}
\dontrun{fbvall(sealevel, nsloc1 = (-40:40)/100)}
}

\keyword{distribution}

