\name{uvmle}
\alias{ffrechet}
\alias{fgumbel}
\alias{frweibull}
\alias{fgev}
\title{Maximum-likelihood Fitting of Univariate Extreme Value Distributions}
\description{
  Maximum-likelihood estimates for the Frechet, Gumbel, reversed Weibull
  and generalized extreme value distributions, including linear
  modelling of the location parameter, and allowing any of the
  parameters to be held fixed if desired.
}
\usage{
ffrechet(x, start, \dots, nsloc = NULL, std.err = TRUE, method = "BFGS")
fgumbel(x, start, \dots, nsloc = NULL, std.err = TRUE, method = "BFGS") 
frweibull(x, start, \dots, nsloc = NULL, std.err = TRUE, method = "BFGS") 
fgev(x, start, \dots, nsloc = NULL, std.err = TRUE, method = "BFGS")
}
\arguments{
  \item{x}{A numeric vector, which may contain missing values.}
  \item{start}{A named list giving the parameters to be optimized with
    initial values (see \bold{Details}).
    If \code{start} is omitted the routine attempts to find good
    starting values using moment estimators.}
  \item{\dots}{Additional parameters, either for the distribution or for
    \code{optim}. If parameters of the distribution are included they
    will be held fixed (see \bold{Details}).}
  \item{nsloc}{A data frame with the same number of rows as the
    length of \code{x}, for linear modelling of the location
    parameter (see \bold{Details}).
    The data frame is treated as a covariate matrix (excluding the
    intercept).
    A numeric vector can be given as an alternative to a single column
    data frame.}
  \item{std.err}{Logical; if \code{TRUE} (the default), the ``standard
    errors'' are returned.}
  \item{method}{The optimization method (see \code{\link{optim}} for
    details).}
}
\details{
  Maximization of the log-likelihood is performed.
  The \code{std.err} component of the returned list is taken from the
  observed information, calculated by a numerical approximation.
  The ``standard errors'' must be interpreted with caution because
  the usual asymptotic properties of maximum likelihood estimators
  may not hold (Smith, 1985).
  
  For stationary fitting, the parameters can be passed (either as
  named components of \code{start} or as fixed values) using
  \code{loc}, \code{scale} and \code{shape}.

  For non-stationary fitting, the non-stationary parameters can be passed
  using the column names of the data frame \code{nsloc} with the prefix
  ``loc''. The intercept of the linear model can be passed as
  \code{loc}. If \code{nsloc1} is a vector it is converted into
  a one column data frame with column name ``trend'', so
  that the associated parameter can be passed as \code{loc1trend}.

  For non-stationary fitting it is recommended that the covariates
  within the linear model for the location parameter are (at least
  approximately) centered and scaled, particularly if automatic starting
  values are used, since the starting values for all the associated
  parameters are taken to be zero. 
}

\value{
  A list with with components
  \item{estimate}{a vector containing the maximum likelihood estimates.}
  \item{std.err}{a vector containing the ``standard errors'' (if
    \eqn{\code{std.err} = \code{TRUE}}).}
  \item{deviance}{the deviance at the maximum likelihood estimates.}
  \item{counts}{a two-element integer vector giving the number of calls
    to the log-likelihood function and the number gradient approximations
    respectively.}
}
\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}

\seealso{\code{\link{optim}}}

\examples{
data <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
fgev(data, control = list(trace = 1))
fgev(data, shape = 0)
fgev(data, scale = 1, shape = 0)
fgev(data, start = list(loc = 0), scale = 1, shape = 0)
fgev(data, nsloc = (-49:50)/100)
fgev(data, nsloc = data.frame(trend = (-49:50)/100, random =
runif(100, min=-.5, max=.5)))
fgev(data, nsloc = data.frame(trend = (-49:50)/100, random =
runif(100, min=-.5, max=.5)), locrandom = 0)
}
\keyword{distribution}

