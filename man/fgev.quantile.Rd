\name{fgev.quantile}
\alias{fgev.quantile}
\title{Maximum-likelihood Fitting of the Generalized Extreme Value Distribution}
\description{
  Maximum-likelihood fitting for the generalized extreme value distribution,
  including linear modelling of the location parameter, and allowing any
  of the parameters to be held fixed if desired.
}
\usage{
fgev.quantile(x, start, \dots, prob, nsloc = NULL, std.err = TRUE,
    corr = FALSE, method = "Nelder-Mead", warn.inf = TRUE)
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
  \item{prob}{The probability in the \emph{upper} tail that corresponds
    to the \code{quantile} parameter (see \bold{Details}).}
  \item{nsloc}{A data frame with the same number of rows as the
    length of \code{x}, for linear modelling of the location
    parameter (see \bold{Details}).
    The data frame is treated as a covariate matrix (excluding the
    intercept).
    A numeric vector can be given as an alternative to a single column
    data frame.}
  \item{std.err}{Logical; if \code{TRUE} (the default), the ``standard
    errors'' are returned.}
  \item{corr}{Logical; if \code{TRUE}, the ``correlation matrix'' is
    returned.}
  \item{method}{The optimization method (see \code{\link{optim}} for
    details).}
  \item{warn.inf}{Logical; if \code{TRUE} (the default), a warning is
    given if the negative log-likelihood is infinite when evaluated at
    the starting values.}
}
\details{
  Maximization of the log-likelihood is performed.
  The ``standard errors'' and the ``correlation matrix'' in the returned
  object are taken from the observed information, calculated by a
  numerical approximation.
  They must be interpreted with caution because the usual asymptotic
  properties of maximum likelihood estimators may not hold (Smith,
  1985).

  \code{fgev.quant} performs the same fit as \code{\link{fgev}},
  but under a different parameterization.
  For stationary fitting using \code{fgev}, the likelihood is optimized
  over \eqn{(a,b,s)}, where \eqn{a} is the location parameter, \eqn{b} is
  the scale and \eqn{s} is the shape.
  For stationary fitting using \code{fgev.quant}, the likelihood is
  optimized over \eqn{(z_p,b,s)}, where \eqn{\code{prob} = p} and
  \deqn{z_p = a - b/s (1 - (-\log(1 - p))^s)}{
        z_p = a - b/s (1 - (-log(1 - p))^s)}
  is such that \eqn{G(z_p) = 1 - p}, where \eqn{G} is the
  generalized extreme value distribution function.
  For non-stationary fitting the parameter \eqn{z_p} is again given by the
  equation above but \eqn{a} becomes the intercept of the linear model
  for the location parameter.
  
  For stationary fitting, the parameters can be passed (either as
  named components of \code{start} or as fixed values) using
  \code{quantile}, \code{scale} and \code{shape}.

  For non-stationary fitting, the non-stationary parameters can be passed
  using the column names of the data frame \code{nsloc} with the prefix
  ``loc''. If \code{nsloc1} is a vector it is converted into
  a one column data frame with column name ``trend'', so
  that the associated parameter can be passed as \code{loctrend}.

  For non-stationary fitting it is recommended that the covariates
  within the linear model for the location parameter are (at least
  approximately) centered and scaled, particularly if automatic starting
  values are used, since the starting values for all the associated
  parameters are taken to be zero. 
}

\value{
  Returns an object of class \code{"evd"}.

  The generic accessor functions \code{\link{fitted}} (or
  \code{\link{fitted.values}}), \code{\link{std.errors}} and
  \code{\link{deviance}} extract various features of the returned
  object.

  The functions \code{profile} and \code{profile2d} are
  used to obtain deviance profiles.
  In particular, profiles of the quantile \eqn{z_p} can be calculated
  and plotted.
  The function \code{anova} compares nested models.
  The function \code{plot} produces diagnostic plots.
  
  An object of class \code{"evd"} is a list containing at most
  the following components
  \item{estimate}{A vector containing the maximum likelihood estimates.}
  \item{std.err}{A vector containing the ``standard errors''.}
  \item{fixed}{A vector containing the parameters that have been set
    to fixed values within the optimization.}
  \item{param}{A vector containing all parameters (optimized and fixed).}
  \item{deviance}{The deviance at the maximum likelihood estimates.}
  \item{corr}{The ``correlation matrix''.}
  \item{convergence,counts,message}{Components taken from the
    list returned by \code{\link{optim}}.}
  \item{call}{The call of the current function.}
  \item{data}{The data passed to the argument \code{x}.}
  \item{tdata}{The data, transformed to stationarity (for non-stationary
    models).}
  \item{nsloc}{The argument \code{nsloc}.}
  \item{n}{The length of \code{x}.}
  \item{model}{A character string describing the fitted model.}
}

\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}

\seealso{\code{\link{anova.evd}}, \code{\link{optim}},
  \code{\link{plot.evd}}}, \code{\link{profile.evd}}},
  \code{\link{profile2d.evd}}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev.quantile(uvdata, prob = 0.1)
M2 <- fgev.quantile(uvdata, prob = 0.01)
\dontrun{M1P <- profile(M1, which = "quantile")}
\dontrun{M2P <- profile(M2, which = "quantile")}
\dontrun{plot(M1P)}
\dontrun{plot(M2P)}
}
\keyword{models}

