\name{fgev}
\alias{fgev}
\alias{fitted.evd}
\alias{std.errors}
\alias{std.errors.evd}
\alias{print.evd}
\alias{logLik.evd}

\title{Maximum-likelihood Fitting of the Generalized Extreme Value
  Distribution}
\description{
  Maximum-likelihood fitting for the generalized extreme value distribution,
  including linear modelling of the location parameter, and allowing any
  of the parameters to be held fixed if desired.
}
\usage{
fgev(x, start, \dots, nsloc = NULL, prob = NULL, std.err = TRUE,
    corr = FALSE, method = "BFGS", warn.inf = TRUE)
}
\arguments{
  \item{x}{A numeric vector, which may contain missing values.}
  \item{start}{A named list giving the initial values for the
    parameters over which the likelihood is to be maximized.
    If \code{start} is omitted the routine attempts to find good
    starting values using moment estimators.}
  \item{\dots}{Additional parameters, either for the GEV model
    or for the optimization function \code{optim}. If parameters
    of the model are included they will be held fixed at the
    values given (see \bold{Examples}).}
  \item{nsloc}{A data frame with the same number of rows as the
    length of \code{x}, for linear modelling of the location
    parameter.
    The data frame is treated as a covariate matrix (excluding the
    intercept).
    A numeric vector can be given as an alternative to a single column
    data frame.}
  \item{prob}{Controls the parameterization of the model (see
    \bold{Details}). Should be either \code{NULL} (the default),
    or a probability in the closed interval [0,1].}
  \item{std.err}{Logical; if \code{TRUE} (the default), the standard
    errors are returned.}
  \item{corr}{Logical; if \code{TRUE}, the correlation matrix is
    returned.}
  \item{method}{The optimization method (see \code{\link{optim}} for
    details).}
  \item{warn.inf}{Logical; if \code{TRUE} (the default), a warning is
    given if the negative log-likelihood is infinite when evaluated at
    the starting values.}
}
\details{
  If \code{prob} is \code{NULL} (the default):
  
  For stationary models the parameter names are \code{loc}, \code{scale}
  and \code{shape}, for the location, scale and shape parameters
  respectively.
  For non-stationary models, the parameter names are \code{loc},
  \code{loc}\emph{x1}, \dots, \code{loc}\emph{xn}, \code{scale} and
  \code{shape}, where \emph{x1}, \dots, \emph{xn} are the column names
  of \code{nsloc}, so that \code{loc} is the intercept of the
  linear model, and \code{loc}\emph{x1}, \dots, \code{loc}\emph{xn}
  are the \code{ncol(nsloc)} coefficients.
  If \code{nsloc} is a vector it is converted into a single column
  data frame with column name \code{trend}, and hence the associated
  trend parameter is named \code{loctrend}.

  If \eqn{\code{prob} = p} is a probability:
  
  The fit is performed using a different parameterization.
  Let \eqn{a}, \eqn{b} and \eqn{s} denote the location, scale
  and shape parameters of the GEV distribution.
  For stationary models, the distribution is parameterized
  using \eqn{(z_p,b,s)}, where
  \deqn{z_p = a - b/s (1 - (-\log(1 - p))^s)}{
        z_p = a - b/s (1 - (-log(1 - p))^s)}
  is such that \eqn{G(z_p) = 1 - p}, where \eqn{G} is the
  GEV distribution function.
  \eqn{\code{prob} = p} is therefore the probability in the upper
  tail corresponding to the quantile \eqn{z_p}.
  If \code{prob} is zero, then \eqn{z_p} is the upper end point
  \eqn{a - b/s}, and \eqn{s} is restricted to the negative
  (Weibull) axis.
  If \code{prob} is one, then \eqn{z_p} is the lower end point
  \eqn{a - b/s}, and \eqn{s} is restricted to the positive
  (Frechet) axis.
  The parameter names are \code{quantile}, \code{scale}
  and \code{shape}, for \eqn{z_p}, \eqn{b} and \eqn{s}
  respectively.
  
  For non-stationary models the parameter \eqn{z_p} is again given by
  the equation above, but \eqn{a} becomes the intercept of the linear
  model for the location parameter, so that \code{quantile} replaces
  (the intercept) \code{loc}, and hence the parameter names are
  \code{quantile}, \code{loc}\emph{x1}, \dots, \code{loc}\emph{xn},
  \code{scale} and \code{shape}, where \emph{x1}, \dots, \emph{xn} are
  the column names of \code{nsloc}.

  In either case:
  
  For non-stationary fitting it is recommended that the covariates
  within the linear model for the location parameter are (at least
  approximately) centered and scaled (i.e.\ that the columns of
  \code{nsloc} are centered and scaled), particularly if automatic
  starting values are used, since the starting values for the
  associated parameters are then zero.
}

\value{
  Returns an object of class \code{c("gev","evd")}.

  The generic accessor functions \code{\link{fitted}} (or
  \code{\link{fitted.values}}), \code{\link{std.errors}},
  \code{\link{deviance}}, \code{\link{logLik}} and
  \code{\link{AIC}} extract various features of the
  returned object.

  The functions \code{profile} and \code{profile2d} are
  used to obtain deviance profiles for the model parameters.
  In particular, profiles of the quantile \eqn{z_p} can be
  calculated and plotted when \eqn{\code{prob} = p}.
  The function \code{anova} compares nested models.
  The function \code{plot} produces diagnostic plots.
  
  An object of class \code{c("gev","evd")} is a list containing
  at most the following components
  \item{estimate}{A vector containing the maximum likelihood estimates.}
  \item{std.err}{A vector containing the standard errors.}
  \item{fixed}{A vector containing the parameters of the model that
    have been held fixed.}
  \item{param}{A vector containing all parameters (optimized and fixed).}
  \item{deviance}{The deviance at the maximum likelihood estimates.}
  \item{corr}{The correlation matrix.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}}.}
  \item{data}{The data passed to the argument \code{x}.}
  \item{tdata}{The data, transformed to stationarity (for non-stationary
    models).}
  \item{nsloc}{The argument \code{nsloc}.}
  \item{n}{The length of \code{x}.}
  \item{prob}{The argument \code{prob}.}
  \item{loc}{The location parameter. If \code{prob} is \code{NULL}
    (the default), this will also be an element of \code{param}.}
  \item{call}{The call of the current function.}
}
\section{Warning}{
  The standard errors and the correlation matrix in the returned
  object are taken from the observed information, calculated by a
  numerical approximation.
  They must be interpreted with caution when the shape parameter
  is less than \eqn{-0.5}, because the usual asymptotic
  properties of maximum likelihood estimators do not then
  hold (Smith, 1985).
}
\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}

\seealso{\code{\link{anova.evd}}, \code{\link{optim}},
  \code{\link{plot.gev}}, \code{\link{profile.evd}},
  \code{\link{profile2d.evd}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
trend <- (-49:50)/100
M1 <- fgev(uvdata, nsloc = trend, control = list(trace = 1))
M2 <- fgev(uvdata)
M3 <- fgev(uvdata, shape = 0)
M4 <- fgev(uvdata, scale = 1, shape = 0)
anova(M1, M2, M3, M4)
\dontrun{plot(M2)}
\dontrun{M2P <- profile(M2)}
\dontrun{plot(M2P)}

rnd <- runif(100, min = -.5, max = .5)
fgev(uvdata, nsloc = data.frame(trend = trend, random = rnd))
fgev(uvdata, nsloc = data.frame(trend = trend, random = rnd), locrandom = 0)

uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata, prob = 0.1)
M2 <- fgev(uvdata, prob = 0.01)
\dontrun{M1P <- profile(M1, which = "quantile")}
\dontrun{M2P <- profile(M2, which = "quantile")}
\dontrun{plot(M1P)}
\dontrun{plot(M2P)}
}
\keyword{models}

