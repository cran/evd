\name{fpot}
\alias{fpot}
\alias{print.pot}
\title{Peaks Over Threshold Modelling using the Generalized Pareto
  or Point Process Representation}
\description{
  Maximum-likelihood fitting for peaks over threshold modelling,
  using the Generalized Pareto or Point Process representation,
  allowing any of the parameters to be held fixed if desired.
}
\usage{
fpot(x, threshold, model = c("gpd", "pp"), start, npp = length(x),
    cmax = FALSE, r = 1, ulow = -Inf, rlow = 1, mper = NULL, \dots,
    std.err = TRUE, corr = FALSE, method = "BFGS", warn.inf = TRUE)
}
\arguments{
  \item{x}{A numeric vector, which may contain missing values.}
  \item{threshold}{The threshold.}
  \item{model}{The model; either \code{"gpd"} (the default) or
    \code{"pp"}, for the Generalized Pareto or Point Process
    representations respectively.}
  \item{start}{A named list giving the initial values for the
    parameters over which the likelihood is to be maximized.
    If \code{start} is omitted the routine attempts to find good
    starting values using moment estimators.}
  \item{npp}{The data should contain \code{npp} observations per
    ``period'', where the return level plot produced by
    \code{plot.pot} will represent return periods in units of
    ``periods''. By default \code{npp = length(x)}, so that the
    ``period'' is the period of time over which the entire
    data set is collected. It may often be useful to change
    this default so that more sensible units are used.
    For example, if yearly periodic units are required, use
    \code{npp = 365.25} for daily data and \code{npp = 52.18}
    for weekly data. The argument only makes a difference to
    the actual fit if \code{mper} is not \code{NULL} or if
    \code{model = "pp"} (see \bold{Details}).}
  \item{cmax}{Logical; if \code{FALSE} (the default), the model
    is fitted using all exceedences over the threshold. If
    \code{TRUE}, the model is fitted using cluster maxima, using
    clusters of exceedences derived from \code{clusters}.}
  \item{r, ulow, rlow}{Arguments used for the identification of
    clusters of exceedences (see \code{\link{clusters}}).
    Ignored if \code{cmax} is \code{FALSE} (the default).}
  \item{mper}{Controls the parameterization of the generalized
    Pareto model. Should be either \code{NULL} (the default),
    or a positive number (see \bold{Details}).
    If \code{mper} is not \code{NULL} and \code{model = "pp"},
    an error is returned.}
  \item{\dots}{Additional parameters, either for the model
    or for the optimization function \code{optim}. If parameters
    of the model are included they will be held fixed at the
    values given (see \bold{Examples}).}
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
  The exeedances over the threshold \code{threshold} (if \code{cmax} is
  \code{FALSE}) or the maxima of the clusters of exeedances (if
  \code{cmax} is \code{TRUE}) are (if \code{model = "gpd"}) fitted to a
  generalized Pareto distribution (GPD) with location \code{threshold}.
  If \code{model = "pp"} the exceedances are fitted to a
  non-homogeneous Poisson process (Coles, 2001).
  
  If \code{mper} is \code{NULL} (the default), the parameters of
  the model (if \code{model = "gpd"}) are \code{scale} and
  \code{shape}, for the scale and shape parameters of the GPD.
  If \code{model = "pp"} the parameters are \code{loc}, \code{scale}
  and \code{shape}. Under \code{model = "pp"} the parameters can be
  interpreted as parameters of the Generalized Extreme Value
  distribution, fitted to the maxima of \code{npp} random variables.
  In this case, the value of \code{npp} should be reasonably large.

  For both characterizations, the shape parameters are
  equivalent. The scale parameter under the generalized Pareto
  characterization is equal to \eqn{b + s(u - a)}, where \eqn{a},
  \eqn{b} and \eqn{s} are the location, scale and shape parameters
  under the Point Process characterization, and where \eqn{u} is
  the threshold.

  If \eqn{\code{mper} = m} is a positive value, then
  the generalized Pareto model is reparameterized so that the
  parameters are \code{rlevel} and \code{shape}, where
  \code{rlevel} is the \eqn{m} ``period'' return level, where
  ``period'' is defined via the argument \code{npp}.

  The \eqn{m} ``period'' return level is defined as follows.
  Let \eqn{G} be the fitted generalized Pareto distribution
  function, with location \eqn{\code{threshold} = u}, so that
  \eqn{1 - G(z)} is the fitted probability of an exceedance
  over \eqn{z > u} given an exceedance over \eqn{u}.
  The fitted probability of an exceedance over \eqn{z > u} is
  therefore \eqn{p(1 - G(z))}, where \eqn{p} is the estimated
  probabilty of exceeding \eqn{u}, which is given by the empirical
  proportion of exceedances.
  The \eqn{m} ``period'' return level \eqn{z_m} satisfies
  \eqn{p(1 - G(z_m)) = 1/(mN)}, where \eqn{N} is the number
  of points per period (multiplied by the estimate of the
  extremal index, if cluster maxima are fitted).
  In other words, \eqn{z_m} is the quantile of the fitted model
  that corresponds to the upper tail probability \eqn{1/(mN)}.
  If \code{mper} is infinite, then \eqn{z_m} is the upper end point,
  given by \code{threshold} minus \eqn{\code{scale}/\code{shape}},
  and the shape parameter is then restricted to be negative.
}

\value{
  Returns an object of class \code{c("pot","uvevd","pot")}.

  The generic accessor functions \code{\link{fitted}} (or
  \code{\link{fitted.values}}), \code{\link{std.errors}},
  \code{\link{deviance}}, \code{\link{logLik}} and
  \code{\link{AIC}} extract various features of the
  returned object.

  The function \code{profile} can be
  used to obtain deviance profiles for the model parameters.
  In particular, profiles of the \eqn{m} \code{period}
  return level \eqn{z_m} can be calculated and plotted when
  \eqn{\code{mper} = m}.
  The function \code{anova} compares nested models.
  The function \code{plot} produces diagnostic plots.
  
  An object of class \code{c("pot","uvevd","pot")} is a list containing
  the following components
  \item{estimate}{A vector containing the maximum likelihood estimates.}
  \item{std.err}{A vector containing the standard errors.}
  \item{fixed}{A vector containing the parameters of the model that
    have been held fixed.}
  \item{param}{A vector containing all parameters (optimized and fixed).}
  \item{deviance}{The deviance at the maximum likelihood estimates.}
  \item{corr}{The correlation matrix.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}}.}
  \item{threshold, r, ulow, rlow, npp}{The arguments of the same
    name.}
  \item{nhigh}{The number of exceedences (if \code{cmax} is
    \code{FALSE}) or the number of clusters of exceedences (if
    \code{cmax} is \code{TRUE}).}
  \item{nat, pat}{The number and proportion of exceedences.}
  \item{extind}{The estimate of the extremal index (i.e.
    \code{nhigh} divided by \code{nat}). If \code{cmax} is
    \code{FALSE}, this is \code{NULL}.}
  \item{data}{The data passed to the argument \code{x}.}
  \item{exceedances}{The exceedences, or the maxima of the clusters
    of exceedences.}
  \item{mper}{The argument \code{mper}.}
  \item{scale}{The scale parameter for the fitted generalized Pareto
    distribution. If \code{mper} is \code{NULL} and \code{model = "gpd"}
    (the defaults), this will also be an element of \code{param}.}
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
  \code{\link{plot.uvevd}}, \code{\link{profile.evd}},
  \code{\link{profile2d.evd}}, \code{\link{mrlplot}},
  \code{\link{tcplot}}}

\examples{
uvdata <- rgpd(100, loc = 0, scale = 1.1, shape = 0.2)
M1 <- fpot(uvdata, 1)
M2 <- fpot(uvdata, 1, shape = 0)
anova(M1, M2)
par(mfrow = c(2,2))
plot(M1)
\dontrun{M1P <- profile(M1)}
\dontrun{plot(M1P)}

M1 <- fpot(uvdata, 1, mper = 10)
M2 <- fpot(uvdata, 1, mper = 100)
\dontrun{M1P <- profile(M1, which = "rlevel", conf=0.975, mesh=0.1)}
\dontrun{M2P <- profile(M2, which = "rlevel", conf=0.975, mesh=0.1)}
\dontrun{plot(M1P)}
\dontrun{plot(M2P)}
}
\keyword{models}

