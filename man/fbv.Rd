\name{fbv}
\alias{fbvlog}
\alias{fbvalog}
\alias{fbvhr}
\alias{fbvneglog}
\alias{fbvaneglog}
\alias{fbvbilog}
\alias{fbvnegbilog}
\alias{fbvct}
\title{Maximum-likelihood Fitting of Bivariate Extreme Value Distributions}
\description{
  Fit models for bivariate extreme value distributions;
  symmetric/asymmetric logistic, symmetric/asymmetric negative
  logistic, bilogistic, negative bilogistic, Husler-Reiss and
  Coles-Tawn models, including linear modelling of the marginal
  location parameters, and allowing any of the parameters to be held
  fixed if desired. 
}
\usage{
fbvlog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
    corr = FALSE, method = "BFGS", warn.inf = TRUE)
fbvalog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
    corr = FALSE, method = "BFGS", warn.inf = TRUE)
fbvhr(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
    corr = FALSE, method = "BFGS", warn.inf = TRUE)
fbvneglog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
    corr = FALSE, method = "BFGS", warn.inf = TRUE)
fbvaneglog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
    corr = FALSE, method = "BFGS", warn.inf = TRUE)
fbvbilog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE, 
    corr = FALSE, method = "BFGS", warn.inf = TRUE) 
fbvnegbilog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE, 
    corr = FALSE, method = "BFGS", warn.inf = TRUE)
fbvct(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE, 
    corr = FALSE, method = "BFGS", warn.inf = TRUE) 
}
\arguments{
  \item{x}{A matrix or data frame with two columns, which may contain
    missing values.}
  \item{start}{A named list giving the parameters to be optimized with
    starting values (see \bold{Details}).
    If \code{start} is omitted the routine attempts to find good
    starting values using marginal maximum likelihood estimators.}
  \item{\dots}{Additional parameters, either for the distribution or for
    \code{optim}. If parameters of the distribution are included they
    will be held fixed (see \bold{Details}).}
  \item{nsloc1, nsloc2}{A data frame with the same number of rows as
    \code{x}, for linear modelling of the location
    parameter on the first/second margin (see \bold{Details}).
    The data frames are treated as covariate matrices, excluding the
    intercept.
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
  properties of maximum likelihood estimators may not hold (Smith, 1985).

  For stationary fitting, the named components of \code{start} can be
  any of \code{dep}, \code{asy} (a vector of length two), \code{alpha},
  \code{beta}, \code{mar1} (a vector of length three) and \code{mar2}
  (a vector of length three).
  Alternatively, \code{mar1} can be passed as three separate arguments;
  \code{loc1}, \code{scale1} and \code{shape1}.
  \code{mar2} can be similarly passed as three separate arguments;
  \code{loc2}, \code{scale2} and \code{shape2}.
  The asymmetry parameters for the asymmetric logistic and negative
  asymmetric logistic models can also be passed separately as
  \code{asy1} and \code{asy2}.
  Associating a separate name with each parameter allows any parameter
  subset to be fixed at specified values.
  All parameters to be fixed must be passed individually (e.g. to fix
  the parameters on the first margin all of \code{loc1}, \code{scale1}
  and \code{shape1} must be specified; using \code{mar1} to specify all
  three simultaneously results in an error).

  For non-stationary fitting, the non-stationary parameters on the
  first margin can be passed (either in \code{start} or as fixed values)
  using the column names of the data frame \code{nsloc1} with the prefix
  ``loc1''. The intercept of the linear model can be passed as
  \code{loc1}. If \code{nsloc1} is a vector it is converted into
  a one column data frame with column name ``trend'', so
  that the associated parameter can be passed as \code{loc1trend}.
  Parameters on the second margin can be passed similarly.

  For non-stationary fitting it is recommended that the covariates
  within the linear models for the location parameters are (at least
  approximately) centered and scaled, particularly if automatic starting
  values are used, since the starting values for all the associated
  parameters are taken to be zero.  
}
\value{
  Returns an object of class \code{c("bvevd","evd")}.

  The generic accessor functions \code{\link{fitted}} (or
  \code{\link{fitted.values}}), \code{\link{std.errors}} and
  \code{\link{deviance}} extract various features of the returned
  object.

  The functions \code{profile} and \code{profile2d} can be
  used to obtain deviance profiles.
  The function \code{anova} compares nested models.
  The function \code{plot} produces diagnostic plots.

  An object of class \code{c("bvevd","evd")} is a list containing
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
  \item{nsloc1,nsloc2}{The arguments \code{nsloc1} and \code{nsloc2}.}
  \item{n}{The number of rows in \code{x}.}
  \item{model}{A character string describing the fitted model.}
}
\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}

\seealso{\code{\link{anova.evd}}, \code{\link{optim}},
  \code{\link{plot.bvevd}}, \code{\link{profile.evd}},
  \code{\link{profile2d.evd}}}

\examples{
bvdata <- rbvlog(100, dep = 0.6, mar1 = c(1.2,1.4,0.4),
mar2 = c(1.2,1.4,0.4))
M1 <- fbvlog(bvdata)
M2 <- fbvlog(bvdata, dep = 0.75)
anova(M1, M2)
\dontrun{plot(M1)}
\dontrun{plot(M1, mar = 1)}
\dontrun{plot(M1, mar = 2)}
\dontrun{plot(M2)}
\dontrun{M1P <- profile(M1, which = "dep")}
\dontrun{plot(M1P)}
trend <- (-49:50)/100
rnd <- runif(100, min = -.5, max = .5)
fbvlog(bvdata, nsloc1 = trend)
fbvlog(bvdata, nsloc1 = trend, nsloc2 = data.frame(trend = trend,
random = rnd))
fbvlog(bvdata, nsloc1 = trend, nsloc2 = data.frame(trend = trend,
random = rnd), loc2random = 0)
}

\keyword{models}

