\name{fbvevd}
\alias{fbvevd}
\alias{print.bvevd}

\title{Maximum-likelihood Fitting of Bivariate Extreme Value Distributions}
\description{
  Fit models for one of eight parametric bivariate extreme value
  distributions, including linear modelling of the marginal
  location parameters, and allowing any of the parameters to be held
  fixed if desired. 
}
\synopsis{
fbvevd(x, model = c("log", "alog", "hr", "neglog", "aneglog", "bilog",
    "negbilog", "ct"), start, \dots, nsloc1 = NULL, nsloc2 = NULL,
    std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS",
    warn.inf = TRUE)
}
\usage{
fbvevd(x, model = "log", start, \dots, nsloc1 = NULL, nsloc2 = NULL,
    std.err = TRUE, dsm = TRUE, corr = FALSE, method = "BFGS",
    warn.inf = TRUE)
}
\arguments{
  \item{x}{A matrix or data frame with two columns, which may contain
    missing values.}
  \item{model}{The specified model; a character string. Must be
    either \code{"log"} (the default), \code{"alog"}, \code{"hr"},
    \code{"neglog"}, \code{"aneglog"}, \code{"bilog"},
    \code{"negbilog"} or \code{"ct"} (or any unique partial match),
    for the logistic, asymmetric logistic, Husler-Reiss, negative
    logistic, asymmetric negative logistic, bilogistic, negative
    bilogistic and Coles-Tawn models respectively. The definition
    of each model is given in \code{\link{rbvevd}}.}
  \item{start}{A named list giving the initial values for the
    parameters over which the likelihood is to be maximized.
    If \code{start} is omitted the routine attempts to find good
    starting values using marginal maximum likelihood estimators.}
  \item{\dots}{Additional parameters, either for the bivariate extreme
    value model or for the optimization function \code{optim}. If
    parameters of the model are included they will be held fixed at
    the values given (see \bold{Examples}).}
  \item{nsloc1, nsloc2}{A data frame with the same number of rows as
    \code{x}, for linear modelling of the location parameter on the
    first/second margin (see \bold{Details}).
    The data frames are treated as covariate matrices, excluding the
    intercept. A numeric vector can be given as an alternative to a
    single column data frame.}
  \item{std.err}{Logical; if \code{TRUE} (the default), the standard
    errors are returned.}
  \item{dsm}{Logical; if \code{TRUE} (the default), summaries of the
    dependence structure are returned.}
  \item{corr}{Logical; if \code{TRUE}, the correlation matrix is
    returned.}
  \item{method}{The optimization method (see \code{\link{optim}} for
    details).}
  \item{warn.inf}{Logical; if \code{TRUE} (the default), a warning is
    given if the negative log-likelihood is infinite when evaluated at
    the starting values.}
}
\details{
  The dependence parameter names are one or more of \code{dep},
  \code{asy1}, \code{asy2}, \code{alpha} and \code{beta}, depending on
  the model selected (see \code{\link{rbvevd}}). The marginal parameter
  names are \code{loc1}, \code{scale1} and \code{shape1} for the first
  margin, and \code{loc2}, \code{scale2} and \code{shape2} for the
  second margin.
  If \code{nsloc1} is not \code{NULL}, so that a linear model is
  implemented for the first marginal location parameter, the parameter
  names for the first margin are \code{loc1}, \code{loc1}\emph{x1},
  \dots, \code{loc1}\emph{xn}, \code{scale} and \code{shape}, where
  \emph{x1}, \dots, \emph{xn} are the column names of \code{nsloc1},
  so that \code{loc1} is the intercept of the linear model, and
  \code{loc1}\emph{x1}, \dots, \code{loc1}\emph{xn} are the
  \code{ncol(nsloc1)} coefficients.
  When \code{nsloc2} is not \code{NULL}, the parameter names for the
  second margin are constructed similarly.
  
  It is recommended that the covariates within the linear models for
  the location parameters are (at least approximately) centered and
  scaled (i.e.\ that the columns of \code{nsloc1} and \code{nsloc2}
  are centered and scaled), particularly if automatic starting values
  are used, since the starting values for the associated parameters are
  then zero.

  If \code{dsm} is \code{TRUE}, three values are returned which
  summarize the dependence structure, based on the fitted
  dependence function \eqn{A} (see \code{\link{abvpar}}).
  Two are measures of the strength of dependence.
  The first (Dependence One) is given by \eqn{2(1-A(1/2))}.
  The second (Dependence Two) is the integral of \eqn{4(1 - A(x))},
  taken over \eqn{0\leq x\leq1}{0 <= x <= 1}.
  Both measures are zero at independence and one at complete dependence.
  
  The third value (Asymmetry) is a measure of asymmetry, given by
  the integral of
  \eqn{4(A(x) - A(1-x))/(3 - 2\sqrt2)}{
  4(A(x) - A(1-x))/(3 - 2 sqrt(2))},
  taken over \eqn{0 \leq x \leq 0.5}{0 <= x <= 0.5}.
  This lies in the closed interval [-1,1] (conjecture), with
  larger absolute values representing stronger asymmetry.
  For the logistic, Husler-Reiss and negative logistic models
  \eqn{A(x) = A(1-x)} for all \eqn{0 \leq x \leq 0.5}{0 <= x <= 0.5},
  so the value will be zero.

  For numerical reasons the parameters of each model are subject the
  artificial constraints given in Table 1 of the User's Guide.  
}
\value{
  Returns an object of class \code{c("bvevd","evd")}.

  The generic accessor functions \code{\link{fitted}} (or
  \code{\link{fitted.values}}), \code{\link{std.errors}},
  \code{\link{deviance}}, \code{\link{logLik}} and
  \code{\link{AIC}} extract various features of the
  returned object.

  The functions \code{profile} and \code{profile2d} can be
  used to obtain deviance profiles.
  The function \code{anova} compares nested models, and the
  function \code{AIC} compares non-nested models.
  The function \code{plot} produces diagnostic plots.

  An object of class \code{c("bvevd","evd")} is a list containing
  the following components
  \item{estimate}{A vector containing the maximum likelihood estimates.}
  \item{std.err}{A vector containing the standard errors.}
  \item{fixed}{A vector containing the parameters that have been set
    to fixed values within the optimization.}
  \item{param}{A vector containing all parameters (optimized and fixed).}
  \item{deviance}{The deviance at the maximum likelihood estimates.}
  \item{dep.summary}{A vector of three values, summarizing the
    dependence structure of the fitted model (see \bold{Details}).}
  \item{corr}{The correlation matrix.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}}.}
  \item{call}{The call of the current function.}
  \item{data}{The data passed to the argument \code{x}.}
  \item{tdata}{The data, transformed to stationarity (for non-stationary
    models).}
  \item{nsloc1, nsloc2}{The arguments \code{nsloc1} and \code{nsloc2}.}
  \item{n}{The number of rows in \code{x}.}
  \item{model}{The argument \code{model}.}
}
\section{Warning}{
  The standard errors and the correlation matrix in the returned
  object are taken from the observed information, calculated by a
  numerical approximation.
  They must be interpreted with caution when either of the
  marginal shape parameters are less than \eqn{-0.5}, because
  the usual asymptotic properties of maximum likelihood estimators
  do not then hold (Smith, 1985).
}
\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}

\seealso{\code{\link{anova.evd}}, \code{\link{optim}},
  \code{\link{plot.bvevd}}, \code{\link{profile.evd}},
  \code{\link{profile2d.evd}}, \code{\link{rbvevd}}}

\examples{
bvdata <- rbvevd(100, dep = 0.6, model = "log", mar1 = c(1.2,1.4,0.4))
M1 <- fbvevd(bvdata, model = "log")
M2 <- fbvevd(bvdata, model = "log", dep = 0.75)
anova(M1, M2)
\dontrun{plot(M1)}
\dontrun{plot(M1, mar = 1)}
\dontrun{plot(M1, mar = 2)}
\dontrun{plot(M2)}
\dontrun{M1P <- profile(M1, which = "dep")}
\dontrun{plot(M1P)}

trend <- (-49:50)/100
rnd <- runif(100, min = -.5, max = .5)
fbvevd(bvdata, model = "log", nsloc1 = trend)
fbvevd(bvdata, model = "log", nsloc1 = trend, nsloc2 = data.frame(trend
= trend,  random = rnd))
fbvevd(bvdata, model = "log", nsloc1 = trend, nsloc2 = data.frame(trend
= trend, random = rnd), loc2random = 0)
}

\keyword{models}
