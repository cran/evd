\name{fbvevd}
\alias{fbvevd}
\alias{print.bvevd}

\title{Maximum-likelihood Fitting of Bivariate Extreme Value Distributions}
\description{
  Fit models for one of nine parametric bivariate extreme value
  distributions, including linear modelling of the marginal
  location parameters, and allowing any of the parameters to be held
  fixed if desired. 
}
\usage{
fbvevd(x, model = c("log", "alog", "hr", "neglog", "aneglog", "bilog",
    "negbilog", "ct", "amix"), start, \dots, sym = FALSE,
    nsloc1 = NULL, nsloc2 = NULL, cshape = cscale, cscale = cloc,
    cloc = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS",
    warn.inf = TRUE)
}
\arguments{
  \item{x}{A matrix or data frame, ordinarily with two columns,
    which may contain missing values. A data frame may also
    contain a third column of mode \code{logical}, which
    itself may contain missing values (see \bold{More Details}).}
  \item{model}{The specified model; a character string. Must be
    either \code{"log"} (the default), \code{"alog"}, \code{"hr"},
    \code{"neglog"}, \code{"aneglog"}, \code{"bilog"},
    \code{"negbilog"}, \code{"ct"} or \code{"amix"} (or any unique
    partial match), for the logistic, asymmetric logistic,
    Husler-Reiss, negative logistic, asymmetric negative logistic,
    bilogistic, negative bilogistic, Coles-Tawn and asymmetric
    mixed models respectively. The definition of each model is
    given in \code{\link{rbvevd}}.}
  \item{start}{A named list giving the initial values for the
    parameters over which the likelihood is to be maximized.
    If \code{start} is omitted the routine attempts to find good
    starting values using marginal maximum likelihood estimators.}
  \item{\dots}{Additional parameters, either for the bivariate extreme
    value model or for the optimization function \code{optim}. If
    parameters of the model are included they will be held fixed at
    the values given (see \bold{Examples}).}
  \item{sym}{Logical; if \code{TRUE}, the dependence structure
    of the models \code{"alog"}, \code{"aneglog"} or \code{"ct"} are
    constrained to be symmetric (see \bold{Details}). For all other
    models, the argument is ignored (and a warning is given).}
  \item{nsloc1, nsloc2}{A data frame with the same number of rows as
    \code{x}, for linear modelling of the location parameter on the
    first/second margin (see \bold{Details}).
    The data frames are treated as covariate matrices, excluding the
    intercept. A numeric vector can be given as an alternative to a
    single column data frame.}
  \item{cshape}{Logical; if \code{TRUE}, a common shape parameter is
    fitted to each margin.}
  \item{cscale}{Logical; if \code{TRUE}, a common scale parameter is
    fitted to each margin, and the default value of \code{cshape}
    is then \code{TRUE}, so that under this default common scale
    and shape parameters are fitted.}
  \item{cloc}{Logical; if \code{TRUE}, a common location parameter is
    fitted to each margin, and the default values of \code{cshape}
    and \code{cscale} are then \code{TRUE}, so that under these
    defaults common marginal parameters are fitted.}
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
  The dependence parameter names are one or more of \code{dep},
  \code{asy1}, \code{asy2}, \code{alpha} and \code{beta}, depending on
  the model selected (see \code{\link{rbvevd}}).  The marginal parameter
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
  scaled (i.e. that the columns of \code{nsloc1} and \code{nsloc2}
  are centered and scaled), particularly if automatic starting values
  are used, since the starting values for the associated parameters are
  then zero. If \code{cloc} is \code{TRUE}, both \code{nsloc1} and
  \code{nsloc2} must be identical, since a common linear model is
  then implemented on both margins.

  If \code{cshape} is true, the models are constrained so that
  \code{shape2 = shape1}. The parameter \code{shape2} is then
  taken to be specified, so that e.g. the common shape
  parameter can only be fixed at zero using \code{shape1 = 0},
  since using \code{shape2 = 0} gives an error. Similar
  comments apply for \code{cscale} and \code{cloc}.
  
  If \code{sym} is \code{TRUE}, the asymmetric logistic and
  asymmetric negative logistic models are constrained so that
  \code{asy2 = asy1}, and the Coles-Tawn model is constrained
  so that \code{beta = alpha}. The parameter \code{asy2} or
  \code{beta} is then taken to be specified, so that e.g.
  the parameters \code{asy1} and \code{asy2} can only
  be fixed at \code{0.8} using \code{asy1 = 0.8}, since
  using \code{asy2 = 0.8} gives an error.
  
  Bilogistic and negative bilogistic models constrained to
  symmetry are logistic and negative logistic models
  respectively. The (symmetric) mixed model (e.g. Tawn, 1998)
  can be obtained as a special case of the asymmetric
  logistic or asymmetric mixed models (see \bold{Examples}).

  The value \code{Dependence} given in the printed output
  is \eqn{2(1-A(1/2))}, where \eqn{A} is the estimated dependence
  function (see \code{\link{abvevd}}). It measures the strength of
  dependence, and lies in the interval [0,1]; at independence and
  complete dependence it is zero and one respectively (Coles,
  Heffernan and Tawn, 1999). See \code{\link{chiplot}} for
  further information.
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
  \item{fixed}{A vector containing the parameters that have been fixed
    at specific values within the optimization.}
  \item{fixed2}{A vector containing the parameters that have been set
    to be equal to other model parameters.}
  \item{param}{A vector containing all parameters (those optimized,
    those fixed to specific values, and those set to be
    equal to other model parameters).}
  \item{deviance}{The deviance at the maximum likelihood estimates.}
  \item{dep.summary}{The estimate of \eqn{2(1-A(1/2))}.}
  \item{corr}{The correlation matrix.}
  \item{var.cov}{The variance covariance matrix.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}}.}
  \item{data}{The data passed to the argument \code{x}.}
  \item{tdata}{The data, transformed to stationarity (for non-stationary
    models).}
  \item{nsloc1, nsloc2}{The arguments \code{nsloc1} and \code{nsloc2}.}
  \item{n}{The number of rows in \code{x}.}
  \item{sym}{The argument \code{sym}.}
  \item{cmar}{The vector \code{c(cloc, cscale, cshape)}.}
  \item{model}{The argument \code{model}.}
  \item{call}{The call of the current function.}
}
\section{More Details}{
  If \code{x} is a data frame with a third column of mode
  \code{logical}, then the model is fitted using the likelihood
  derived by Stephenson and Tawn (2004). This is appropriate
  when each bivariate data point comprises componentwise maxima
  from some underlying bivariate process, and where the
  corresponding logical value denotes whether or not the maxima
  were caused by the same event within that process.
  
  Under this scheme the diagnostic plots that are produced
  using \code{plot} are somewhat different to those described
  in \code{\link{plot.bvevd}}: the density, dependence function
  and quantile curves plots contain fitted functions for
  observations where the logical case is unknown, and the
  conditional P-P plots condition on both the logical case and
  the given margin (which requires numerical integration at each
  data point).
}
\section{Artificial Constraints}{
  For numerical reasons parameters are subject to artificial
  constraints. Specifically, these constraints are: marginal
  scale parameters not less than 0.01; \code{dep} not less
  than [0.1] [0.2] [0.05] in [logistic] [Husler-Reiss]
  [negative logistic] models; \code{dep} not greater
  than [10] [5] in [Husler-Reiss] [negative logistic] models;
  \code{asy1} and \code{asy2} not less than 0.001;
  \code{alpha} and \code{beta} not less than [0.1] [0.1]
  [0.001] in [bilogistic] [negative bilogistic] [Coles-Tawn]
  models; \code{alpha} and \code{beta} not greater than [0.999]
  [20] [30] in [bilogistic] [negative bilogistic] [Coles-Tawn]
  models.
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
  Coles, S. G., Heffernan, J. and Tawn, J. A. (1999)
  Dependence measures for extreme value analyses.
  \emph{Extremes}, \bold{2}, 339--365.
  
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.

  Stephenson, A. G. and Tawn, J. A. (2004)
  Exploiting Occurence Times in Likelihood Inference for
  Componentwise Maxima.
  \emph{Biometrika} \bold{92}(1), 213--217.

  Tawn, J. A. (1988)
  Bivariate extreme value theory: models and estimation.
  \emph{Biometrika}, \bold{75}, 397--415.
}

\seealso{\code{\link{anova.evd}}, \code{\link{optim}},
  \code{\link{plot.bvevd}}, \code{\link{profile.evd}},
  \code{\link{profile2d.evd}}, \code{\link{rbvevd}}}

\examples{
bvdata <- rbvevd(100, dep = 0.6, model = "log", mar1 = c(1.2,1.4,0.4))
M1 <- fbvevd(bvdata, model = "log")
M2 <- fbvevd(bvdata, model = "log", dep = 0.75)
anova(M1, M2)
par(mfrow = c(2,2))
plot(M1)
plot(M1, mar = 1)
plot(M1, mar = 2)
\dontrun{par(mfrow = c(1,1))}
\dontrun{M1P <- profile(M1, which = "dep")}
\dontrun{plot(M1P)}

trend <- (-49:50)/100
rnd <- runif(100, min = -.5, max = .5)
fbvevd(bvdata, model = "log", nsloc1 = trend)
fbvevd(bvdata, model = "log", nsloc1 = trend, nsloc2 = data.frame(trend
= trend,  random = rnd))
fbvevd(bvdata, model = "log", nsloc1 = trend, nsloc2 = data.frame(trend
= trend, random = rnd), loc2random = 0)

bvdata <- rbvevd(100, dep = 1, asy = c(0.5,0.5), model = "anegl")
anlog <- fbvevd(bvdata, model = "anegl")
mixed <- fbvevd(bvdata, model = "anegl", dep = 1, sym = TRUE)
anova(anlog, mixed)
amixed <- fbvevd(bvdata, model = "amix")
mixed <- fbvevd(bvdata, model = "amix", beta = 0)
anova(amixed, mixed)
}

\keyword{models}

