\name{fbvpot}
\alias{fbvpot}
\alias{print.bvpot}

\title{Maximum-likelihood Fitting of Bivariate Extreme Value
  Distributions to Threshold Exceedances}
\description{
  Fit models for one of nine parametric bivariate extreme-value
  distributions using threshold exceedances, allowing any of the
  parameters to be held fixed if desired.
}
\synopsis{
fbvpot(x, threshold, model = c("log", "bilog", "alog", "neglog",
    "negbilog", "aneglog", "ct", "hr", "amix"), likelihood =
    c("censored", "poisson"), start, \dots, sym = FALSE, cshape =
    cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method =
    "BFGS", warn.inf = TRUE)
}
\usage{
fbvpot(x, threshold, model = "log", start, \dots, sym = FALSE, cshape =
   cscale, cscale = FALSE, std.err = TRUE, corr = FALSE, method = "BFGS")
}
\arguments{
  \item{x}{A matrix or data frame with two columns. If this
    contains missing values, those values are treated as if
    they fell below the corresponding marginal threshold.}
  \item{threshold}{A vector of two thresholds.}
  \item{model}{The specified model; a character string. Must be
    either \code{"log"} (the default), \code{"alog"}, \code{"hr"},
    \code{"neglog"}, \code{"aneglog"}, \code{"bilog"},
    \code{"negbilog"}, \code{"ct"} or \code{"amix"} (or any unique
    partial match), for the logistic, asymmetric logistic,
    Husler-Reiss, negative logistic, asymmetric negative logistic,
    bilogistic, negative bilogistic, Coles-Tawn and asymmetric mixed
    models respectively. The definition of each model is given in
    \code{\link{rbvevd}}.}
  \item{start}{A named list giving the initial values for all of the
    parameters in the model. If \code{start} is omitted the routine
    attempts to find good starting values using marginal maximum
    likelihood estimators.}
  \item{\dots}{Additional parameters, either for the bivariate extreme
    value model or for the optimization function \code{optim}. If
    parameters of the model are included they will be held fixed at
    the values given (see \bold{Examples}).}
  \item{sym}{Logical; if \code{TRUE}, the dependence structure
    of the models \code{"alog"}, \code{"aneglog"} or \code{"ct"} are
    constrained to be symmetric (see \bold{Details}). For all other
    models, the argument is ignored (and a warning is given).}
  \item{cshape}{Logical; if \code{TRUE}, a common shape parameter is
    fitted to each margin.}
  \item{cscale}{Logical; if \code{TRUE}, a common scale parameter is
    fitted to each margin, and the default value of \code{cshape}
    is then \code{TRUE}, so that under this default common marginal
    parameters are fitted.}
  \item{std.err}{Logical; if \code{TRUE} (the default), the standard
    errors are returned.}
  \item{corr}{Logical; if \code{TRUE}, the correlation matrix is
    returned.}
  \item{method}{The optimization method (see \code{\link{optim}} for
    details).}
}
\details{
  The bivariate peaks over threshold models are fitted by
  maximizing the censored likelihood as given in e.g. Section
  8.3.1 of Coles(2001). Declustering of the series is not
  currently implemented.
  
  The dependence parameter names are one or more of \code{dep},
  \code{asy1}, \code{asy2}, \code{alpha} and \code{beta}, depending on
  the model selected (see \code{\link{rbvevd}}). 
  The marginal parameter names are \code{scale1} and \code{shape1}
  for the first margin, and \code{scale2} and \code{shape2} for the
  second margin.

  If \code{cshape} is true, the models are constrained so that
  \code{shape2 = shape1}. The parameter \code{shape2} is then
  taken to be specified, so that e.g. the common shape
  parameter can only be fixed at zero using \code{shape1 = 0},
  since using \code{shape2 = 0} gives an error. Similar
  comments apply for \code{cscale}.
  
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
  logistic or asymmetric mixed models (see \bold{fbvevd}).

  For numerical reasons the parameters of each model are subject the
  artificial constraints given in \code{\link{fbvevd}}.  
}
\value{
  Returns an object of class \code{c("bvpot","evd")}.

  The generic accessor functions \code{\link{fitted}} (or
  \code{\link{fitted.values}}), \code{\link{std.errors}},
  \code{\link{deviance}}, \code{\link{logLik}} and
  \code{\link{AIC}} extract various features of the
  returned object.

  The functions \code{profile} and \code{profile2d} can be
  used to obtain deviance profiles.
  The function \code{anova} compares nested models, and the
  function \code{AIC} compares non-nested models.
  There is currently no plot method available.

  An object of class \code{c("bvpot","evd")} is a list containing
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
  \item{dep.summary}{A value summarizing the strength of
    dependence in the fitted model (see \bold{fbvevd}).}
  \item{corr}{The correlation matrix.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}}.}
  \item{data}{The data passed to the argument \code{x}.}
  \item{threshold}{The argument \code{threshold}.}
  \item{n}{The number of rows in \code{x}.}
  \item{nat}{The vector of length three containing the number of
    exceedances on the first, second and both margins respectively.}
  \item{sym}{The argument \code{sym}.}
  \item{cmar}{The vector \code{c(cscale, cshape)}.}
  \item{model}{The argument \code{model}.}
  \item{call}{The call of the current function.}
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
  Coles, S. G. (2001)
  \emph{An Introduction to Statistical Modelling of Extreme
  Values}, London: Springer--Verlag.

  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}
\author{Chris Ferro and Alec Stephenson}

\seealso{\code{\link{abvevd}}, \code{\link{anova.evd}},
  \code{\link{fbvevd}}, \code{\link{optim}}, \code{\link{rbvevd}}}

\examples{
bvdata <- rbvevd(1000, dep = 0.5, model = "log")
u <- apply(bvdata, 2,  quantile, probs = 0.9)
M1 <- fbvpot(bvdata, u, model = "log")
M2 <- fbvpot(bvdata, u, "log", dep = 0.5)
anova(M1, M2)
}

\keyword{models}

