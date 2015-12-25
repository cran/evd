\name{fextreme}
\alias{fextreme}
\title{Maximum-likelihood Fitting of Maxima and Minima}
\description{
  Maximum-likelihood fitting for the distribution of the
  maximum/minimum of a given number of independent variables from
  a specified distribution.
}
\usage{
fextreme(x, start, densfun, distnfun, \dots, distn, mlen = 1, largest =
    TRUE, std.err = TRUE, corr = FALSE, method = "Nelder-Mead")
}
\arguments{
  \item{x}{A numeric vector.}
  \item{start}{A named list giving the initial values for the
    parameters over which the likelihood is to be maximized.}
  \item{densfun, distnfun}{Density and distribution function of the
    specified distribution.}
  \item{\dots}{Additional parameters, either for the specified
    distribution or for the optimization function \code{optim}.
    If parameters of the distribution are included they will be
    held fixed at the values given (see \bold{Examples}). If
    parameters of the distribution are not included either here
    or as a named component in \code{start} they will be held
    fixed at the default values specified in the corresponding
    density and distribution functions (assuming they exist; an
    error will be generated otherwise).}
  \item{distn}{A character string, optionally specified as an alternative
    to \code{densfun} and \code{distnfun} such that the density and
    distribution functions are formed upon the addition of the
    prefixes \code{d} and \code{p} respectively.}
  \item{mlen}{The number of independent variables.}
  \item{largest}{Logical; if \code{TRUE} (default) use maxima,
    otherwise minima.}
  \item{std.err}{Logical; if \code{TRUE} (the default), the standard
    errors are returned.}
  \item{corr}{Logical; if \code{TRUE}, the correlation matrix is
    returned.}
  \item{method}{The optimization method (see \code{\link{optim}} for
    details).}
}
\details{
  Maximization of the log-likelihood is performed.
  The estimated standard errors are taken from the observed information,
  calculated by a numerical approximation. 

  If the density and distribution functions are user defined, the order
  of the arguments must mimic those in R base (i.e. data first,
  parameters second).
  Density functions must have \code{log} arguments.
}
\value{
  Returns an object of class \code{c("extreme","evd")}.

  The generic accessor functions \code{\link{fitted}} (or
  \code{\link{fitted.values}}), \code{\link{std.errors}},
  \code{\link{deviance}}, \code{\link{logLik}} and
  \code{\link{AIC}} extract various features of the
  returned object.
  The function \code{anova} compares nested models.
  
  An object of class \code{c("extreme","evd")} is a list containing
  at most the following components
  \item{estimate}{A vector containing the maximum likelihood estimates.}
  \item{std.err}{A vector containing the standard errors.}
  \item{deviance}{The deviance at the maximum likelihood estimates.}
  \item{corr}{The correlation matrix.}
  \item{var.cov}{The variance covariance matrix.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}}.}
  \item{call}{The call of the current function.}
  \item{data}{The data passed to the argument \code{x}.}
  \item{n}{The length of \code{x}.}
}

\seealso{\code{\link{anova.evd}}, \code{\link{forder}}, \code{\link{optim}}}

\examples{
uvdata <- rextreme(100, qnorm, mean = 0.56, mlen = 365)
fextreme(uvdata, list(mean = 0, sd = 1), distn = "norm", mlen = 365)
fextreme(uvdata, list(rate = 1), distn = "exp", mlen = 365, 
  method = "Brent", lower=0.01, upper=10)
fextreme(uvdata, list(scale = 1), shape = 1, distn = "gamma", mlen = 365,
  method = "Brent", lower=0.01, upper=10)
fextreme(uvdata, list(shape = 1, scale = 1), distn = "gamma", mlen = 365)
}
\keyword{models}


