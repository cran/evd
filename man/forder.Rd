\name{forder}
\alias{forder}
\title{Maximum-likelihood Fitting of Order Statistics}
\description{
  Maximum-likelihood fitting for the distribution of a selected order
  statistic of a given number of independent variables from a
  specified distribution.
}
\usage{
forder(x, start, densfun, distnfun, \dots, distn, mlen = 1, j = 1, 
largest = TRUE, std.err = TRUE, corr = FALSE, method = "Nelder-Mead")
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
    distribution and functions are formed upon the addition of the
    prefixes \code{d} and \code{p} respectively.}
  \item{mlen}{The number of independent variables.}
  \item{j}{The order statistic, taken as the \code{j}th largest
    (default) or smallest of \code{mlen}, according to the value of
    \code{largest}.}
  \item{largest}{Logical; if \code{TRUE} (default) use the \code{j}th
    largest order statistic, otherwise use the \code{j}th smallest.}
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
  This class is defined in \code{\link{fextreme}}.

  The generic accessor functions \code{\link{fitted}} (or
  \code{\link{fitted.values}}), \code{\link{std.errors}},
  \code{\link{deviance}}, \code{\link{logLik}} and
  \code{\link{AIC}} extract various features of the
  returned object.
  The function \code{anova} compares nested models.
}

\seealso{\code{\link{anova.evd}}, \code{\link{fextreme}},
  \code{\link{optim}}}

\examples{
uvd <- rorder(100, qnorm, mean = 0.56, mlen = 365, j = 2)
forder(uvd, list(mean = 0, sd = 1), distn = "norm", mlen = 365, j = 2)
forder(uvd, list(rate = 1), distn = "exp", mlen = 365, j = 2)
forder(uvd, list(scale = 1), shape = 1, distn = "gamma", mlen = 365, j = 2)
forder(uvd, list(shape = 1, scale = 1), distn = "gamma", mlen = 365, j = 2)
}
\keyword{models}


