\name{ordermle}
\alias{forder}
\title{Maximum-likelihood Fitting of Order Statistics}
\description{
  Maximum-likelihood estimates for the distribution of a given order
  statistic of a sample from any specified distribution function.
}
\usage{
forder(x, start, densfun, distnfun, \dots, distn, mlen = 1, j = 1, 
largest = TRUE, std.err = TRUE, method = "Nelder-Mead")
}
\arguments{
  \item{x}{A numeric vector.}
  \item{start}{A named list giving the parameters to be optimized with
    initial values.}
  \item{densfun, distnfun}{Density and distribution function of the
    specified distribution.}
  \item{\dots}{Additional parameters, either for the specified
    distribution or for \code{optim}. If parameters of the distribution
    are included they will be held fixed. If parameters of the
    distribution are not included either here or as a named component
    in \code{start} they will be held fixed at the default values
    specified in the corresponding density and distribution functions
    (assuming they exist; an error will be generated otherwise).}
  \item{distn}{A character string, optionally specified as an alternative
    to \code{densfun} and \code{distnfun} such that the density and
    distribution and functions are formed upon the addition of the
    prefixes \code{d} and \code{p} respectively.}
  \item{mlen}{The size of the random sample.}
  \item{j}{The order statistic, taken as the \code{j}th largest
    (default) or smallest of \code{mlen}, according to the value of
    \code{largest}.}
  \item{largest}{Logical; if \code{TRUE} (default) use the \code{j}th
    largest order statistic, otherwise use the \code{j}th smallest.}
  \item{std.err}{Logical; if \code{TRUE} (the default), the standard
    errors are returned.}
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

  Although \code{j} can be any integer in \eqn{[1,\code{mlen}]} it
  is preferable, for computational reasons, to specify
  \code{j} in the interval \eqn{[1,\code{ceiling}(\code{mlen}/2)]},
  which can always be achieved using the argument \code{largest}.
}
\value{
  A list with with components
  \item{estimate}{a vector containing the maximum likelihood estimates.}
  \item{std.err}{a vector containing the standard errors (if
    \eqn{\code{std.err} = \code{TRUE}}).}
  \item{deviance}{the deviance at the maximum likelihood estimates.}
  \item{counts}{a two-element integer vector giving the number of calls
    to the log-likelihood function and the number gradient approximations
    respectively.}
}

\seealso{\code{\link{fext}}, \code{\link{optim}}}

\examples{
data <- rorder(100, qnorm, mean = 0.56, mlen = 365, j = 2)
forder(data, list(mean = 0, sd = 1), distn = "norm", mlen = 365, j = 2)
forder(data, list(rate = 1), distn = "exp", mlen = 365, j = 2)
forder(data, list(scale = 1), shape = 1, distn = "gamma", mlen = 365,
j = 2)
forder(data, list(shape = 1, scale = 1), distn = "gamma", mlen = 365,
j = 2)
}
\keyword{distribution}


