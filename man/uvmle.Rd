\name{uvmle}
\alias{ffrechet}
\alias{fgumbel}
\alias{frweibull}
\alias{fgev}
\title{Maximum-likelihood Fitting of Univariate Extreme Value Distributions}
\description{
  Maximum-likelihood estimates for the Frechet, Gumbel, reversed Weibull
  and generalized extreme value distributions, allowing parameters to be
  held fixed if desired. 
}
\usage{
ffrechet(x, start, \dots)
fgumbel(x, start, \dots) 
frweibull(x, start, \dots) 
fgev(x, start, \dots)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{start}{a named list giving the parameters to be optimized with
    initial values.}
  \item{\dots}{additional parameters, either for the distribution or for
    \code{optim}. If parameters of the distribution are included they
    will be held fixed. If parameters of the distribution are not included
    either here or as a named component in \code{start} they will be held
    fixed at the default values specified in the corresponding density
    function.}
}
\details{
  Maximization of the log-likelihood is performed.
  The estimated standard errors are taken from the observed information,
  calculated by a numerical approximation. 
}
\note{
  These functions are simple wrappers to \code{fdensfun} which is based
  on \code{fitdistr} in package \code{MASS}.
}
\value{
  A list with with components
  \item{estimate}{a vector containing the maximum likelihood estimates.}
  \item{std.err}{a vector containing the standard errors.}
  \item{deviance}{the deviance at the maximum likelihood estimates.}
  \item{counts}{a two-element integer vector giving the number of calls
    to the log-likelihood function and the number gradient approximations
    respectively.}
}

\seealso{\code{\link{optim}}}

\examples{
data <- rgev(1000, loc = 0.13, scale = 1.1, shape = 0.2)
fgev(data, start = list(loc=0,scale=1,shape=0), method = "BFGS", 
control = list(trace = 1))
fgev(data, start = list(loc=0,scale=1), shape = 0, method = "BFGS")
fgev(data, start = list(loc=0), scale = 1, shape = 0)
}
\keyword{distribution}

