\name{bvmle}
\alias{fbvlog}
\alias{fbvalog}
\alias{fbvhr}
\alias{fbvneglog}
\alias{fbvaneglog}
\title{Maximum-likelihood Fitting of Bivariate Extreme Value Distributions}
\description{
  Maximum-likelihood estimates for the logistic (symmetric and asymmetric),
  negative logistic (symmetric and asymmetric) and Husler-Reiss models,
  allowing any of the parameters (including marginal parameters) to be held
  fixed if desired. 
}
\usage{
fbvlog(x, start, \dots)
fbvalog(x, start, \dots)
fbvhr(x, start, \dots)
fbvneglog(x, start, \dots)
fbvaneglog(x, start, \dots)
}
\arguments{
  \item{x}{a matrix or data frame with two columns.}
  \item{start}{a named list giving the parameters to be optimized with
    initial values (see \bold{Details}).}
  \item{\dots}{additional parameters, either for the distribution or for
    \code{optim}. If parameters of the distribution are included they
    will be held fixed (see \bold{Details}).
    If parameters of the distribution are not included either
    here or as a named component in \code{start} they will be held fixed
    at the default values specified in the corresponding density
    function}
}
\details{
  Maximization of the log-likelihood is performed.
  The estimated standard errors are taken from the observed information,
  calculated by a numerical approximation. 

  The named components of \code{start} can be any or possibly all of
  \code{dep}, \code{asy} (a vector of length two), \code{mar1} and
  \code{mar2} (vectors of length three).
  Alternatively, \code{mar1} can be passed as three separate arguments;
  \code{loc1}, \code{scale1} and \code{shape1}.
  \code{mar2} can be similarly passed as three separate arguments;
  \code{loc2}, \code{scale2} and \code{shape2}.
  The asymmetry parameters can also be passed separately as \code{asy1}
  and \code{asy2}.
  Associating a separate name with each parameter allows any parameter
  subset to be fixed at specified values.
  All parameters to be fixed must be passed individually (e.g. to fix
  the parameters on the first margin all of \code{loc1}, \code{scale1}
  and \code{shape1} must be specified; using \code{mar1} to specify all
  three simultaneously results in an error).
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
bvdata <- rbvlog(100, dep = 0.6, mar1 = c(1.2,1.4,0.4),
mar2 = c(1.2,1.4,0.4))
fbvlog(bvdata, start = list(mar1 = c(2,1,0), mar2 = c(1,1,0), dep = 0.75), 
method="BFGS")
fbvlog(bvdata, start = list(mar1 = c(2,1,0), mar2 = c(1,1,0)), dep = 1, 
method="BFGS")
fgev(bvdata[,1], start = list(loc=1,scale=1,shape=0),
method="BFGS")$estimate
fgev(bvdata[,2], start = list(loc=1,scale=1,shape=0),
method="BFGS")$estimate
}

\keyword{distribution}

