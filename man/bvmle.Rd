\name{bvmle}
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
  Maximum-likelihood estimates for the symmetric/asymmetric logistic,
  symmetric/asymmetric negative logistic, bilogistic, negative
  bilogistic, Husler-Reiss and Coles-Tawn models, including linear
  modelling of the marginal location parameters, and allowing any of
  the parameters to be held fixed if desired. 
}
\usage{
fbvlog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
    method = "BFGS")
fbvalog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
    method = "BFGS")
fbvhr(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
    method = "BFGS")
fbvneglog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
    method = "BFGS")
fbvaneglog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE,
    method = "BFGS")
fbvbilog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE, 
    method = "BFGS") 
fbvnegbilog(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE, 
    method = "BFGS")
fbvct(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE, 
    method = "BFGS") 
}
\arguments{
  \item{x}{A matrix or data frame with two columns, which may contain
    missing values.}
  \item{start}{A named list giving the parameters to be optimized with
    initial values (see \bold{Details}).
    If \code{start} is omitted the routine attempts to find good
    starting values using marginal maximum likelihood estimators.}
  \item{\dots}{Additional parameters, either for the distribution or for
    \code{optim}. If parameters of the distribution are included they
    will be held fixed (see \bold{Details}).}
  \item{nsloc1, nsloc2}{A data frame with the same number of rows as the
    length of \code{x}, for linear modelling of the location
    parameter on the first/second margin (see \bold{Details}).
    The data frames are treated as covariate matrices, excluding the
    intercept.
    A numeric vector can be given as an alternative to a single column
    data frame.}
  \item{std.err}{Logical; if \code{TRUE} (the default), the ``standard
    errors'' are returned.}
  \item{method}{The optimization method (see \code{\link{optim}} for
    details).}
}
\details{
  Maximization of the log-likelihood is performed.
  The \code{std.err} component of the returned list is taken from the
  observed information, calculated by a numerical approximation.
  The ``standard errors'' must be interpreted with caution because
  the usual asymptotic properties of maximum likelihood estimators
  may not hold (Smith, 1985).

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
  A list with with components
  \item{estimate}{A vector containing the maximum likelihood estimates.}
  \item{std.err}{A vector containing the ``standard errors'' (if
    \eqn{\code{std.err} = \code{TRUE}}).}
  \item{deviance}{The deviance at the maximum likelihood estimates.}
  \item{counts}{A two-element integer vector giving the number of calls
    to the log-likelihood function and the number gradient approximations
    respectively.}
}
\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}

\seealso{\code{\link{optim}}}

\examples{
bvdata <- rbvlog(100, dep = 0.6, mar1 = c(1.2,1.4,0.4),
mar2 = c(1.2,1.4,0.4))
fbvlog(bvdata)
fbvlog(bvdata, dep = 1)
fgev(bvdata[,1])$estimate
fgev(bvdata[,2])$estimate
fbvlog(bvdata, nsloc1 = (-49:50)/100)
fbvlog(bvdata, nsloc1 = (-49:50)/100, nsloc2 = data.frame(trend =
(-49:50)/100, random = runif(100, min=-.5, max=.5)))
fbvlog(bvdata, nsloc1 = (-49:50)/100, nsloc2 = data.frame(trend =
(-49:50)/100, random = runif(100, min=-.5, max=.5)), loc2random = 0)
}

\keyword{distribution}

