\name{fgumbelx}
\alias{fgumbelx}

\title{Maximum-likelihood Fitting of the Maximum of Two Gumbel 
  Distributions}
\description{
  Maximum-likelihood fitting for the maximum of two gumbel distributions,
  allowing any of the parameters to be held fixed if desired.
}
	
\usage{
fgumbelx(x, start, \dots, nsloc1 = NULL, nsloc2 = NULL, std.err = TRUE, 
    corr = FALSE, method = "BFGS", warn.inf = TRUE)
}
\arguments{
  \item{x}{A numeric vector, which may contain missing values.}
  \item{start}{A named list giving the initial values for the
    parameters over which the likelihood is to be maximized.
    If \code{start} is omitted the routine attempts to find good
    starting values using moment estimators.}
  \item{\dots}{Additional parameters, either for the fitted model
    or for the optimization function \code{optim}. If parameters
    of the model are included they will be held fixed at the
    values given (see \bold{Examples}).}
  \item{nsloc1}{A data frame with the same number of rows as the
    length of \code{x}, for linear modelling of the location
	parameter of the first Gumbel distribution. This is not 
	recommended as the model is already complex.}
  \item{nsloc2}{A data frame with the same number of rows as the
    length of \code{x}, for linear modelling of the location
	parameter of the second Gumbel distribution. This is not 
	recommended as the model is already complex.}
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

  For stationary models the parameter names are \code{loc1}, \code{scale1},
  \code{loc2} and \code{scale2} for the location and scale parameters of
  two Gumbel distributions, where \code{loc2} must be greater or equal to 
  \code{loc1}. 

  The likelihood may have multiple local optima and therefore may be 
  difficult to fit properly; the default starting values use a moment 
  based approach, however it is recommended that the user specify multiple 
  different starting values and experiment with different optimization 
  methods.
}

\value{
  Returns an object of class \code{c("gumbelx","evd")}.

  The generic accessor functions \code{\link{fitted}} (or
  \code{\link{fitted.values}}), \code{\link{std.errors}},
  \code{\link{deviance}}, \code{\link{logLik}} and
  \code{\link{AIC}} extract various features of the
  returned object.

  Profile likelihoods and diagnostic plots are not currently
  implemented for this model.
		
  An object of class \code{c("gumbelx","evd")} is a list
  containing at most the following components
  \item{estimate}{A vector containing the maximum likelihood estimates.}
  \item{std.err}{A vector containing the standard errors.}
  \item{fixed}{A vector containing the parameters of the model that
    have been held fixed.}
  \item{param}{A vector containing all parameters (optimized and fixed).}
  \item{deviance}{The deviance at the maximum likelihood estimates.}
  \item{corr}{The correlation matrix.}
  \item{var.cov}{The variance covariance matrix.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}}.}
  \item{data}{The data passed to the argument \code{x}.}
  \item{nsloc1}{The argument \code{nsloc1}.}
  \item{nsloc2}{The argument \code{nsloc2}.}
  \item{n}{The length of \code{x}.}
  \item{call}{The call of the current function.}
}
\section{Warning}{
  This function is experimental and involves optimizing over a potentially
  complex surface.
}

\seealso{\code{\link{fgev}}, \code{\link{optim}}}

\examples{
uvdata <- pmax(rgumbel(100, loc = 0, scale = 1), rgumbel(100, loc = 1, scale = 1))
fgumbelx(uvdata, loc1 = 0, scale1 = 1)
}
\keyword{models}

