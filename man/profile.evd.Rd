\name{profile.evd}
\alias{profile.evd}

\title{Method for Profiling EVD Objects}
\description{
  Calculate profile traces for fitted models.
}
\usage{
\method{profile}{evd}(fitted, which = names(fitted$estimate), conf = 0.999,
    mesh = fitted$std.err[which]/4, xmin = rep(-Inf, length(which)),
    xmax = rep(Inf, length(which)), convergence = FALSE, method = "BFGS",
    control = list(maxit = 500), \dots)
}
\arguments{
  \item{fitted}{An object of class \code{"evd"}.}
  \item{which}{A character vector giving the model parameters that
    are to be profiled. By default, all parameters are profiled.}
  \item{conf}{Controls the range over which the parameters are profiled.
    The profile trace is constructed so that (assuming the usual
    asymptotic properties hold) profile confidence intervals with
    confidence coefficients \code{conf} or less can be derived from it.}
  \item{mesh}{A numeric vector containing one value for each
    parameter in \code{which}. The values represent the
    distance between the points profiled. By default \code{mesh} is
    one quarter of the standard errors. If the fitted object does not
    contain standard errors the argument must be specified.
    The argument should also be specified when an estimator is
    on or close to a parameter boundary, since the approximated
    ``standard error'' will then be close to zero.}
  \item{xmin, xmax}{Numeric vectors containing one value for each
    parameter in \code{which}. Each value represents the theoretical
    lower/upper bound of the corresponding parameter.
    The arguments are needed only when a parameter has a
    lower/upper bound at which the likelihood is non-zero. Do not
    use these arguments to specify plotting ranges in a subsequent
    plot (as they are used in the calculation of profile confidence
    intervals); to do this use \code{xlim} in the call to \code{plot}.}
  \item{convergence}{Logical; print convergence code after each
    optimization? (A warning is given for each non-zero convergence
    code, irrespective of the value of \code{convergence}.)}
  \item{method}{The optimization method.}
  \item{control}{Passed to \code{optim}. See \code{\link{optim}} for
    details.}
  \item{\dots}{Ignored.}
}
\value{
  An object of class \code{"profile.evd"}, which is a list with an
  element for each parameter being profiled. The elements are
  matrices. The first column contains the values of the profiled
  parameter. The second column contains profile deviances. The
  remaining columns contain the constrained maximum likelihood
  estimates for the remaining model parameters. For calculation of 
  profile confidence intervals, use the \code{\link{confint.profile.evd}} 
  function.
}

\seealso{\code{\link{confint.profile.evd}}, \code{\link{profile2d.evd}}, 
  \code{\link{plot.profile.evd}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{M1P <- profile(M1)}
\dontrun{par(mfrow = c(2,2))}
\dontrun{cint <- plot(M1P)}
\dontrun{cint}
}
\keyword{models}



