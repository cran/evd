\name{profile.evd}
\alias{profile.evd}

\title{Method for Profiling evd Objects}
\description{
  Investigates behavior of objective function near the solution
  represented by \code{fitted}.
}
\usage{
\method{profile}{evd}(fitted, which = names(fitted$estimate), conf = 0.999,
    mesh = fitted$std.err[which]/2, xmin = rep(-Inf, length(which)),
    xmax = rep(Inf, length(which)), convergence = FALSE, method = "BFGS",
    control = list(maxit = 500), \dots)
}
\arguments{
  \item{fitted}{The original fitted model object.}
  \item{which}{The original model parameters that are to be profiled.
    By default, all parameters are profiled.}
  \item{conf}{Controls the range over which the parameters are profiled.
    The confidence coefficient for the confidence interval defined by
    this range (assuming the usual asymptotic properties) will be at
    least as large as \code{conf}, unless theoretical limits are set
    using \code{xmin} or/and \code{xmax}.}
  \item{mesh}{A numeric vector containing one value for each
    parameter in \code{which}. The values represent the
    distance between the points profiled. By default half the
    standard errors are used. If the fitted object does not contain
    standard errors the argument must be specified.
    The argument should also be specified when an estimator is
    on or close to a parameter boundary, since the approximated
    ``standard error'' will then be close to zero.}
  \item{xmin, xmax}{Numeric vectors containing one value for each
    parameter in \code{which}. Each value represents the theoretical
    lower/upper bound of the corresponding parameter.
    The arguments are typically needed only when a parameter has a
    theoretical lower/upper bound that may be contained in the
    range over which the parameters are profiled (as specified via
    \code{conf}). Do not use these arguments to specify plotting
    ranges in a subsequent plot; to do this use \code{xlim} in
    the call to \code{plot}.}
  \item{convergence}{Logical; print convergence code after each
    optimization?}
  \item{method}{The optimization method.}
  \item{control}{Passed to \code{optim}. See \code{\link{optim}} for
    details.}
  \item{\dots}{Ignored.}
}
\value{
  An object of class \code{"profile.evd"}; a list with an element for
  each parameter being profiled. The
  elements are matrices. The first column contains the values
  of the profiled parameter. The second column contains profile
  deviances. The remaining columns contain the constrained maximum
  likelihood estimates for the remaining model parameters.
}

\seealso{\code{\link{profile}}, \code{\link{profile2d}},
  \code{\link{profile2d.evd}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{M1P <- profile(M1)}
\dontrun{pcint(M1P)}
\dontrun{plot(M1P)}
}
\keyword{models}



