\name{rl}
\alias{rl}

\title{A Return Level Plot for an evd Object}
\description{
  A return level plot for an evd object.
}
\usage{
pp(x, ci = TRUE, main = "Quantile Plot", xlab = "Model", 
    ylab = "Empirical", \dots)
}
\arguments{
  \item{x}{An object of class \code{"evd"}.}
  \item{ci}{Logical; if \code{TRUE} (the default), plot simulated
    95\% confidence intervals.}
  \item{main}{Title of plot.}
  \item{xlab,ylab}{Labels for x and y axes.}
  \item{\dots}{Other plot parameters.}
}
\details{
  Let \eqn{G} be the generalized extreme value distribution function,
  with location, scale and shape parameters \eqn{a}, \eqn{b}
  and \eqn{s} respectively.
  Let \eqn{z_t} be defined by \eqn{G(z_t) = 1 - 1/t}.
  In common terminology, \eqn{z_t} is the return level
  associated with the return period \eqn{t}.

  Let \eqn{y_t = -1/\log(1 - 1/t)}{y_t = -1/log(1 - 1/t)}.
  It follows that
  \deqn{z_t = a + b(y_t^s - 1)/s.}{
    z_t = a + b((y_t)^s - 1)/s.}
  When \eqn{s = 0}, \eqn{z_t} is defined by continuity, so that
  \deqn{z_t = a + b\log(y_t).}{
    z_t = a + b log(y_t).}
  The curve within the return level plot is \eqn{z_t} plotted
  against \eqn{y_t} on a logarithmic scale, using maximum likelihood
  estimates of \eqn{(a,b,s)}. If the estimator of \eqn{s} is zero, the
  curve will be linear.
  For large values of \eqn{t}, \eqn{y_t} is approximately equal
  to \eqn{t}.

  The points on the plot are  
  \deqn{\{(-1/\log(p_i), z_i), i = 1,\ldots,m\}}{
    {(-1/log(p_i), z_i), i = 1,\ldots,m}}
  where \eqn{p_1,\ldots,p_m} are plotting points defined by
  \code{\link{ppoints}}, and \eqn{z_1,\ldots,z_m} are the data
  used in the fitted model, sorted into ascending order.
  For a good fit the points should lie ``close'' to the curve
  defined by \eqn{(z_t,\log(y_t))}{(z_t,\log(y_t))}.

  For non-stationary models the data are transformed to stationarity.
  The plot then corresponds to the distribution obtained when all
  covariates are zero.
}

\seealso{\code{\link{dens}}, \code{\link{plot.evd}},
  \code{\link{ppoints}}, \code{\link{pp}}, \code{\link{qq}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{rl(M1)}
}
\keyword{hplot}