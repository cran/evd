\name{plot.gev}
\alias{plot.gev}

\title{Plot Diagnostics for a GEV Object}
\description{
  Four plots (selectable by \code{which}) are currently provided:
  a P-P plot, a Q-Q plot, a density plot and a return level plot.
}
\usage{
\method{plot}{gev}(x, which = 1:4, main = c("Probability Plot",
    "Quantile Plot", "Density Plot", "Return Level Plot"),
    ask = nb.fig < length(which) && dev.interactive(),
    ci = TRUE, adjust = 1, jitter = FALSE, nplty = 2, \dots)
}
\arguments{
  \item{x}{An object of class \code{"gev"}.}
  \item{which}{If a subset of the plots is required, specify a
    subset of the numbers \code{1:4}.}
  \item{main}{Title of each plot.}
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{ci}{Logical; if \code{TRUE} (the default), plot simulated
    95\% confidence intervals for the P-P, Q-Q and return level
    plots.}
  \item{adjust, jitter, nplty}{Arguments to the density plot.
    The density of the fitted model is plotted with a rug plot and
    (optionally) a non-parameteric estimate. The argument
    \code{adjust} controls the smoothing bandwidth for the
    non-parametric estimate (see \code{\link{density}}).
    \code{jitter} is logical; if \code{TRUE}, the (possibly
    transformed) data are jittered to produce the rug plot.
    This need only be used if the data contains repeated
    values. \code{nplty} is the line type of the non-parametric
    estimate. To omit the non-parametric estimate set \code{nplty}
    to zero.}
  \item{\dots}{Other parameters to be passed through to plotting
    functions.}
}
\details{
  The following discussion assumes that the fitted model is
  stationary. For non-stationary models the data are transformed
  to stationarity. The plot then corresponds to the distribution
  obtained when all covariates are zero.

  The P-P plot consists of the points
  \deqn{\{(G_n(z_i), G(z_i)), i = 1,\ldots,m\}}{
    {(G_n(z_i), G(z_i)), i = 1,\ldots,m}}
  where \eqn{G_n} is the empirical
  distribution function (defined using \code{\link{ppoints}}),
  G is the model based estimate of the generalized extreme
  value distribution, and \eqn{z_1,\ldots,z_m} are the data
  used in the fitted model, sorted into ascending order.

  The Q-Q plot consists of the points
  \deqn{\{(G^{-1}(p_i), z_i), i = 1,\ldots,m\}}{
    {(G^{-1}(p_i), z_i), i = 1,\ldots,m}}
  where \eqn{G^{-1}} is the model based
  estimate of the generalized extreme value quantile function,
  \eqn{p_1,\ldots,p_m} are plotting points defined by
  \code{\link{ppoints}}, and \eqn{z_1,\ldots,z_m} are the data
  used in the fitted model, sorted into ascending order.

  The return level plot is defined as follows.
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
  estimates of \eqn{(a,b,s)}. If the estimate of \eqn{s} is zero, the
  curve will be linear.
  For large values of \eqn{t}, \eqn{y_t} is approximately equal
  to the return period \eqn{t}. It is usual practice to label the
  x-axis as the return period.  

  The points on the plot are  
  \deqn{\{(-1/\log(p_i), z_i), i = 1,\ldots,m\}}{
    {(-1/log(p_i), z_i), i = 1,\ldots,m}}
  where \eqn{p_1,\ldots,p_m} are plotting points defined by
  \code{\link{ppoints}}, and \eqn{z_1,\ldots,z_m} are the data
  used in the fitted model, sorted into ascending order.
  For a good fit the points should lie ``close'' to the curve
  defined by \eqn{(z_t,\log(y_t))}{(z_t,\log(y_t))}.
}

\seealso{\code{\link{plot.bvevd}}, \code{\link{density}},
  \code{\link{jitter}}, \code{\link{rug}}, \code{\link{ppoints}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{par(mfrow = c(2,2))}
\dontrun{plot(M1)}
}
\keyword{hplot}
