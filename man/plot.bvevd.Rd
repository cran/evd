\name{plot.bvevd}
\alias{plot.bvevd}

\title{Plot Diagnostics for a Bivariate EVD Object}
\description{
  Four plots (selectable by \code{which}) are currently provided:
  two conditional P-P plots (conditioning on each margin),
  a density plot and a dependence function plot.
  Plot diagnostics for the generalized extreme value margins
  (selectable by \code{mar} and \code{which}) are also available.
}
\usage{
\method{plot}{bvevd}(x, mar = 0, which = 1:4, main = c("Conditional Plot One", 
    "Conditional Plot Two", "Density Plot", "Dependence Function"), 
    ask = nb.fig < length(which) && dev.interactive(), ci = TRUE, 
    jitter = FALSE, grid = 50, nplty = 2, blty = 3, method = "cfg",
    convex = FALSE, wf = function(t) t, \dots)
}
\arguments{
  \item{x}{An object of class \code{"bvevd"}.}
  \item{mar}{If \code{mar = 1} or \code{mar = 2} diagnostics
    are given for the first or second genereralized extreme
    value margin respectively. The values of the remaining
    parameters are then passed to the plot method
    \code{\link{plot.gev}}.}
  \item{which}{If a subset of the plots is required, specify a
    subset of the numbers \code{1:4}.}
  \item{main}{Title of each plot.}
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{ci}{Logical; if \code{TRUE} (the default), plot simulated
    95\% confidence intervals for the conditional P-P plots.}
  \item{jitter, grid}{Arguments for the density plot. The (possibly
    transformed) data is plotted with a contour plot of the
    bivariate density of the fitted model. The density is evaluated
    at \code{grid^2} points. If \code{jitter} is \code{TRUE}, the
    data are jittered. This need only be used if the data contains
    repeated values.}
  \item{nplty, blty, method, convex, wf}{Arguments to the dependence
    function plot. The dependence function for the fitted model is
    plotted and (optionally) compared to a non-parameteric estimate.
    See \code{\link{abvnonpar}} for a definition of the dependence
    function, and for a description of the arguments \code{method},
    \code{modify} and \code{wf}, which alter the behaviour of the
    non-parametric estimator. \code{nplty} is the line type of the
    non-parametric estimate. To omit the non-parametric estimate set
    \code{nplty} to zero. \code{blty} is the line type of the
    triangular border. To omit the border estimate set \code{blty}
    to zero.}
  \item{\dots}{Other arguments to be passed through to plotting
    functions.}
}
\details{
  The following discussion assumes that the fitted model is
  stationary. For non-stationary models the data are transformed
  to stationarity. The plot then corresponds to the distribution
  obtained when all covariates are zero.
  
  A conditional P-P plot is a P-P plot for the condition
  distribution function of a bivariate evd object.
  Let \eqn{G(.|.)} be the conditional distribution of
  the first margin given the second, under the fitted model.
  Let \eqn{z_1,\ldots,z_m} be the data used in the fitted model,
  where \eqn{z_j = (z_{1j}, z_{2j})} for \eqn{j = 1,\ldots,m}.
  The plot that (by default) is labelled Conditional Plot Two,
  conditioning on the second margin, consists of the points
  \deqn{\{(p_i, c_i), i = 1,\ldots,m\}}{
    {(p_i, c_i), i = 1,\ldots,m}}
  where \eqn{p_1,\ldots,p_m} are plotting points defined by
  \code{\link{ppoints}} and \eqn{c_i} is the \eqn{i}th largest
  value from the sample
  \eqn{\{G(z_{j1}|z_{j2}), j = 1,\ldots,m\}.}{
    {G(z_{j1}|z_{j2}), j = 1,\ldots,m}.}
  The margins are reversed for Conditional Plot One, so that
  \eqn{G(.|.)} is the conditional distribution of the second
  margin given the first. 
}

\seealso{\code{\link{plot.gev}}, \code{\link{contour}},
  \code{\link{jitter}}, \code{\link{abvnonpar}}}

\examples{
bvdata <- rbvevd(100, dep = 0.6, model = "log")
M1 <- fbvevd(bvdata, model = "log")
\dontrun{par(mfrow = c(2,2))}
\dontrun{plot(M1)}
\dontrun{plot(M1, mar = 1)}
\dontrun{plot(M1, mar = 2)}
}
\keyword{hplot}
