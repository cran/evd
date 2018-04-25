\name{plot.bvevd}
\alias{plot.bvevd}

\title{Plot Diagnostics for a Bivariate EVD Object}
\description{
  Six plots (selectable by \code{which}) are currently provided:
  two conditional P-P plots (1,2), conditioning on each margin, a
  density plot (3), a dependence function plot (4), a quantile
  curves plot (5) and a spectral density plot (6).
  Plot diagnostics for the generalized extreme value margins
  (selectable by \code{mar} and \code{which}) are also available.
}
\usage{
\method{plot}{bvevd}(x, mar = 0, which = 1:6, main, ask = nb.fig <
    length(which) && dev.interactive(), ci = TRUE, cilwd = 1,
    a = 0, grid = 50, legend = TRUE, nplty = 2, blty = 3, method = "cfg",
    convex = FALSE, rev = FALSE, p = seq(0.75, 0.95, 0.05),
    mint = 1, half = FALSE, \dots)
}
\arguments{
  \item{x}{An object of class \code{"bvevd"}.}
  \item{mar}{If \code{mar = 1} or \code{mar = 2} diagnostics
    are given for the first or second genereralized extreme
    value margin respectively.}
  \item{which}{A subset of the numbers \code{1:6} selecting
    the plots to be shown. By default all are plotted.}
  \item{main}{Title of each plot. If given, should be a
    character vector with the same length as \code{which}.}
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{ci}{Logical; if \code{TRUE} (the default), plot simulated
    95\% confidence intervals for the conditional P-P plots.}
  \item{cilwd}{Line width for confidence interval lines.}
  \item{a}{Passed through to \code{ppoints} for empirical estimation.
    Larger values give less probability for extreme events.}
  \item{grid}{Argument for the density plot. The (possibly
    transformed) data is plotted with a contour plot of the
    bivariate density of the fitted model. The density is evaluated
    at \code{grid^2} points.}
  \item{legend}{If \code{legend} is \code{TRUE} and if the
    fitted data contained a third column of mode \code{logical},
    then a legend is included in the density and quantile curve
    plots.}
  \item{method, convex, rev}{Arguments to the dependence function
    plot. The dependence function for the fitted model is plotted and
    (optionally) compared to a non-parameteric estimate. See
    \code{\link{abvnonpar}} for a description of the arguments.}
  \item{nplty, blty}{Line types for the dependence function plot.
    \code{nplty} is the line type of the non-parametric estimate.
    To omit the non-parametric estimate set \code{nplty} to zero.
    \code{blty} is the line type of the triangular border. To omit
    the border estimate set \code{blty} to zero.}
  \item{p, mint}{Arguments to the quantile curves plot. See
    \code{\link{qcbvnonpar}} for a description of the plot and
    the arguments.}
  \item{half}{Argument to the spectral density plot. See
    \code{\link{hbvevd}}.}
  \item{\dots}{Other arguments to be passed through to plotting
    functions.}
}
\details{
  In all plots we assume that the fitted model is
  stationary. For non-stationary models the data are transformed
  to stationarity. The plot then corresponds to the distribution
  obtained when all covariates are zero. In particular, the
  density and quanitle curves plots will not plot the original
  data for non-stationary models.
  
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

\seealso{\code{\link{plot.uvevd}}, \code{\link{contour}},
  \code{\link{jitter}}, \code{\link{abvnonpar}},
  \code{\link{qcbvnonpar}}}

\examples{
bvdata <- rbvevd(100, dep = 0.6, model = "log")
M1 <- fbvevd(bvdata, model = "log")
\dontrun{par(mfrow = c(2,2))}
\dontrun{plot(M1, which = 1:5)}
\dontrun{plot(M1, mar = 1)}
\dontrun{plot(M1, mar = 2)}
}
\keyword{hplot}
