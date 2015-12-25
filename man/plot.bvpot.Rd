\name{plot.bvpot}
\alias{plot.bvpot}

\title{Plot Diagnostics for a Bivariate POT EVD Object}
\description{
  Four plots (selectable by \code{which}) are currently provided:
  a density plot (1), a dependence function plot (2), a quantile
  curves plot (3) and a spectral density plot (4).
  Plot diagnostics for the generalized Pareto peaks-over-threshold
  margins (selectable by \code{mar} and \code{which}) are also
  available.
}
\usage{
\method{plot}{bvpot}(x, mar = 0, which = 1:4, main, ask = nb.fig <
    length(which) && dev.interactive(), grid = 50, above = FALSE,
    levels = NULL, tlty = 1, blty = 3, rev = FALSE, p = seq(0.75,
    0.95, 0.05), half = FALSE, \dots)
}
\arguments{
  \item{x}{An object of class \code{"bvpot"}.}
  \item{mar}{If \code{mar = 1} or \code{mar = 2} diagnostics
    are given for the first or second generalized Pareto
    margin respectively.}
  \item{which}{A subset of the numbers \code{1:4} selecting
    the plots to be shown. By default all are plotted.}
  \item{main}{Title of each plot. If given, should be a
    character vector with the same length as \code{which}.}
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{grid, levels}{Arguments for the density plot. The
    data is plotted with a contour plot of the bivariate density
    of the fitted model in the tail region. The density is evaluated
    at \code{grid^2} points, and contours are plotted at the values
    given in the numeric vector \code{levels}. If \code{levels} is
    \code{NULL} (the default), the routine attempts to find sensible
    values.}
  \item{above}{Logical; if \code{TRUE}, only data points above
    both marginal thresholds are plotted.}
  \item{tlty}{Line type for the lines identifying the thresholds.}
  \item{rev, blty}{Arguments to the dependence function
    plot. See \code{\link{abvevd}}.}
  \item{p}{Lower tail probabilities for the quantile curves plot.
    The plot is of the same type as given by the function
    \code{\link{qcbvnonpar}}, but applied to the parametric
    bivariate threshold model.}
  \item{half}{Argument to the spectral density plot. See
    \code{\link{hbvevd}}.}
  \item{\dots}{Other arguments to be passed through to plotting
    functions.}
}

\seealso{\code{\link{plot.bvevd}}, \code{\link{contour}},
  \code{\link{abvnonpar}}, \code{\link{qcbvnonpar}}, \code{\link{hbvevd}}}

\examples{
bvdata <- rbvevd(500, dep = 0.6, model = "log")
M1 <- fbvpot(bvdata, threshold = c(0,0), model = "log")
\dontrun{plot(M1)}
\dontrun{plot(M1, mar = 1)}
\dontrun{plot(M1, mar = 2)}
}
\keyword{hplot}
