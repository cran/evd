\name{plot.evd}
\alias{plot.evd}

\title{Plot Diagnostics for an evd Object}
\description{
  Four plots (selectable by \code{which}) are currently provided:
  a P-P plot, a Q-Q plot, a density plot and a return level plot.
}
\usage{
\method{plot}{evd}(x, which = 1:4, main = c("Probability Plot",
    "Quantile Plot", "Density Plot", "Return Level Plot"),
    ask = nb.fig < length(which) && dev.interactive(),
    ci = TRUE, adjust = 1, jitter = FALSE, nplty = 2, \dots)
}
\arguments{
  \item{x}{An object of class \code{"evd"}.}
  \item{which}{If a subset of the plots is required, specify a
    subset of the numbers \code{1:4}.}
  \item{main}{Title of each plot.}
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{ci}{Logical; if \code{TRUE} (the default), plot simulated
    95\% confidence intervals for the P-P, Q-Q and return level
    plots.}
  \item{adjust,jitter,nplty}{Arguments to the density plot.
    See \code{\link{dens}} for details.}
  \item{\dots}{Other parameters to be passed through to plotting functions.}
}
\details{
  The plot method calls the functions \code{\link{pp}},
  \code{\link{qq}}, \code{\link{dens}} and \code{\link{rl}}.
}

\seealso{\code{\link{dens}}, \code{\link{plot.bvevd}}, \code{\link{pp}},
  \code{\link{qq}}, \code{\link{rl}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{plot(M1)}
}
\keyword{hplot}