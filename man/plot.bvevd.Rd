\name{plot.bvevd}
\alias{plot.bvevd}

\title{Plot Diagnostics for a Bivariate evd Object}
\description{
  Four plots (selectable by \code{which}) are currently provided:
  two conditional P-P plots (conditioning on each margin),
  a density plot and a dependence function plot.
  Plot diagnostics for the generalized extreme value margins
  (selectable by \code{mar} and \code{which}) are also provided.
}
\usage{
\method{plot}{bvevd}(x, mar = 0, which = 1:4, main = c("Conditional Plot One", 
    "Conditional Plot Two", "Density Plot", "Dependence Function"), 
    ask = nb.fig < length(which) && dev.interactive(), ci = TRUE, 
    nlevels = 10, levels, jitter = FALSE, nplty = 2, method = "cfg", 
    modify = 0, wf = function(t) t, \dots)
}
\arguments{
  \item{x}{An object of class \code{"bvevd"}.}
  \item{mar}{If \code{mar = 1} or \code{mar = 2} diagnostics
    are given for the first or second genereralized extreme
    value margin. The values of the remaining parameters are
    then passed to the plot method \code{\link{plot.evd}}.}
  \item{which}{If a subset of the plots is required, specify a
    subset of the numbers \code{1:4}.}
  \item{main}{Title of each plot.}
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{ci}{Logical; if \code{TRUE} (the default), plot simulated
    95\% confidence intervals for the conditional P-P plots.}
  \item{nlevels,levels,jitter}{Arguments to the density plot.
    See \code{\link{bvdens}} for details.}
  \item{nplty,method,modify,wf}{Arguments to the dependence function
    plot. See \code{\link{bvdp}} for details.}
  \item{\dots}{Other parameters to be passed through to plotting
    functions.}
}
\details{
  The plot method calls the functions \code{\link{bvcpp}},
  \code{\link{bvdens}} and \code{\link{bvdp}}.
}

\seealso{\code{\link{bvcpp}}, \code{\link{bvdens}},
  \code{\link{bvdp}}, \code{\link{plot.evd}}}

\examples{
bvdata <- rbvlog(100, dep = 0.6)
M1 <- fbvlog(bvdata)
\dontrun{plot(M1)}
\dontrun{plot(M1, mar = 1)}
\dontrun{plot(M1, mar = 2)}
}
\keyword{hplot}
