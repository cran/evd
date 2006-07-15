\name{tcplot}
\alias{tcplot}
\title{Threshold Choice Plot}
\description{
  Plots of parameter estimates at various thresholds for peaks
  over threshold modelling, using the Generalized Pareto or Point
  Process representation. 
}
\usage{
tcplot(data, tlim, model = c("gpd","pp"), pscale = FALSE, cmax =
    FALSE, r = 1, ulow = -Inf, rlow = 1, nt = 25, which = 1:npar,
    conf = 0.95, lty = 1, lwd = 1, type = "b", cilty = 1, vci =
    TRUE, xlab, xlim, ylabs, ylims, ask = nb.fig < length(which) &&
    dev.interactive(), \dots)
}
\arguments{
  \item{data}{A numeric vector.}
  \item{tlim}{A numeric vector of length two, giving the limits for
    the thresholds at which the model is fitted.}
  \item{model}{The model; either \code{"gpd"} (the default) or
    \code{"pp"}, for the Generalized Pareto or Point Process
    representations respectively.}
  \item{pscale}{If \code{TRUE}, then the x-axis gives the threshold
    exceedance probability rather than the threshold itself.}
  \item{cmax}{Logical; if \code{FALSE} (the default), the models
    are fitted using all exceedences over the thresholds. If
    \code{TRUE}, the models are fitted using cluster maxima, using
    clusters of exceedences derived from \code{clusters}.}
  \item{r, ulow, rlow}{Arguments used for the identification of
    clusters of exceedences (see \code{\link{clusters}}).
    Ignored if \code{cmax} is \code{FALSE} (the default).}
  \item{nt}{The number of thresholds at which the model is fitted.}
  \item{which}{If a subset of the plots is required, specify a
    subset of the numbers \code{1:npar}, where \code{npar} is
    the number of parameters, so that \code{npar = 2} when
    \code{model = "gpd"} (the default) and \code{npar = 3}
    when \code{model = "pp"}.}
  \item{conf}{The (pointwise) confidence coefficient for the
    plotted confidence intervals. Use zero to suppress.}
  \item{lty, lwd}{The line type and width of the line connecting
    the parameter estimates.}
  \item{type}{The form taken by the line connecting the parameter
    estimates and the points denoting these estimates. Possible
    values include \code{"b"} (the default) for points joined by
    lines, \code{"o"} for overplotted points and lines, and
    \code{"l"} for an unbroken line with no points.}
  \item{cilty}{The line type of the lines depicting the
    confidence intervals.}
  \item{vci}{If \code{TRUE} (the default), confidence intervals
    are plotted using vertical lines.}
  \item{xlab, xlim}{Label and limits for the x-axis; if given,
    these arguments apply to every plot.}
  \item{ylabs, ylims}{A vector of y-axis labels and a matrix of
    y-axis limits. If given, \code{ylabs} should have the same
    length as \code{which}, and \code{ylims} should have two
    columns and \code{length(which)} rows. If the length of
    \code{which} is one, then \code{ylims} can be a vector of
    length two.} 
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{\dots}{Other arguments to be passed to the model fit
    function \code{fpot}.}
}
\details{
  For each of the \code{nt} thresholds a peaks over threshold model
  is fitted using the function \code{fpot}. When \code{model = "gpd"}
  (the default), the maximum likelihood estimates for the shape and the
  modified scale parameter (modified by subtracting the shape multiplied
  by the threshold) are plotted against the thresholds.
  When \code{model = "pp"} the maximum likelihood estimates for the
  location, scale and shape parameters are plotted against the
  thresholds. (The modified scale parameter in the \code{"gpd"} case
  is equivalent to the scale parameter in the \code{"pp"} case.)
  If the threshold \code{u} is a valid threshold to be used for peaks
  over threshold modelling, the parameter estimates depicted should
  be approximately constant above \code{u}.
}
\value{
  A list is invisibly returned. Each component is a matrix with
  three columns giving parameter estimates and confidence limits.
}
\author{Stuart Coles and Alec Stephenson}

\seealso{\code{\link{fpot}}, \code{\link{mrlplot}},
    \code{\link{clusters}}}

\examples{
data(portpirie)
tlim <- c(3.6, 4.2)
\dontrun{tcplot(portpirie, tlim)}
\dontrun{tcplot(portpirie, tlim, nt = 100, lwd = 3, type = "l")}
\dontrun{tcplot(portpirie, tlim, model = "pp")}
}
\keyword{hplot}

