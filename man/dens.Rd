\name{dens}
\alias{dens}

\title{A Density Plot for an evd Object}
\description{
  A density plot for an evd object. The density of the fitted
  model is plotted with a rug plot and (optionally) a
  non-parameteric estimate.
}
\usage{
dens(x, adjust = 1, nplty = 2, jitter = FALSE, main = "Density Plot", 
    xlab = "Quantile", ylab = "Density", \dots)
}
\arguments{
  \item{x}{An object of class \code{"evd"}.}
  \item{adjust}{Adjust the smoothing bandwidth for the non-parametric
     estimate. See \code{\link{density}} for details.}
  \item{nplty}{The line type of the non-parametric estimate.
    Use zero to suppress.}
  \item{jitter}{Logical; if \code{TRUE}, the (transformed) data
    are jittered to produce the rug plot. See \code{\link{jitter}}.
    This need only be used if the data contains repeated
    values.}
  \item{main}{Title of plot.}
  \item{xlab,ylab}{Labels for x and y axes.}
  \item{\dots}{Other plot parameters.}
}

\details{
  For non-stationary models the data are transformed to stationarity.
  The plot then corresponds to the distribution obtained when all
  covariates are zero.
}

\seealso{\code{\link{density}}, \code{\link{jitter}},
  \code{\link{plot.evd}}, \code{\link{pp}}, \code{\link{qq}},
  \code{\link{rl}}, \code{\link{rug}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{dens(M1)}
}
\keyword{hplot}



