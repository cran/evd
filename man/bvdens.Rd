\name{bvdens}
\alias{bvdens}

\title{A Density Plot for a Bivariate evd Object}
\description{
  The density of the fitted model is plotted using
  \code{\link{contour}}.
}
\usage{
bvdens(x, jitter = FALSE, nlevels = 10, levels, main = "Density Plot",
    xlab = "", ylab = "", \dots)
}
\arguments{
  \item{x}{An object of class \code{"bvevd"}.}
  \item{jitter}{Logical; if \code{TRUE}, the (transformed) data
    are jittered. See \code{\link{jitter}}.
    This need only be used if the data contains repeated
    values.}
  \item{nlevels,levels}{Arguments passed to \code{\link{contour}}.}
  \item{main}{Title of plot.}
  \item{xlab,ylab}{Labels for x and y axes.}
  \item{\dots}{Other parameters to be passed to \code{contour}.}
}

\details{
  For non-stationary models the data are transformed to stationarity.
  The plot then corresponds to the distribution obtained when all
  covariates are zero.
}

\seealso{\code{\link{bvcpp}}, \code{\link{bvdp}}, \code{\link{contour}},
  \code{\link{jitter}}, \code{\link{plot.bvevd}}}

\examples{
bvdata <- rbvlog(100, dep = 0.6)
M1 <- fbvlog(bvdata)
\dontrun{bvdens(M1)}
}
\keyword{hplot}
