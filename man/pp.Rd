\name{pp}
\alias{pp}

\title{A P-P Plot for an evd Object}
\description{
  A P-P plot for an evd object.
}
\usage{
pp(x, ci = TRUE, main = "Probability Plot", xlab = "Empirical", 
    ylab = "Model", \dots)
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
  The P-P plot constists of the points
  \deqn{\{(G_n(z_i), G(z_i)), i = 1,\ldots,m\}}{
    {(G_n(z_i), G(z_i)), i = 1,\ldots,m}}
  where \eqn{G_n} is the empirical
  distribution function (defined using \code{\link{ppoints}}),
  G is the model based estimate of the generalized extreme
  value distribution, and \eqn{z_1,\ldots,z_m} are the data
  used in the fitted model, sorted into ascending order.

  For non-stationary models the data are transformed to stationarity.
  The plot then corresponds to the distribution obtained when all
  covariates are zero.
}

\seealso{\code{\link{dens}}, \code{\link{plot.evd}},
  \code{\link{ppoints}}, \code{\link{qq}}, \code{\link{rl}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{pp(M1)}
}
\keyword{hplot}