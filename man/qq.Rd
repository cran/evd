\name{qq}
\alias{qq}

\title{A Q-Q Plot for an evd Object}
\description{
  A Q-Q plot for an evd object.
}
\usage{
pp(x, ci = TRUE, main = "Quantile Plot", xlab = "Model", 
    ylab = "Empirical", \dots)
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
  The Q-Q plot consists of the points
  \deqn{\{(G^{-1}(p_i), z_i), i = 1,\ldots,m\}}{
    {(G^{-1}(p_i), z_i), i = 1,\ldots,m}}
  where \eqn{G^{-1}} is the model based
  estimate of the generalized extreme value quantile function,
  \eqn{p_1,\ldots,p_m} are plotting points defined by
  \code{\link{ppoints}}, and \eqn{z_1,\ldots,z_m} are the data
  used in the fitted model, sorted into ascending order.

  For non-stationary models the data are transformed to stationarity.
  The plot then corresponds to the distribution obtained when all
  covariates are zero.
}

\seealso{\code{\link{dens}}, \code{\link{plot.evd}},
  \code{\link{ppoints}}, \code{\link{pp}}, \code{\link{rl}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{qq(M1)}
}
\keyword{hplot}