\name{bvcpp}
\alias{bvcpp}

\title{A Conditional P-P Plot for a Bivariate evd Object}
\description{
  A P-P plot for the condition distribution function of a
  bivariate evd object.
}
\usage{
bvcpp(x, mar = 2, ci = TRUE, main = "Conditional Probability Plot", 
    xlab = "Empirical", ylab = "Model", \dots)
}
\arguments{
  \item{x}{An object of class \code{"bvevd"}.}
  \item{mar}{The margin that is conditioned on; one (first) or two
    (second, the default).}
    \item{ci}{Logical; if \code{TRUE} (the default), plot simulated
    95\% confidence intervals.}
  \item{main}{Title of plot.}
  \item{xlab,ylab}{Labels for x and y axes.}
  \item{\dots}{Other plot parameters.}
}

\details{
  Let \eqn{G(.|.)} be the conditional distribution of
  the first margin given the second, under the fitted model.
  Let \eqn{z_1,\ldots,z_m} be the data used in the fitted model,
  where \eqn{z_j = (z_{1j}, z_{2j})} for \eqn{j = 1,\ldots,m}.
  
  The conditional P-P plot with \eqn{\code{mar} = 2} constists
  of the points
  \deqn{\{(p_i, c_i), i = 1,\ldots,m\}}{
    {(p_i, c_i), i = 1,\ldots,m}}
  where \eqn{p_1,\ldots,p_m} are plotting points defined by
  \code{\link{ppoints}} and \eqn{c_i} is the \eqn{i}th largest
  value from the sample
  \eqn{\{G(z_{j1}|z_{j2}), j = 1,\ldots,m\}.}{
    {G(z_{j1}|z_{j2}), j = 1,\ldots,m}.}
  When \eqn{\code{mar} = 1} the margins are reversed, so that
  \eqn{G(.|.)} is the conditional distribution of the second margin
  given the first.
  
  For non-stationary models the data are transformed to stationarity.
  The plot then corresponds to the distribution obtained when all
  covariates are zero.
}

\seealso{\code{\link{bvdens}}, \code{\link{bvdp}},
  \code{\link{plot.bvevd}}}

\examples{
bvdata <- rbvlog(100, dep = 0.6)
M1 <- fbvlog(bvdata)
\dontrun{bvcpp(M1)}
}
\keyword{hplot}