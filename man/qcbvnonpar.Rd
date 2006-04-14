\name{qcbvnonpar}
\alias{qcbvnonpar}
\title{Non-parametric Estimates for Bivariate Quantile Curves}
\description{
  Calculate or plot non-parametric estimates for quantile curves
  of bivariate extreme value distributions.
}
\synopsis{
qcbvnonpar(p = seq(0.75, 0.95, 0.05), data, epmar = FALSE, nsloc1 =
    NULL, nsloc2 = NULL, mint = 1, method = c("cfg", "pickands",
    "deheuvels", "halltajvidi", "tdo"), convex = FALSE, madj = 0,
    kmar = NULL, plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1,
    xlim = range(data[,1], na.rm = TRUE), ylim = range(data[,2],
    na.rm = TRUE), xlab = colnames(data)[1], ylab = colnames(data)[2],
    \dots)
}
\usage{
qcbvnonpar(p = seq(0.75, 0.95, 0.05), data, epmar = FALSE, nsloc1 =
    NULL, nsloc2 = NULL, mint = 1, method = c("cfg", "pickands"),
    convex = FALSE, madj = 0, plot = FALSE, add = FALSE, lty = 1,
    lwd = 1, col = 1, xlim, ylim, xlab, ylab, \dots)
}
\arguments{
  \item{p}{A vector of lower tail probabilities. One quantile curve is
    calculated or plotted for each probability.}
  \item{data}{A matrix or data frame with two columns, which may contain
    missing values.}
  \item{epmar}{If \code{TRUE}, an empirical transformation of the
    marginals is performed in preference to marginal parametric
    GEV estimation, and the \code{nsloc} arguments are ignored.}
  \item{nsloc1, nsloc2}{A data frame with the same number of rows as
    \code{data}, for linear modelling of the location parameter on the
    first/second margin. The data frames are treated as covariate
    matrices, excluding the intercept. A numeric vector can be given
    as an alternative to a single column data frame.}
  \item{mint}{An integer \eqn{m}. Quantile curves are plotted or
    calculated using the lower tail probabilities \eqn{p^m}.}
  \item{method}{Arguments for the non-parametric estimate of the
    dependence function. See \code{\link{abvnonpar}}.}
  \item{convex, madj}{Other arguments for the non-parametric
    estimate of the dependence function.}
  \item{plot}{Logical; if \code{TRUE} the data is plotted along
    with the quantile curves. If \code{plot} and \code{add} are
    \code{FALSE} (the default), the arguments following \code{add}
    are ignored.}
  \item{add}{Logical; add quantile curves to an existing data plot?
    The existing plot should have been created using either
    \code{qcbvnonpar} or \code{\link{plot.bvevd}}, the latter of
    which can plot quantile curves for parametric fits.}
  \item{lty, lwd}{Line types and widths.}
  \item{col}{Line colour.}
  \item{xlim, ylim}{x and y-axis limits.}
  \item{xlab, ylab}{x and y-axis labels.}
  \item{\dots}{Other high-level graphics parameters to be passed to
    \code{plot}.}
}
\details{
  Let G be a fitted bivariate distribution function with
  margins \eqn{G_1} and \eqn{G_2}. A quantile curve for a fitted
  distribution function G at lower tail probability p is defined
  by
  \deqn{Q(G, p) = \{(y_1,y_1):G(y_1,y_2) = p\}.}{
    Q(G, p) = {(y_1,y_1):G(y_1,y_2) = p}.}
  
  For bivariate extreme value distributions, it consists
  of the points
  \deqn{\left\{G_1^{-1}(p_1),G_2^{-1}(p_2))\right\}}{
    {G_1^{-1}(p_1),G_2^{-1}(p_2))}}
  where \eqn{p_1 = p^{t/A(t)}} and \eqn{p_2 = p^{(1-t)/A(t)}},
  with \eqn{A} being the estimated dependence function defined
  in \code{\link{abvevd}}, and where \eqn{t} lies in the interval
  \eqn{[0,1]}.

  By default the margins \eqn{G_1} and \eqn{G_2} are modelled using
  estimated generalized extreme value distributions.
  For non-stationary generalized extreme value margins the plotted
  data are transformed to stationarity, and the plot corresponds
  to the distribution obtained when all covariates are zero.

  If \code{epmar} is \code{TRUE}, empirical transformations
  are used in preference to generalized extreme value models.
  Note that the marginal empirical quantile functions are
  evaluated using \code{\link{quantile}}, which linearly
  interpolates between data points, hence the curve will not
  be a step function.

  The idea behind the argument \eqn{\code{mint} = m} is that if
  G is fitted to a dataset of componentwise maxima, and the
  underlying observations are \emph{iid} distributed according
  to F, then if \eqn{m} is the size of the blocks over which the
  maxima were taken, approximately \eqn{F^m = G}, leading
  to \eqn{Q(F, p) = Q(G, p^m)}.
}

\value{
  \code{qcbvnonpar} calculates or plots non-parametric quantile
  curve estimates for bivariate extreme value distributions.
  If \code{p} has length one it returns a two column matrix
  giving points on the curve, else it returns a list of
  such matrices.
}

\seealso{\code{\link{abvevd}}, \code{\link{abvnonpar}},
  \code{\link{plot.bvevd}}}

\examples{
bvdata <- rbvevd(100, dep = 0.7, model = "log")
qcbvnonpar(c(0.9,0.95), data = bvdata, plot = TRUE)
qcbvnonpar(c(0.9,0.95), data = bvdata, epmar = TRUE, plot = TRUE)
}
\keyword{nonparametric}








