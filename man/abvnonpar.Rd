\name{abvnonpar}
\alias{abvnonpar}
\title{Non-parametric Estimates for Dependence Functions of the
    Bivariate Extreme Value Distribution}
\description{
  Calculate or plot non-parametric estimates for the dependence function
  \eqn{A} of the bivariate extreme value distribution.
}
\synopsis{
abvnonpar(x = 0.5, data, nsloc1 = NULL, nsloc2 = NULL,
    method = c("cfg", "pickands", "deheuvels", "halltajvidi", "tdo"),
    convex = FALSE, wf = function(t) t, kmar = NULL, plot = FALSE,
    add = FALSE, lty = 1, lwd = 1, col = 1, blty = 3, xlim = c(0, 1),
    ylim = c(0.5, 1), xlab = "", ylab = "", \dots)
}
\usage{
abvnonpar(x = 0.5, data, nsloc1 = NULL, nsloc2 = NULL,
    method = c("cfg", "pickands", "deheuvels", "halltajvidi", "tdo"),
    convex = FALSE, wf = function(t) t, plot = FALSE,
    add = FALSE, lty = 1, lwd = 1, col = 1, blty = 3, xlim = c(0, 1),
    ylim = c(0.5, 1), xlab = "", ylab = "", \dots)
}
\arguments{
  \item{x}{A vector of values at which the dependence function is
    evaluated (ignored if plot or add is \code{TRUE}). \eqn{A(1/2)}
    is returned by default since it is often a useful summary of
    dependence.}
  \item{data}{A matrix or data frame with two columns, which may contain
    missing values.}
  \item{nsloc1, nsloc2}{A data frame with the same number of rows as
    \code{data}, for linear modelling of the location
    parameter on the first/second margin.
    The data frames are treated as covariate matrices, excluding the
    intercept.
    A numeric vector can be given as an alternative to a single column
    data frame.}
  \item{method}{The estimation method (see \bold{Details}). The
    \code{"cfg"} method is used by default.}
  \item{convex}{Logical; take the convex minorant?}
  \item{wf}{The weight function used in the \code{"cfg"}
    method (see \bold{Details}). The function must be vectorized.}
  \item{plot}{Logical; if \code{TRUE} the function is plotted. The
    x and y values used to create the plot are returned invisibly.
    If \code{plot} and \code{add} are \code{FALSE} (the default),
    the arguments following \code{add} are ignored.}
  \item{add}{Logical; add to an existing plot? The existing plot
    should have been created using either \code{abvnonpar} or
    \code{\link{abvpar}}, the latter of which plots (or calculates)
    the dependence function for a number of parametric models.}
  \item{lty, blty}{Function and border line types. Set \code{blty}
    to zero to omit the border.}
  \item{lwd}{Line width.}
  \item{col}{Line colour.}
  \item{xlim, ylim}{x and y-axis limits.}
  \item{xlab, ylab}{x and y-axis labels.}
  \item{\dots}{Other high-level graphics parameters to be passed to
    \code{plot}.}
}
\details{
  The dependence function  \eqn{A(\cdot)}{A()} of the bivariate
  extreme value distribution is defined in \code{\link{abvpar}}.
  Non-parametric estimates are constructed as follows.
  Suppose \eqn{(z_{i1},z_{i2})} for \eqn{i=1,\ldots,n} are \eqn{n}
  bivariate observations that are passed using the \code{data}
  argument.
  The marginal parameters are estimated (under the assumption of
  independence) and the data is transformed using
  \deqn{y_{i1} = \{1+\hat{s}_1(z_{i1}-\hat{a}_1)/
    \hat{b}_1\}_{+}^{-1/\hat{s}_1}}{
        y_{i1} = {1 + s'_1(z_{i1}-a'_1)/b'_1}^(-1/s'_1)}
  and
  \deqn{y_{i2} = \{1+\hat{s}_2(z_{i2}-\hat{a}_2)/
    \hat{b}_2\}_{+}^{-1/\hat{s}_2}}{
        y_{i2} = {1 + s'_2(z_{i2}-a'_2)/b'_2}^(-1/s'_2)}
  for \eqn{i = 1,\ldots,n}, where
  \eqn{(\hat{a}_1,\hat{b}_1,\hat{s}_1)}{(a'_1,b'_1,s'_1)} and
  \eqn{(\hat{a}_2,\hat{b}_2,\hat{s}_2)}{(a'_2,b'_2,s'_2)}
  are the maximum likelihood estimates for the location, scale
  and shape parameters on the first and second margins.
  If \code{nsloc1} or \code{nsloc2} are given, the location
  parameters may depend on \eqn{i} (see \code{\link{fgev}}).

  Five different estimators of the dependence function can be
  implemented.
  They are defined (on \eqn{0 \leq w \leq 1}{0 <= w <= 1}) as
  follows.

  \code{method = "cfg"} (Caperaa, Fougeres and Genest, 1997)
  \deqn{A_c(w) = \exp\left\{ \{1-p(w)\} \int_{0}^{w}
    \frac{H(x) - x}{x(1-x)} \, \mbox{d}x - p(w) \int_{w}^{1}
    \frac{H(x) - x}{x(1-x)} \, \mbox{d}x \right\}}{
    A_c(w) = exp{ [1-p(w)] integral_0^w (H(x) - x)/[x(1-x)] dx -
      p(w) integral_w^1 (H(x) - x)/[x(1-x)] dx }}
  
  \code{method = "pickands"} (Pickands, 1981)
  \deqn{A_p(w) = n\left\{\sum_{i=1}^n \min\left(\frac{y_{i1}}{w},
    \frac{y_{i2}}{1-w}\right)\right\}^{-1}}{
    A_p(w) = n / {sum_{i=1}^n min[y_{i1}/w, y_{i2}/(1-w)]}}

  \code{method = "deheuvels"} (Deheuvels, 1991)
  \deqn{A_d(w) = n\left\{\sum_{i=1}^n \min\left(\frac{y_{i1}}{w},
    \frac{y_{i2}}{1-w}\right) - w\sum_{i=1}^n y_{i1} - (1-w)
    \sum_{i=1}^n y_{i2} + n\right\}^{-1}}{
        A_d(w) = n / {sum_{i=1}^n min[y_{i1}/w, y_{i2}(1-w)] -
	w sum_{i=1}^n y_{i1} - (1-w) sum_{i=1}^n y_{i2} + n}}

  \code{method = "halltajvidi"} (Hall and Tajvidi, 2000)
  \deqn{A_h(w) = n\left\{\sum_{i=1}^n \min\left(\frac{y_{i1}}{\bar{y}_1
      w},\frac{y_{i2}}{\bar{y}_2 (1-w)}\right)\right\}^{-1}}{
    A_h(w) = n (sum_{i=1}^n min[y_{i1}/(by_1 w),
    y_{i2}/(by_2 (1-w))])^{-1}}

  \code{method = "tdo"} (Tiago de Oliveira, 1997)
  \deqn{A_t(w) = 1 - \frac{1}{1 + \log n} \sum_{i=1}^n
    \min\left(\frac{w}{1+ny_{i1}},\frac{1-w}{1+ny_{i2}}\right)}{
    A_t(w) = 1 - 1/(1 + log n) sum_{i=1}^n
    min[w/(1 + n y_{i1}), (1 - w)/(1 + n y_{i2})]}
    
  In the estimator \eqn{A_h(\cdot)}{A_h()},
  \eqn{\bar{y}_j = n^{-1}\sum_{i=1}^n y_{ij}}{
    by_j = (sum_{i=1}^n y_{ij})/n} for \eqn{j = 1,2}.
  In the estimator \eqn{A_c(\cdot)}{A_c()}, \eqn{H(x)} is the
  empirical distribution function of \eqn{x_1,\ldots,x_n}, where
  \eqn{x_i = y_{i1} / (y_{i1} + y_{i2})} for \eqn{i = 1,\ldots,n},
  and \eqn{p(w)} is any bounded function on \eqn{[0,1]}, which
  can be specified using the argument \code{wf}.
  By default \code{wf} is the identity function.

  Let \eqn{A_n(\cdot)}{A_n()} be any estimator of \eqn{A(\cdot)}{A()}.
  The constraint \eqn{A_n(0) = A_n(1) = 1} is satisfied by
  \eqn{A_d(\cdot)}{A_d()}, \eqn{A_t(\cdot)}{A_t()} and
  \eqn{A_h(\cdot)}{A_h()}, and by \eqn{A_c(\cdot)}{A_c()} when
  \eqn{p(0) = 0} and \eqn{p(1) = 1}.
  None of the estimators satisfy
  \eqn{\max(w,1-w) \leq A_n(w) \leq 1}{max(w,1-w) <= A_n(w) <= 1}
  for all \eqn{0\leq w \leq1}{0 <= w <= 1}. An obvious modification is
  \deqn{A_n^{'}(w) = \min(1, \max\{A_n(w), w, 1-w\}).}{
    A'_n(w) = min(1, max{A_n(w), w, 1-w}).}
  This modification is always implemented.

  \eqn{A_t(w)}{A_t(w)} is the only estimator that is convex.
  Convex estimators can be derived from other methods by taking
  the convex minorant, which can be achieved by setting \code{convex}
  to \code{TRUE}.
}
\note{
  Appendix A of the User's Guide contains a short simulation
  study that compares the estimators defined above.

  The estimators \eqn{A_p(\cdot)}{A_p()}, \eqn{A_d(\cdot)}{A_d()}
  and \eqn{A_h(\cdot)}{A_h()} are very similar, and may not be
  distinguishable when plotted.

  I have been asked to point out that Hall, P. and Tajvidi, N.
  (2000) suggest putting a constrained smoothing spline on their
  estimator, but this is not done here.
}

\value{
  \code{abvnonpar} calculates or plots a non-parametric estimate of
  the dependence function of the bivariate extreme value distribution.
}
\references{
  Caperaa, P. Fougeres, A.-L. and Genest, C. (1997)
  A non-parametric estimation procedure for bivariate
  extreme value copulas.
  \emph{Biometrika}, \bold{84}, 567--577.

  Deheuvels, P. (1991)
  On the limiting behaviour of the Pickands estimator 
  for bivariate extreme-value distributions.
  \emph{Statist. Probab. Letters}, \bold{12}, 429--439.

  Hall, P. and Tajvidi, N. (2000)
  Distribution and dependence-function estimation for
  bivariate extreme-value distributions.
  \emph{Bernoulli}, \bold{6}, 835--844.

  Pickands, J. (1981)
  Multivariate extreme value distributions.
  \emph{Proc. 43rd Sess. Int. Statist. Inst.}, \bold{49}, 859--878.

  Tiago de Oliveira, J. (1997)
  \emph{Statistical Analysis of Extremes.} Pendor.
}

\seealso{\code{\link{abvpar}}, \code{\link{atvnonpar}},
  \code{\link{fgev}}}

\examples{
bvdata <- rbvevd(100, dep = 0.7, model = "log")
abvnonpar(seq(0, 1, length = 10), data = bvdata, convex = TRUE)
abvnonpar(data = bvdata, method = "d", plot = TRUE)

M1 <- fitted(fbvevd(bvdata, model = "log"))
abvpar(dep = M1["dep"], model = "log", plot = TRUE)
abvnonpar(data = bvdata, add = TRUE, lty = 2)
}
\keyword{nonparametric}








