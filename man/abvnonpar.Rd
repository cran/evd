\name{abvnonpar}
\alias{abvnonpar}
\title{Non-parametric Estimates for Dependence Functions of the
    Bivariate Extreme Value Distribution}
\description{
  Calculate or plot non-parametric estimates for the dependence function
  \eqn{A} of the bivariate extreme value distribution.
}
\synopsis{
abvnonpar(x = 0.5, data, epmar = FALSE, nsloc1 = NULL, nsloc2 = NULL,
    method = c("cfg", "pickands", "deheuvels", "halltajvidi", "tdo"),
    convex = FALSE, rev = FALSE, madj = 0, kmar = NULL, plot = FALSE,
    add = FALSE, lty = 1, lwd = 1, col = 1, blty = 3, blwd = 1, xlim =
    c(0, 1), ylim = c(0.5, 1), xlab = "t", ylab = "A(t)", \dots)
}
\usage{
abvnonpar(x = 0.5, data, epmar = FALSE, nsloc1 = NULL, nsloc2 = NULL,
    method = c("cfg", "pickands"), convex = FALSE, rev = FALSE,
    madj = 0, plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1,
    blty = 3, blwd = 1, xlim = c(0, 1), ylim = c(0.5, 1), xlab = "t",
    ylab = "A(t)", \dots)
}
\arguments{
  \item{x}{A vector of values at which the dependence function is
    evaluated (ignored if plot or add is \code{TRUE}). \eqn{A(1/2)}
    is returned by default since it is often a useful summary of
    dependence.}
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
  \item{method}{The estimation method (see \bold{Details}). Either
    \code{"cfg"} (the default) or \code{"pickands"}.}
  \item{convex}{Logical; take the convex minorant?}
  \item{rev}{Logical; reverse the dependence function? This is
    equivalent to evaluating the function at \code{1-x}.}
  \item{madj}{Performs marginal adjustments for the \code{"pickands"}
    method (see \bold{Details}).}
  \item{plot}{Logical; if \code{TRUE} the function is plotted. The
    x and y values used to create the plot are returned invisibly.
    If \code{plot} and \code{add} are \code{FALSE} (the default),
    the arguments following \code{add} are ignored.}
  \item{add}{Logical; add to an existing plot? The existing plot
    should have been created using either \code{abvnonpar} or
    \code{\link{abvevd}}, the latter of which plots (or calculates)
    the dependence function for a number of parametric models.}
  \item{lty, blty}{Function and border line types. Set \code{blty}
    to zero to omit the border.}
  \item{lwd, blwd}{Function and border line widths.}
  \item{col}{Line colour.}
  \item{xlim, ylim}{x and y-axis limits.}
  \item{xlab, ylab}{x and y-axis labels.}
  \item{\dots}{Other high-level graphics parameters to be passed to
    \code{plot}.}
}
\details{
  The dependence function  \eqn{A(\cdot)}{A()} of the bivariate
  extreme value distribution is defined in \code{\link{abvevd}}.
  Non-parametric estimates are constructed as follows.
  Suppose \eqn{(z_{i1},z_{i2})} for \eqn{i=1,\ldots,n} are \eqn{n}
  bivariate observations that are passed using the \code{data}
  argument.
  If \code{epmar} is \code{FALSE} (the default), then
  the marginal parameters of the GEV margins are estimated
  (under the assumption of independence) and the data is
  transformed using
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

  Two different estimators of the dependence function can be
  implemented.
  They are defined (on \eqn{0 \leq w \leq 1}{0 <= w <= 1}) as
  follows.

  \code{method = "cfg"} (Caperaa, Fougeres and Genest, 1997)
  \deqn{\log(A_c(w)) = 1/n \left\{ sum_{i=1}^n \log(max[(1-w)y_{i1},
    wy_{i1}]) - (1-w)sum_{i=1}^n y_{i1} - w sum_{i=1}^n y_{i2}
    \right\}}{log(A_c(w)) = 1/n { sum_{i=1}^n \log (max[(1-w)y_{i1},
    wy_{i1}]) - (1-w)sum_{i=1}^n y_{i1} - w sum_{i=1}^n y_{i2} }}
  
  \code{method = "pickands"} (Pickands, 1981)
  \deqn{A_p(w) = n\left\{\sum_{i=1}^n \min\left(\frac{y_{i1}}{w},
    \frac{y_{i2}}{1-w}\right)\right\}^{-1}}{
    A_p(w) = n / {sum_{i=1}^n min[y_{i1}/w, y_{i2}/(1-w)]}}

  Two variations on the estimator \eqn{A_p(\cdot)}{A_p()} are
  also implemented. If the argument \code{madj = 1}, an adjustment
  given in Deheuvels (1991) is applied. If the argument
  \code{madj = 2}, an adjustment given in Hall and Tajvidi (2000)
  is applied. These are marginal adjustments; they are only
  useful when empirical marginal estimation is used.

  Let \eqn{A_n(\cdot)}{A_n()} be any estimator of \eqn{A(\cdot)}{A()}.
  None of the estimators satisfy
  \eqn{\max(w,1-w) \leq A_n(w) \leq 1}{max(w,1-w) <= A_n(w) <= 1}
  for all \eqn{0\leq w \leq1}{0 <= w <= 1}. An obvious modification is
  \deqn{A_n^{'}(w) = \min(1, \max\{A_n(w), w, 1-w\}).}{
    A'_n(w) = min(1, max{A_n(w), w, 1-w}).}
  This modification is always implemented.

  Convex estimators can be derived by taking the convex minorant,
  which can be achieved by setting \code{convex} to \code{TRUE}.
}
\note{
  I have been asked to point out that Hall and Tajvidi (2000)
  suggest putting a constrained smoothing spline on their
  modified Pickands estimator, but this is not done here.
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

  Pickands, J. (1981)
  Multivariate extreme value distributions.
  \emph{Proc. 43rd Sess. Int. Statist. Inst.}, \bold{49}, 859--878.

  Deheuvels, P. (1991)
  On the limiting behaviour of the Pickands estimator 
  for bivariate extreme-value distributions.
  \emph{Statist. Probab. Letters}, \bold{12}, 429--439.

  Hall, P. and Tajvidi, N. (2000)
  Distribution and dependence-function estimation for
  bivariate extreme-value distributions.
  \emph{Bernoulli}, \bold{6}, 835--844.
}

\seealso{\code{\link{abvevd}}, \code{\link{amvnonpar}},
  \code{\link{fgev}}}

\examples{
bvdata <- rbvevd(100, dep = 0.7, model = "log")
abvnonpar(seq(0, 1, length = 10), data = bvdata, convex = TRUE)
abvnonpar(data = bvdata, method = "p", plot = TRUE)

M1 <- fitted(fbvevd(bvdata, model = "log"))
abvevd(dep = M1["dep"], model = "log", plot = TRUE)
abvnonpar(data = bvdata, add = TRUE, lty = 2)
}
\keyword{nonparametric}








