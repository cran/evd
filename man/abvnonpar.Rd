\name{abvnonpar}
\alias{abvnonpar}
\title{Non-parametric Estimates for the Dependence Function}
\description{
  Calculate or plot non-parametric estimates for the dependence function
  of the bivariate extreme value distribution.
}
\usage{
abvnonpar(x = 0.5, data, nsloc1 = NULL, nsloc2 = NULL,
    method = c("cfg", "deheuvels", "pickands"), modify = 0,
    wf = function(t) t, plot = FALSE, border = TRUE, add = FALSE,
    lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1),
    xlab = "", ylab = "", \dots) 
}
\arguments{
  \item{x}{A vector of values at which the dependence function is
    evaluated (ignored if plot is \code{TRUE}).}
  \item{data}{A matrix or data frame with two columns, which may contain
    missing values.}
  \item{nsloc1,nsloc2}{A data frame with the same number of rows as the
    length of \code{x}, for linear modelling of the location
    parameter on the first/second margin (see \bold{Details}).
    The data frames are treated as covariate matrices (excluding the
    intercept).
    A numeric vector can be given as an alternative to a single column
    data frame.}
  \item{method}{The estimation method (see \bold{Details}).}
  \item{modify}{An integer. Can be used to modify the estimation method
    (see \bold{Details}).}
  \item{wf}{The weight function used in the (default) ``cfg''
    method (see \bold{Details}). The function must be vectorized.}
  \item{plot}{Logical; if \code{TRUE} the function is plotted and
    the values used to create the plot are returned invisibly.}
  \item{border}{Logical; if \code{TRUE} a border representing the
    maximal domain is added to the plot.}
  \item{add}{Logical; add to an existing plot?}
  \item{lty, blty}{Function and border line types.}
  \item{xlim, ylim}{x and y-axis limits.}
  \item{xlab, ylab}{x and y-axis labels.}
  \item{\dots}{Other high-level graphics parameters to be passed to
    \code{plot}.}
}
\details{
  Any bivariate extreme value distribution can be written as
  \deqn{G(z_1,z_2) = \exp\left[-(y_1+y_2)A\left(
    \frac{y_1}{y_1+y_2}\right)\right]}{
    G(z_1,z_2) = exp{-(y_1 + y_2) A[y_1/(y_1+y_2)]}}
  for some function \eqn{A(\cdot)}{A()} defined on \eqn{[0,1]}, where    
  \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    y_i = {1 + s_i(z_i-a_i)/b_i}^(-1/s_i)}
  for \eqn{1 + s_i(z_i-a_i)/b_i > 0} and
  \eqn{i = 1,2}, with the marginal
  parameters given by \eqn{(a_i,b_i,s_i)},
  \eqn{b_i > 0}.

  \eqn{A(\cdot)}{A()} is called (by some authors) the dependence
  function.
  It follows that \eqn{A(0) = A(1) = 1},
  and that \eqn{A(\cdot)}{A()} is
  a convex function with \eqn{\max(w,1-w) \leq A(w)\leq 1}{max(w,1-w)
    <= A(w) <= 1} for all \eqn{0\leq w\leq1}{0 <= w <= 1}.
  \eqn{A(\cdot)}{A()} does not depend on the marginal parameters.

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
  If non-stationary fitting is implemented using the \code{nsloc1}
  or \code{nsloc2} arguments (see \code{\link{fgev}} or
  \code{\link{fbvlog}}) the marginal location parameters may depend
  on \eqn{i}.

  Three different estimators of the dependence function can be
  implemented.
  They are defined (on \eqn{0 \leq w \leq 1}{0 <= w <= 1}) as
  follows.

  \eqn{\code{method} = ``pickands''} (Pickands, 1981)
  \deqn{A_p(w) = n\left\{\sum_{i=1}^n \min\left(\frac{y_{i1}}{w},
    \frac{y_{i2}}{1-w}\right)\right\}^{-1}}{
    A_p(w) = n / {sum_{i=1}^n min[y_{i1}/w, y_{i2}/(1-w)]}}

  \eqn{\code{method} = ``deheuvels''} (Deheuvels, 1991)
  \deqn{A_d(w) = n\left\{\sum_{i=1}^n \min\left(\frac{y_{i1}}{w},
    \frac{y_{i2}}{1-w}\right) - w\sum_{i=1}^n y_{i1} - (1-w)
    \sum_{i=1}^n y_{i2} + n\right\}^{-1}}{
        A_d(w) = n / {sum_{i=1}^n min[y_{i1}/w, y_{i2}(1-w)] -
	w sum_{i=1}^n y_{i1} - (1-w) sum_{i=1}^n y_{i2} + n}}

  \eqn{\code{method} = ``cfg''}; The Default Method
  (Caperaa, Fougeres and Genest, 1997)
  \deqn{A_c(w) = \exp\left\{ \{1-p(t)\} \int_{0}^{t}
    \frac{H(x) - x}{x(1-x)} \, \mbox{d}x - p(t) \int_{t}^{1}
    \frac{H(x) - x}{x(1-x)} \, \mbox{d}x \right\}}{
    A_c(w) = exp{ [1-p(t)] integral_0^t (H(x) - x)/[x(1-x)] dx -
      p(t) integral_t^1 (H(x) - x)/[x(1-x)] dx }}

  In the estimator \eqn{A_c(\cdot)}{A_c()}, \eqn{H(x)} is the
  empirical distribution function of \eqn{x_1,\ldots,x_n}, where
  \eqn{x_i = y_{i2} / (y_{i1} + y_{i2})} for \eqn{i = 1,\ldots,n},
  and \eqn{p(t)} is any bounded function on \eqn{[0,1]}, which
  can be specified using the argument \code{wf}.
  By default \code{wf} is the identity function.

  Let \eqn{A_n(\cdot)}{A_n()} be any estimator of \eqn{A(\cdot)}{A()}.
  The estimator \eqn{A_d(\cdot)}{A_d()} satisfies
  \eqn{A_n(0) = A_n(1) = 1}.
  \eqn{A_c(\cdot)}{A_c()} satisfies this constraint when
  \eqn{p(0) = 0} and \eqn{p(1) = 1}.

  None of the estimators satisfy
  \eqn{\max(w,1-w) \leq A_n(w) \leq 1}{max(w,1-w) <= A_n(w) <= 1}
  for all \eqn{0\leq w \leq1}{0 <= w <= 1}. An obvious modification is
  \deqn{A_n^{'}(w) = \min(1, \max\{A_n(w), w, 1-w\}).}{
        A'_n(w) = min(1, max{A_n(w), w, 1-w}).}

  Another estimator \eqn{A_n^{''}(w)}{A''_n(w)} can be derived by
  taking the convex hull of \eqn{A_n^{'}(w)}{A'_n(w)}.
  These modifications can be implemented using the \code{modify}
  argument.
  Set \eqn{\code{modify} = 1} to plot or calculate
  \eqn{A_n^{'}(w)}{A'_n(w)}.
  Set \eqn{\code{modify} = 2} to plot or calculate
  \eqn{A_n^{''}(w)}{A''_n(w)}. 

  \eqn{A_n(1/2)} is returned by default since it is often a useful
  summary of dependence.
}

\value{
  \code{abvnonpar} gives a non-parametric estimate of the dependence
  function.
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

  Pickands, J. (1981)
  Multivariate extreme value distributions.
  \emph{Proc. 43rd Sess. Int. Statist. Inst.}, \bold{49}, 859--878.
}

\seealso{\code{\link{abvlog}}, \code{\link{fbvlog}}, \code{\link{fgev}}}

\examples{
bvdata <- rbvlog(100, dep = 0.7)
abvnonpar(seq(0, 1, length = 10), data = bvdata, modify = 2)
\dontrun{abvnonpar(data = bvdata, method = "d", plot = TRUE)}
}
\keyword{distribution}








