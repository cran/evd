\name{bvdp}
\alias{bvdp}

\title{A Dependence Function Plot for a Bivariate evd Object}
\description{
  The dependence function for the fitted model is plotted and
  (optionally) compared to a non-parameteric estimate.
}
\usage{
bvdp(x, method = "cfg", modify = 0, wf = function(t) t, add = FALSE,
    lty = 1, nplty = 2, blty = 3, main = "Dependence Function", 
    xlab = "", ylab = "", \dots)
}
\arguments{
  \item{x}{An object of class \code{"bvevd"}.}
  \item{method,modify,wf}{Arguments passed to \code{\link{abvnonpar}},
      which calculates and plots non-parametric dependence function
      estimates.}
  \item{add}{Logical; add to an existing plot?}
  \item{lty,nplty,blty}{Line types; for the model estmate, the
    non-parametric estimate and the border respectively.
    Use zero to suppress.}  
  \item{main}{Title of plot.}
  \item{xlab,ylab}{Labels for x and y axes.}
  \item{\dots}{Other plot parameters.}
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
  \eqn{i = 1,2}, with the (generalized extreme value) marginal
  parameters given by \eqn{(a_i,b_i,s_i)},
  \eqn{b_i > 0}.

  \eqn{A(\cdot)}{A()} is called (by some authors) the dependence
  function.
  It follows that \eqn{A(0) = A(1) = 1},
  and that \eqn{A(\cdot)}{A()} is
  a convex function with \eqn{\max(w,1-w) \leq A(w)\leq 1}{max(w,1-w)
    <= A(w) <= 1} for all \eqn{0\leq w\leq1}{0 <= w <= 1}.
  \eqn{A(\cdot)}{A()} does not depend on the marginal parameters.
  
  For non-stationary models the data are transformed to stationarity.
  The plot then corresponds to the distribution obtained when all
  covariates are zero.
}

\seealso{\code{\link{abvnonpar}}, \code{\link{abvlog}},
  \code{\link{bvcpp}}, \code{\link{bvdens}}, \code{\link{plot.bvevd}}}

\examples{
bvdata <- rbvlog(100, dep = 0.6)
M1 <- fbvlog(bvdata)
\dontrun{bvdp(M1)}
}
\keyword{hplot}
