\name{bvdepfns}
\alias{abvalog}
\alias{abvaneglog}
\alias{abvhr}
\alias{abvlog}
\alias{abvneglog}
\title{Dependence Functions for Bivariate Extreme Value Distributions}
\description{
  Calculate or plot the dependence function for logistic, negative
  logistic and Husler-Reiss models.
}
\usage{
abvalog(x = 0.5, dep, asy, plot = FALSE, border = TRUE, add = FALSE, 
    lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1), xlab = "", 
    ylab = "", \dots) 
abvaneglog(x = 0.5, dep, asy, plot = FALSE, border = TRUE, add = FALSE, 
    lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1), xlab = "", 
    ylab = "", \dots) 
abvhr(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, 
    lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1), xlab = "", 
    ylab = "", \dots) 
abvlog(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, 
    lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1), xlab = "", 
    ylab = "", \dots) 
abvneglog(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, 
    lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1), xlab = "", 
    ylab = "", \dots) 
}
\arguments{
  \item{x}{a vector values at which the dependence function is
    evaluated (ignored if plot is \code{TRUE}).}
  \item{dep}{dependence parameter.}
  \item{asy}{a vector containing two asymmetry parameters.}
  \item{plot}{logical; if \code{TRUE} the function is plotted and
    the values used to create the plot are returned invisibly.}
  \item{border}{logical; if \code{TRUE} a border representing the
    maximal domain is added to the plot.}
  \item{add}{logical; add to an existing plot?}
  \item{lty, blty}{function and border line types.}
  \item{xlim, ylim}{x and y-axis limits.}
  \item{xlab, ylab}{x and y-axis labels.}
  \item{\dots}{other high-level graphics parameters to be passed to
    \code{plot}.}
}
\details{
  Any bivariate extreme value distribution can be written as
  \deqn{G(z_1,z_2) = \exp\left[-(y_1+y_2)A\left(
    \frac{y_1}{y_1+y_2}\right)\right]}{
    G(z1,z2) = exp{-(y1+y2)A[y1/(y1+y2)]}}
  for some function \eqn{A(\cdot)}{A()} defined on \eqn{[0,1]}, where    
  \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,2}, with the marginal
  parameters given by \eqn{(a_i,b_i,s_i)}{(ai,bi,si)},
  \eqn{b_i > 0}{bi > 0}.

  \eqn{A(\cdot)}{A()} is called (by some authors) the dependence
  function.
  It follows that \eqn{A(0)=A(1)=1}, and that \eqn{A(\cdot)}{A()} is
  a convex function with \eqn{\max(x,1-x) \leq A(x)\leq 1}{max(x,1-x)
    <= A(x) <= 1} for all \eqn{0\leq x\leq1}{0 <= x <= 1}.
  \eqn{A(\cdot)}{A()} does not depend on the marginal parameters.
  \eqn{A(1/2)} is returned by default since it is often a useful
  summary of dependence.
}

\value{
  \code{abvlog} and \code{abvalog} give the dependence function for the
  logistic and asymmetric logistic models respectively.
  \code{abvneglog} and \code{abvaneglog} give the dependence function
  for the negative logistic and asymmetric negative logistic models
  respectively.
  \code{abvhr} gives the dependence function for the Husler-Reiss
  model.
}

\seealso{\code{\link{rbvalog}}, \code{\link{rbvaneglog}},
  \code{\link{rbvhr}}, \code{\link{rbvlog}}, \code{\link{rbvneglog}}}

\examples{
abvhr(dep = 2.7)
abvalog(dep = .3, asy = c(.7,.9))
abvalog(seq(0,1,0.25), dep = .3, asy = c(.7,.9))
}
\keyword{distribution}








