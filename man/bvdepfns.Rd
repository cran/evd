\name{bvdepfns}
\alias{abvalog}
\alias{abvaneglog}
\alias{abvhr}
\alias{abvlog}
\alias{abvneglog}
\alias{abvbilog}
\alias{abvnegbilog}
\alias{abvct}
\title{Dependence Functions for Bivariate Extreme Value Distributions}
\description{
  Calculate or plot the dependence function for the
  symmetric/asymmetric logistic, symmetric/asymmetric negative logistic,
  bilogistic, negative bilogistic, Husler-Reiss and Coles-Tawn models.
}
\usage{
abvalog(x = 0.5, dep, asy = c(1, 1), plot = FALSE, border = TRUE,
    add = FALSE, lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1),
    xlab = "", ylab = "", \dots) 
abvaneglog(x = 0.5, dep, asy = c(1, 1), plot = FALSE, border = TRUE,
    add = FALSE, lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1),
    xlab = "", ylab = "", \dots) 
abvhr(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, 
    lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1), xlab = "", 
    ylab = "", \dots) 
abvlog(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, 
    lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1), xlab = "", 
    ylab = "", \dots) 
abvneglog(x = 0.5, dep, plot = FALSE, border = TRUE, add = FALSE, 
    lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1), xlab = "", 
    ylab = "", \dots)
abvbilog(x = 0.5, alpha, beta, plot = FALSE, border = TRUE, add = FALSE,
    lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1), xlab = "",
    ylab = "", \dots) 
abvnegbilog(x = 0.5, alpha, beta, plot = FALSE, border = TRUE,
    add = FALSE, lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1),
    xlab = "", ylab = "", \dots)
abvct(x = 0.5, alpha, beta, plot = FALSE, border = TRUE, 
    add = FALSE, lty = 1, blty = 3, xlim = c(0, 1), ylim = c(0.5, 1),
    xlab = "", ylab = "", \dots)  
}
\arguments{
  \item{x}{A vector of values at which the dependence function is
    evaluated (ignored if plot is \code{TRUE}).}
  \item{dep}{Dependence parameter.}
  \item{asy}{A vector containing two asymmetry parameters.}
  \item{alpha, beta}{Alpha and beta parameters.}
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
  \code{abvct} gives the dependence function for the Coles-Tawn
  model.
  \code{abvbilog} and \code{abvnegbilog} give the dependence function
  for the bilogistic and negative bilogistic models respectively.
}

\seealso{\code{\link{abvnonpar}}, \code{\link{rbvalog}},
  \code{\link{rbvaneglog}}, \code{\link{rbvhr}}, \code{\link{rbvlog}},
  \code{\link{rbvneglog}}, \code{\link{rbvbilog}},
  \code{\link{rbvnegbilog}}, \code{\link{rbvct}}}

\examples{
abvhr(dep = 2.7)
abvalog(dep = .3, asy = c(.7,.9))
abvalog(seq(0,1,0.25), dep = .3, asy = c(.7,.9))
\dontrun{abvnegbilog(alpha = .3, beta = 1.2, plot = TRUE)}
}
\keyword{distribution}








