\name{abvpar}
\alias{abvpar}
\title{Parametric Dependence Functions of Bivariate Extreme
    Value Models}
\description{
  Calculate or plot the dependence function \eqn{A} for
  eight parametric bivariate extreme value models.
}
\synopsis{
abvpar(x = 0.5, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct"),
     plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1, blty = 3,
     xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", \dots)
}
\usage{
abvpar(x = 0.5, dep, asy = c(1,1), alpha, beta, model = "log",
     plot = FALSE, add = FALSE, lty = 1, lwd = 1, col = 1, blty = 3,
     xlim = c(0,1), ylim = c(0.5,1), xlab = "", ylab = "", \dots)
}
\arguments{
  \item{x}{A vector of values at which the dependence function is
    evaluated (ignored if plot or add is \code{TRUE}). \eqn{A(1/2)}
    is returned by default since it is often a useful summary of
    dependence.}
  \item{dep}{Dependence parameter for the logistic, asymmetric
    logistic, Husler-Reiss, negative logistic and asymmetric
    negative logistic models.}
  \item{asy}{A vector of length two, containing the two asymmetry
    parameters for the asymmetric logistic and asymmetric negative
    logistic models.}
  \item{alpha, beta}{Alpha and beta parameters for the bilogistic,
    negative bilogistic and Coles-Tawn models.}
  \item{model}{The specified model; a character string. Must be
    either \code{"log"} (the default), \code{"alog"}, \code{"hr"},
    \code{"neglog"}, \code{"aneglog"}, \code{"bilog"},
    \code{"negbilog"} or \code{"ct"} (or any unique partial match),
    for the logistic, asymmetric logistic, Husler-Reiss, negative
    logistic, asymmetric negative logistic, bilogistic, negative
    bilogistic and Coles-Tawn models respectively. The definition
    of each model is given in \code{\link{rbvevd}}. If parameter
    arguments are given that do not correspond to the specified
    model those arguments are ignored, with a warning.}
  \item{plot}{Logical; if \code{TRUE} the function is plotted. The
    x and y values used to create the plot are returned invisibly.
    If \code{plot} and \code{add} are \code{FALSE} (the default),
    the arguments following \code{add} are ignored.}
  \item{add}{Logical; add to an existing plot? The existing plot
    should have been created using either \code{abvpar} or
    \code{\link{abvnonpar}}, the latter of which plots (or calculates)
    a non-parametric estimate of the dependence function.}
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
  Any bivariate extreme value distribution can be written as
  \deqn{G(z_1,z_2) = \exp\left[-(y_1+y_2)A\left(
    \frac{y_1}{y_1+y_2}\right)\right]}{
    G(z1,z2) = exp{-(y1+y2)A[y1/(y1+y2)]}}
  for some function \eqn{A(\cdot)}{A()} defined on \eqn{[0,1]}, where    
  \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,2}, with the (generalized extreme value) marginal
  parameters given by \eqn{(a_i,b_i,s_i)}{(ai,bi,si)},
  \eqn{b_i > 0}{bi > 0}.
  If \eqn{s_i = 0}{si = 0} then \eqn{y_i}{yi} is defined by
  continuity.

  \eqn{A(\cdot)}{A()} is called (by some authors) the dependence
  function.
  It follows that \eqn{A(0)=A(1)=1}, and that \eqn{A(\cdot)}{A()} is
  a convex function with \eqn{\max(x,1-x) \leq A(x)\leq 1}{max(x,1-x)
    <= A(x) <= 1} for all \eqn{0\leq x\leq1}{0 <= x <= 1}.
  The lower and upper limits of \eqn{A} are obtained under complete
  dependence and independence respectively.
  \eqn{A(\cdot)}{A()} does not depend on the marginal parameters.
}

\value{
  \code{abvpar} calculates or plots the dependence function
  for one of eight parametric bivariate extreme value models,
  at specified parameter values.
}

\seealso{\code{\link{abvnonpar}}, \code{\link{fbvevd}},
  \code{\link{rbvevd}}, \code{\link{atvpar}}}

\examples{
abvpar(dep = 2.7, model = "hr")
abvpar(seq(0,1,0.25), dep = 0.3, asy = c(.7,.9), model = "alog")
\dontrun{abvpar(alpha = 0.3, beta = 1.2, model = "negbi", plot = TRUE)}

bvdata <- rbvevd(100, dep = 0.7, model = "log")
M1 <- fitted(fbvevd(bvdata, model = "log"))
abvpar(dep = M1["dep"], model = "log", plot = TRUE)
abvnonpar(data = bvdata, add = TRUE, lty = 2)
}
\keyword{distribution}








