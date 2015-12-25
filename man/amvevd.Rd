\name{amvevd}
\alias{amvevd}
\title{Parametric Dependence Functions of Multivariate Extreme
    Value Models}
\description{
  Calculate the dependence function \eqn{A} for the multivariate
  logistic and multivariate asymmetric logistic models; plot the
  estimated function in the trivariate case.
}
\usage{
amvevd(x = rep(1/d,d), dep, asy, model = c("log", "alog"), d = 3, plot =
    FALSE, col = heat.colors(12), blty = 0, grid = if(blty) 150 else 50,
    lower = 1/3, ord = 1:3, lab = as.character(1:3), lcex = 1)
}
\arguments{
  \item{x}{A vector of length \code{d} or a matrix with \code{d}
    columns, in which case the dependence function is evaluated
    across the rows (ignored if plot is \code{TRUE}). The
    elements/rows of the vector/matrix should be positive and should
    sum to one, or else they should have a positive sum, in which
    case the rows are rescaled and a warning is given.
    \eqn{A(1/d,\dots,1/d)} is returned by default since it is often
    a useful summary of dependence.}
  \item{dep}{The dependence parameter(s). For the logistic model,
    should be a single value. For the asymmetric logistic model,
    should be a vector of length \eqn{2^d-d-1}, or a single value,
    in which case the value is used for each of the \eqn{2^d-d-1}
    parameters (see \code{\link{rmvevd}}).}
  \item{asy}{The asymmetry parameters for the asymmetric logistic
    model. Should be a list with \eqn{2^d-1} vector elements
    containing the asymmetry parameters for each separate component
    (see \code{\link{rmvevd}} and \bold{Examples}).}
  \item{model}{The specified model; a character string. Must be
    either \code{"log"} (the default) or \code{"alog"} (or any unique
    partial match), for the logistic and asymmetric logistic models
    respectively. The definition of each model is given in
    \code{\link{rmvevd}}.}
  \item{d}{The dimension; an integer greater than or equal to two.
    The trivariate case \code{d = 3} is the default.}
  \item{plot}{Logical; if \code{TRUE}, and the dimension \code{d} is
    three (the default dimension), the dependence function of a
    trivariate model is plotted. For plotting in the bivariate case,
    use \code{\link{abvevd}}. If \code{FALSE} (the default), the
    following arguments are ignored.}
  \item{col}{A list of colours (see \code{\link{image}}). The first
    colours in the list represent smaller values, and hence
    stronger dependence. Each colour represents an equally spaced
    interval between \code{lower} and one.}
  \item{blty}{The border line type, for the border that surrounds
    the triangular image. By default \code{blty} is zero, so no
    border is plotted. Plotting a border leads to (by default) an
    increase in \code{grid} (and hence computation time), to ensure
    that the image fits within it.}
  \item{grid}{For plotting, the function is evaluated at \code{grid^2}
    points.}
  \item{lower}{The minimum value for which colours are plotted. By
    defualt \eqn{\code{lower} = 1/3} as this is the theoretical
    minimum of the dependence function of the trivariate extreme
    value distribution.}
  \item{ord}{A vector of length three, which should be a permutation
    of the set \eqn{\{1,2,3\}}{{1,2,3}}. The points
    \eqn{(1,0,0)}, \eqn{(0,1,0)} and \eqn{(0,0,1)} (the vertices of
    the simplex) are depicted clockwise from the top in
    the order defined by \code{ord}.The argument alters the way in
    which the function is plotted; it does not change the function
    definition.}
  \item{lab}{A character vector of length three, in which case the
    \code{i}th margin is labelled using the \code{i}th component,
    or \code{NULL}, in which case no labels are given. The actual
    location of the margins, and hence the labels, is defined by
    \code{ord}.}
  \item{lcex}{A numerical value giving the amount by which the
    labels should be scaled relative to the default. Ignored
    if \code{lab} is \code{NULL}.}
}
\details{
  Let \eqn{z = (z_1,\dots,z_d)}{z = (z1,\dots,zd)} and
  \eqn{w = (w_1,\dots,w_d)}{w = (w1,\dots,wd)}.
  Any multivariate extreme value distribution can be written as
  \deqn{G(z) = \exp\left\{- \left\{\sum\nolimits_{j=1}^{d} y_j
    \right\} A\left(\frac{y_1}{\sum\nolimits_{j=1}^{d} y_j}, \dots,
    \frac{y_d}{\sum\nolimits_{j=1}^{d} y_j}\right)\right\}}{
    G(z) = exp{-(y1+\dots+yd) A[y1/(y1+\dots+yd),
      \dots, yd/(y1+\dots+yd)]}}
  for some function \eqn{A} defined on the simplex
  \eqn{S_d = \{w \in R^d_+ : \sum\nolimits_{j=1}^{d}
    w_j = 1\}}{S_d = {w: w1 + \dots + wd = 1}},
  where \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,\dots,d}, and where the (generalized extreme value)
  marginal parameters are given by
  \eqn{(a_i,b_i,s_i)}{(ai,bi,si)}, \eqn{b_i > 0}{bi > 0}.
  If \eqn{s_i = 0}{si = 0} then \eqn{y_i}{yi} is defined by
  continuity.

  \eqn{A} is called (by some authors) the dependence function.
  It follows that \eqn{A(w) = 1} when \eqn{w} is one of the \eqn{d}
  vertices of \eqn{S_d}, and that \eqn{A} is a convex function with
  \eqn{\max(w_1,\dots,w_d) \leq A(w)\leq 1}{
    max(w1,\dots,wd) <= A(w) <= 1} for
  all \eqn{w} in \eqn{S_d}.
  The lower and upper limits of \eqn{A} are obtained under complete
  dependence and mutual independence respectively.
  \eqn{A} does not depend on the marginal parameters.
}

\value{
  A numeric vector of values. If plotting, the smallest evaluated
  function value is returned invisibly.
}

\seealso{\code{\link{amvnonpar}}, \code{\link{abvevd}},
  \code{\link{rmvevd}}, \code{\link{image}}}

\examples{
amvevd(dep = 0.5, model = "log")
s3pts <- matrix(rexp(30), nrow = 10, ncol = 3)
s3pts <- s3pts/rowSums(s3pts)
amvevd(s3pts, dep = 0.5, model = "log")
\dontrun{amvevd(dep = 0.05, model = "log", plot = TRUE, blty = 1)}
amvevd(dep = 0.95, model = "log", plot = TRUE, lower = 0.94)

asy <- list(.4, .1, .6, c(.3,.2), c(.1,.1), c(.4,.1), c(.2,.3,.2))
amvevd(s3pts, dep = 0.15, asy = asy, model = "alog")
amvevd(dep = 0.15, asy = asy, model = "al", plot = TRUE, lower = 0.7)
}
\keyword{distribution}








