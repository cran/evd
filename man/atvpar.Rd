\name{atvpar}
\alias{atvpar}
\title{Parametric Dependence Functions of Trivariate Extreme
    Value Models}
\description{
  Calculate or plot the dependence function \eqn{A} for the trivariate
  logistic and trivariate asymmetric logistic models.
}
\usage{
atvpar(x = rep(1/3,3), dep, asy, model = c("log", "alog"), plot =
    FALSE, col = heat.colors(12), blty = 0, grid = if(blty) 150 else 50,
    lower = 1/3, ord = 1:3, lab = as.character(1:3), lcex = 1)
}
\arguments{
  \item{x}{A vector of length three or a matrix with three columns,
    in which case the dependence function is evaluated across
    the rows (ignored if plot is \code{TRUE}). The elements/rows
    of the vector/matrix should be positive and should sum to one,
    or else they should have a positive sum, in which case the rows
    are rescaled and a warning is given. \eqn{A(1/3,1/3,1/3)} is
    returned by default since it is often a useful summary of
    dependence.}
  \item{dep}{The dependence parameter(s). For the logistic model,
    should be a single value. For the asymmetric logistic model,
    should be a vector of length four, or a single value,
    in which case the value is used for each of the four
    parameters (see \code{\link{rmvevd}}).}
  \item{asy}{The asymmetry parameters for the asymmetric logistic
    model. Should be a list with seven vector elements; three
    of length one, three of length two and one of length three,
    containing the asymmetry parameters for each separate component
    (see \code{\link{rmvevd}} and \bold{Examples}).}
  \item{model}{The specified model; a character string. Must be
    either \code{"log"} (the default) or \code{"alog"} (or any unique
    partial match), for the logistic and asymmetric logistic models
    respectively. The definition of each model is given (for general
    dimensions) in \code{\link{rmvevd}}.}
  \item{plot}{Logical; if \code{TRUE} the function is plotted. The
    minimum (evaluated) value is returned invisibly.
    If \code{FALSE} (the default), the following arguments are
    ignored.}
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
    the order defined by \code{ord}.}
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
  Let \eqn{z = (z_1,z_2,z_3)}{z = (z1,z2,z3)} and
  \eqn{w = (w_1,w_2,w_3)}{w = (w1,w2,w3)}.
  Any trivariate extreme value distribution can be written as
  \deqn{G(z) = \exp\left\{- \left\{\sum\nolimits_{j=1}^{3} y_j
    \right\} A\left(\frac{y_1}{\sum\nolimits_{j=1}^{3} y_j},
    \frac{y_2}{\sum\nolimits_{j=1}^{3} y_j},
    \frac{y_3}{\sum\nolimits_{j=1}^{3} y_j}\right)\right\}}{
    G(z) = exp{-(y1+y2+y3) A[y1/(y1+y2+y3),y2/(y1+y2+y3),
      y3/(y1+y2+y3)]}}
  for some function \eqn{A} defined on the simplex
  \eqn{S_3 = \{w \in R^3_+ : \sum\nolimits_{j=1}^{3}
    w_j = 1\}}{S_3 = {w: w1 + w2 + w3 = 1}},
  where \deqn{y_i = \{1+s_i(z_i-a_i)/b_i\}^{-1/s_i}}{
    yi = {1+si(zi-ai)/bi}^(-1/si)}
  for \eqn{1+s_i(z_i-a_i)/b_i > 0}{1+si(zi-ai)/bi > 0} and
  \eqn{i = 1,2,3}, and where the (generalized extreme value)
  marginal parameters are given by
  \eqn{(a_i,b_i,s_i)}{(ai,bi,si)}, \eqn{b_i > 0}{bi > 0}.
  If \eqn{s_i = 0}{si = 0} then \eqn{y_i}{yi} is defined by
  continuity.

  \eqn{A} is called (by some authors) the dependence function.
  It follows that \eqn{A(1,0,0) = A(0,1,0) = A(0,0,1) = 1}, and
  that \eqn{A} is a convex function with
  \eqn{\max(w_1,w_2,w_3) \leq A(w)\leq 1}{
    max(w1,w2,w3) <= A(w) <= 1} for
  all \eqn{w} in \eqn{S_3}.
  The lower and upper limits of \eqn{A} are obtained under complete
  dependence and mutual independence respectively.
  \eqn{A} does not depend on the marginal parameters.
}

\value{
  \code{atvpar} calculates or plots the dependence function
  for the trivariate logistic and trivariate asymmetric logistic
  models, at specified parameter values.
}

\seealso{\code{\link{atvnonpar}}, \code{\link{abvpar}},
  \code{\link{rmvevd}}, \code{\link{image}}}

\examples{
atvpar(dep = 0.5, model = "log")
s3pts <- matrix(rexp(30), nrow = 10, ncol = 3)
s3pts <- s3pts/rowSums(s3pts)
atvpar(s3pts, dep = 0.5, model = "log")
\dontrun{atvpar(dep = 0.05, model = "log", plot = TRUE, blty = 1)}
atvpar(dep = 0.95, model = "log", plot = TRUE, lower = 0.94)

asy <- list(.4, .1, .6, c(.3,.2), c(.1,.1), c(.4,.1), c(.2,.3,.2))
atvpar(s3pts, dep = 0.15, asy = asy, model = "alog")
atvpar(dep = 0.15, asy = asy, model = "al", plot = TRUE, lower = 0.7)
}
\keyword{distribution}








