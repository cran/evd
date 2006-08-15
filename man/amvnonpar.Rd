\name{amvnonpar}
\alias{amvnonpar}
\title{Non-parametric Estimates for Dependence Functions of the
    Multivariate Extreme Value Distribution}
\description{
  Calculate non-parametric estimates for the dependence function
  \eqn{A} of the multivariate extreme value distribution and plot
  the estimated function in the trivariate case.
}
\synopsis{
amvnonpar(x = rep(1/d,d), data, d = 3, epmar = FALSE, nsloc = NULL,
    madj = 0, kmar = NULL, plot = FALSE, col = heat.colors(12),
    blty = 0, grid = if(blty) 150 else 50, lower = 1/3, ord = 1:3,
    lab = as.character(1:3), lcex = 1)
}
\usage{
amvnonpar(x = rep(1/d,d), data, d = 3, epmar = FALSE, nsloc = NULL,
    madj = 0, plot = FALSE, col = heat.colors(12), blty = 0,
    grid = if(blty) 150 else 50, lower = 1/3, ord = 1:3,
    lab = as.character(1:3), lcex = 1)
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
  \item{data}{A matrix or data frame with \code{d} columns, which may
    contain missing values.}
  \item{d}{The dimension; an integer greater than or equal to two.
    The trivariate case \code{d = 3} is the default.}
  \item{epmar}{If \code{TRUE}, an empirical transformation of the
    marginals is performed in preference to marginal parametric
    GEV estimation, and the \code{nsloc} argument is ignored.}
  \item{nsloc}{A data frame with the same number of rows as \code{data},
    or a list containing \code{d} elements of this type, for linear
    modelling of the marginal location parameters. In the former case,
    the argument is applied to all margins. The data frames are treated
    as covariate matrices, excluding the intercept. Numeric vectors can
    be given as alternatives to single column data frames. A list can
    contain \code{NULL} elements for stationary modelling of selected
    margins.}
  \item{madj}{Performs marginal adjustments. See
    \code{\link{abvnonpar}}.}
  \item{plot}{Logical; if \code{TRUE}, and the dimension \code{d} is
    three (the default dimension), the dependence function of a
    trivariate extreme value distribution is plotted. For plotting in
    the bivariate case, use \code{\link{abvnonpar}}. If \code{FALSE}
    (the default), the following arguments are ignored.}
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
    default \eqn{\code{lower} = 1/3} as this is the theoretical
    minimum of the dependence function of the trivariate extreme
    value distribution.}
  \item{ord}{A vector of length three, which should be a permutation
    of the set \eqn{\{1,2,3\}}{{1,2,3}}. The points \eqn{(1,0,0)},
    \eqn{(0,1,0)} and \eqn{(0,0,1)} (the vertices of the simplex)
    are depicted clockwise from the top in the order defined by
    \code{ord}. The argument alters the way in which the function
    is plotted; it does not change the function definition.}
  \item{lab}{A character vector of length three, in which case the
    \code{i}th margin is labelled using the \code{i}th component,
    or \code{NULL}, in which case no labels are given. By default,
    \code{lab} is \code{as.character(1:3)}. The actual location of
    the margins, and hence the labels, is defined by \code{ord}.}
  \item{lcex}{A numerical value giving the amount by which the
    labels should be scaled relative to the default. Ignored
    if \code{lab} is \code{NULL}.}
}

\note{
  The rows of \code{data} that contain missing values are not used
  in the estimation of the dependence structure, but every non-missing
  value is used in estimating the margins.

  The dependence function of the multivariate extreme value
  distribution is defined in \code{\link{amvevd}}.
  The function \code{\link{amvevd}} calculates and plots dependence
  functions of multivariate logistic and multivariate asymmetric
  logistic models.

  The estimator plotted or calculated is a multivariate extension of
  the bivariate Pickands estimator defined in \code{\link{abvnonpar}}.
}

\value{
  A numeric vector of estimates. If plotting, the smallest evaluated
  estimate is returned invisibly.
}

\seealso{\code{\link{amvevd}}, \code{\link{abvnonpar}},
  \code{\link{fgev}}}

\examples{
s5pts <- matrix(rexp(50), nrow = 10, ncol = 5)
s5pts <- s5pts/rowSums(s5pts)
sdat <- rmvevd(100, dep = 0.6, model = "log", d = 5)
amvnonpar(s5pts, sdat, d = 5)

\dontrun{amvnonpar(data = sdat, plot = TRUE)}
\dontrun{amvnonpar(data = sdat, plot = TRUE, ord = c(2,3,1), lab = LETTERS[1:3])}
\dontrun{amvevd(dep = 0.6, model = "log", plot = TRUE)}
\dontrun{amvevd(dep = 0.6, model = "log", plot = TRUE, blty = 1)}
}
\keyword{nonparametric}








