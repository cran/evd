\name{atvnonpar}
\alias{atvnonpar}
\title{Non-parametric Estimates for Dependence Functions of the
    Trivariate Extreme Value Distribution}
\description{
  Calculate or plot non-parametric estimates for the dependence function
  \eqn{A} of the trivariate extreme value distribution.
}
\synopsis{
atvnonpar(x = rep(1/3,3), data, nsloc1 = NULL, nsloc2 = NULL, nsloc3 = NULL,
    method = c("pickands", "deheuvels", "hall"), kmar = NULL, plot = FALSE,
    col = heat.colors(12), blty = 0, grid = if(blty) 150 else 50,
    lower = 1/3, ord = 1:3, lab = as.character(1:3), lcex = 1)
}
\usage{
atvnonpar(x = rep(1/3,3), data, nsloc1 = NULL, nsloc2 = NULL, nsloc3 = NULL,
    method = c("pickands", "deheuvels", "hall"), plot = FALSE,
    col = heat.colors(12), blty = 0, grid = if(blty) 150 else 50,
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
  \item{data}{A matrix or data frame with three columns, which may
    contain missing values.}
  \item{nsloc1, nsloc2, nsloc3}{A data frame with the same number of
    rows as \code{data}, for linear modelling of the location
    parameter on the first/second/third margin.
    The data frames are treated as covariate matrices, excluding the
    intercept.
    A numeric vector can be given as an alternative to a single column
    data frame.}
  \item{method}{The estimation method; a character string. Must be
    either \code{"pickands"} (the default), \code{"deheuvels"}
    or \code{"hall"} (or any unique partial match). The three
    estimators are very similar, and may not be distinguishable
    when plotted.}
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
    default \eqn{\code{lower} = 1/3} as this is the theoretical
    minimum of the dependence function of the trivariate extreme
    value distribution.}
  \item{ord}{A vector of length three, which should be a permutation
    of the set \eqn{\{1,2,3\}}{{1,2,3}}. The points \eqn{(1,0,0)},
    \eqn{(0,1,0)} and \eqn{(0,0,1)} (the vertices of the simplex)
    are depicted clockwise from the top in the order defined by
    \code{ord}.}
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
  value is used in estimating the generalized extreme value margins.

  The dependence function of the trivariate extreme value
  distribution is defined in \code{\link{atvpar}}.
  The function \code{\link{atvpar}} calculates and plots dependence
  functions of trivariate logistic and trivariate asymmetric
  logistic models.
}

\value{
  \code{atvnonpar} calculates or plots a non-parametric estimate of
  the dependence function of the trivariate extreme value distribution.
}

\seealso{\code{\link{atvpar}}, \code{\link{abvnonpar}},
  \code{\link{fgev}}}

\examples{
s3pts <- matrix(rexp(30), nrow = 10, ncol = 3)
s3pts <- s3pts/rowSums(s3pts)
sdat <- rmvevd(100, dep = 0.6, model = "log", d = 3)
atvnonpar(s3pts, sdat)

\dontrun{atvnonpar(data = sdat, plot = TRUE)}
\dontrun{atvnonpar(data = sdat, plot = TRUE, ord = c(2,3,1), lab = LETTERS[1:3])}
\dontrun{atvpar(dep = 0.6, model = "log", plot = TRUE)}
\dontrun{atvpar(dep = 0.6, model = "log", plot = TRUE, blty = 1)}
}
\keyword{nonparametric}








