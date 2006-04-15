\name{mtransform}
\alias{mtransform}
\title{GEV Transformations}
\description{
  Transforms to exponential margins under the GEV model. 
}
\usage{
mtransform(x, p, inv = FALSE, drp = FALSE)
}
\arguments{
  \item{x}{A matrix with n rows and d columns, or a vector. In
    the latter case, if \code{p} is a list with the same length
    as the vector, it is treated as a matrix with one row. If
    \code{p} is not a list, it is treated as a matrix with one
    column.}
  \item{p}{A vector of length three or a matrix with n rows and
    three columns. It can also be a list of length d, in which
    case each element must be a vector of length three or a matrix
    with n rows and three columns.}
  \item{inv}{Logical; use the inverse transformation?}
  \item{drp}{Logical; return a vector rather than a single row
    matrix?. Note that a single column matrix is always returned
    as a vector.}
}
\details{
  Let \eqn{x_i} denote a vector of observations for
  \eqn{i = 1,\ldots,n}.
  This function implements the transformation
  \deqn{y_{i} = \{1+s_i(x_{i}-a_i)/b_i\}_{+}^{-1/s_i}}
  to each column of the matrix \code{x}.

  The values \eqn{(a_i,b_i,s_i)} are contained in the ith
  row of the n by 3 matrix \code{p}. If \code{p} is a vector
  of length three, the parameters are the same for every
  \eqn{i = 1,\ldots,n}. Alternatively, \code{p} can be a list
  with d elements, in which case the jth element is used to
  transform the jth column of \code{x}.

  This function is mainly for internal use. It is used by
  bivariate and multivariate routines to calculate marginal
  transformations.
}
\value{
  A numeric matrix or vector.
}

\keyword{manip}

