\name{profile2d.evd}
\alias{profile2d}
\alias{profile2d.evd}

\title{Method for Profiling EVD Objects}
\description{
  Calculate joint profile traces for fitted models.
}
\usage{
\method{profile2d}{evd}(fitted, prof, which, pts = 20, convergence =
    FALSE, method = "Nelder-Mead", control = list(maxit = 10000), \dots)
}
\arguments{
  \item{fitted}{An object of class \code{"evd"}.}
  \item{prof}{An object of class \code{"profile.evd"}, created using
    \code{\link{profile.evd}} with argument \code{fitted}.
    The object must contain the (marginal) profile traces for the two
    parameters specified in \code{which}.}
  \item{which}{A character vector of length two containing the original 
     model parameters that are to be jointly profiled.}
   \item{pts}{The number of distinct values used for each profiled
     parameter in \code{which}.
    There are \code{pts^2} optimizations performed in total.}
  \item{convergence}{Logical; print convergence code after each
    optimization? (A warning is given for each non-zero convergence
    code, irrespective of the value of \code{convergence}.)}
  \item{method}{The optimization method.}
  \item{control}{Passed to \code{optim}. See \code{\link{optim}} for
    details.}
  \item{\dots}{Ignored.}
}
\value{
  An object of class \code{"profile2d.evd"}, which is a list with three
  elements.
  The first element, a matrix named \code{trace}, has the same structure
  as the elements of an object of class \code{"profile.evd"}.
  The last two elements give the distinct values used for each profiled
  parameter in \code{which}.
}

\seealso{\code{\link{profile.evd}}, \code{\link{plot.profile2d.evd}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{M1P <- profile(M1)}
\dontrun{M1JP <- profile2d(M1, M1P, which = c("scale", "shape"))}
\dontrun{plot(M1JP)}
}
\keyword{models}
