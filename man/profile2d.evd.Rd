\name{profile2d.evd}
\alias{profile2d.evd}

\title{Method for Profiling evd Objects}
\description{
  Investigates behavior of objective function near the solution
  represented by \code{fitted}.
}
\usage{
\method{profile2d}{evd}(fitted, prof, which,  pts = 20,
    convergence = FALSE, control = list(maxit = 5000), \dots)
}
\arguments{
  \item{fitted}{The original fitted model object.}
  \item{prof}{An object of class \code{"profile.evd"}, created using
    \code{\link{profile.evd}} with argument \code{fitted}.
    The object must contain the (marginal) profiles for the two
    parameters specified in \code{which}.}
  \item{which}{A character vector of length two containing the original 
     model parameters that are to be jointly profiled.}
   \item{pts}{The number of distinct values used for each profiled
     parameter in \code{which}.
    There are \code{pts^2} optimizations performed in total.}
  \item{convergence}{Logical; print convergence code after each
    optimization?}
  \item{control}{Passed to \code{optim}. See \code{\link{optim}} for
    details.}
  \item{\dots}{Ignored.}
}
\value{
  An object of class \code{"profile2d.evd"};
  a list with three elements.
  The first element, a matrix named \code{trace}, has the same structure
  as the elements of an object of class \code{"profile.evd"}.
  The last two elements give the distinct values used for each profiled
  parameter in \code{which}.
}

\seealso{\code{\link{profile}}, \code{\link{profile.evd}},
  \code{\link{profile2d}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{M1P <- profile(M1)}
\dontrun{M1JP <- profile2d(M1, M1P, which = c("scale", "shape"))}
\dontrun{plot(M1JP)}
}
\keyword{models}