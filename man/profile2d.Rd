\name{profile2d}
\alias{profile2d}

\title{Generic Function for Profiling Models}
\description{
  Investigates behavior of objective function near the solution
  represented by \code{fitted}.
  See documentation on method functions for further details.
}
\usage{
profile2d(fitted, ...)
}
\arguments{
  \item{fitted}{The original fitted model object.}
  \item{\dots}{Additional arguments. See documentation on individual
    methods.}
}
\value{
  A list containing information on joint profile likelihood/deviance
  surfaces. See the individual methods for further details.
}

\seealso{\code{\link{profile2d.evd}}}

\keyword{models}