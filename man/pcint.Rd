\name{pcint}
\alias{pcint}

\title{Profile Likelihood Confidence Intervals}
\description{
  Calculates profile likelihood confidence intervals from an
  object of class \code{"profile.evd"}.
}
\usage{
pcint(prof, which = names(prof), ci = 0.95)
}
\arguments{
  \item{prof}{An object of class \code{\link{profile.evd}}.}
  \item{which}{The original model parameters for which profile
    likelihood confidence intervals are calculated.
    By default intervals are calculated for all the parameters
    that have been profiled in \code{prof}.}
  \item{ci}{Confidence coefficient (but see \bold{Warning}).}
}
\value{
  A list with one element for each parameter in \code{which}.
  Each element is a numeric vector of length two giving the
  lower and upper limits of the confidence intervals
  (but see \bold{Warning}).
}
\section{Warning}{
  The calculated intervals may not be confidence intervals with
  (asymptotic) confidence coefficient \code{ci}, because the
  usual asymptotic properties of maximum likelihood estimators
  may not hold!
  The usual asymptotic properties hold when the shape parameter(s)
  is(are) greater than \eqn{-0.5}, and when
  the parameters are not on the edge of the parameter space
  (Smith, 1985). Fortunately, this is usually the case.
}
\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}
  
\seealso{\code{\link{profile.evd}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{M1P <- profile(M1)}
\dontrun{pcint(M1P)}
\dontrun{plot(M1P)}
}
\keyword{models}



