\name{confint.evd}
\alias{confint.evd}
\alias{confint.profile.evd}
\title{Calculate Confidence Intervals}
\description{
  Calculate profile and Wald confidence intervals of parameters
  in fitted models.
}
\usage{
\method{confint}{evd}(object, parm, level = 0.95, \dots)
\method{confint}{profile.evd}(object, parm, level = 0.95, \dots)
}
\arguments{
  \item{object}{Either a fitted model object (of class \code{evd})
    for Wald confidence intervals, or a profile trace (of class
    \code{profile.evd}) for profile likelihood confidence intervals.} 
  \item{parm}{A character vector of parameters; a confidence
    interval is calculated for each parameter. If missing, then
    intervals are returned for all parameters in the fitted model
    or profile trace.}
  \item{level}{A single number giving the confidence level.}
  \item{\dots}{Not used.}
}
\value{
  A matrix with two columns giving lower and upper confidence
  limits.
  
  For profile confidence intervals, this function assumes that the 
  profile trace is unimodal. If the profile trace is not unimodal then
  the function will give spurious results.  
}

\seealso{\code{\link{profile.evd}}}

\examples{
m1 <- fgev(portpirie)
confint(m1)
\dontrun{pm1 <- profile(m1)}
\dontrun{plot(pm1)}
\dontrun{confint(pm1)}
}
\keyword{manip}

