\name{plot.profile.evd}
\alias{plot.profile.evd}

\title{Plot Profile Devainces}
\description{
  Displays profile deviances from a model profiled with
  \code{\link{profile.evd}}.
}
\usage{
\method{plot}{profile.evd}(x, which = names(x), main = NULL,
    ask = nb.fig < length(which) && dev.interactive(), ci = 0.95,
    clty = 2, \dots)
}
\arguments{
  \item{x}{An object of class \code{"profile.evd"}.}
  \item{which}{Parameters for which the profile deviance is
    plotted. By default all profiled parameters in \code{x}
    are used.}
  \item{main}{Title of each plot; a character vector, the
    same length as \code{which}.}
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{ci}{By default the plots include a horizonal line
    representing a profile confidence interval with confidence
    coefficient \code{ci} (but see \bold{Warning}).}
  \item{clty}{The line type of the horizontal line that represents
    a profile confidence interval. Use zero to suppress.}
  \item{\dots}{Other graphics parameters.}
}
\details{
  The profile deviance function is minus twice the logarithm of
  the profile likelihood. 
}
\section{Warning}{
  The interval represented by the horizontal line may not be
  a confidence interval with (asymptotic) confidence coefficient
  \code{ci}, because the usual asymptotic properties of maximum
  likelihood estimators may not hold!
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

\seealso{\code{\link{plot.profile2d.evd}}, \code{\link{profile.evd}},
  \code{\link{profile2d.evd}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{M1P <- profile(M1)}
\dontrun{pcint(M1P)}
\dontrun{plot(M1P)}
}
\keyword{hplot}