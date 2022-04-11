\name{plot.profile.evd}
\alias{plot.profile.evd}

\title{Plot Profile Log-likelihoods}
\description{
  Displays profile log-likelihoods from a model profiled with
  \code{\link{profile.evd}}.
}
\usage{
\method{plot}{profile.evd}(x, which = names(x), main = NULL,
    ask = nb.fig < length(which) && dev.interactive(), ci = 0.95,
    clty = 2, \dots)
}
\arguments{
  \item{x}{An object of class \code{"profile.evd"}.}
  \item{which}{A character vector giving the parameters for which
    the profile deviance is plotted, and for which profile confidence
    intervals are calculated. By default all profiled parameters in
    \code{x} are used.}
  \item{main}{Title of each plot; a character vector, the
    same length as \code{which}.}
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{ci}{A numeric vector. For each parameter in \code{which}
    profile confidence intervals are calculated, for each confidence
    coefficient in \code{ci} (but see \bold{Warning}).
    The intervals are returned invisibly as a list of
    vectors/matrices. Each plot then (by default) includes
    horizonal lines that represent each interval.}
  \item{clty}{The line type of the horizontal lines that represent
    the profile confidence intervals. To omit the lines set
    \code{clty} to zero.}
  \item{\dots}{Other graphics parameters.}
}
\value{
  Profile devainces are plotted for each parameter in
  \code{which}. For calculation of profile confidence intervals,
  use the \code{\link{confint.profile.evd}} function.
}
\section{Warning}{
  The profile confidence intervals may not have confidence coefficient
  \code{ci}, because the usual asymptotic properties of maximum
  likelihood estimators may not hold.
  For the GEV model, the usual asymptotic properties hold when the
  shape parameter is greater than \eqn{-0.5} (Smith, 1985).
}
\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}

\seealso{\code{\link{confint.profile.evd}}, \code{\link{plot.profile2d.evd}}, 
  \code{\link{profile.evd}}, \code{\link{profile2d.evd}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
M1 <- fgev(uvdata)
\dontrun{M1P <- profile(M1)}
\dontrun{par(mfrow = c(2,2))}
\dontrun{cint <- plot(M1P, ci = c(0.95, 0.99))}
\dontrun{cint}
}
\keyword{hplot}
