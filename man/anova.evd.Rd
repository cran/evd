\name{anova.evd}
\alias{anova.evd}

\title{Compare Nested EVD Objects}
\description{
  Compute an analysis of deviance table for two or more nested
  evd objects.
}
\usage{
\method{anova}{evd}(object, object2, \dots)
}
\arguments{
  \item{object}{An object of class \code{"evd"}.}
  \item{object2}{An object of class \code{"evd"} that
    represents a model nested within \code{object}.}
  \item{\dots}{Further successively nested objects.}
}
\value{
  An object of class \code{c("anova", "data.frame")}, with one
  row for each model, and the following five columns
  \item{M.Df}{The number of parameters.}
  \item{Deviance}{The deviance.}
  \item{Df}{The number of parameters of the model in the previous row
    minus the number of parameters.}
  \item{Chisq}{The deviance minus the deviance of the model
    in the previous row.}
  \item{Pr(>chisq)}{The p-value calculated by comparing the quantile
    \code{Chisq} with a chi-squared distribution on \code{Df} degrees
    of freedom.}
}
\section{Warning}{
  Circumstances may arise such that the asymptotic distribution of the
  test statistic is not chi-squared.
  In particular, this occurs when the nested
  model is constrained at the edge of the parameter space. It is up
  to the user recognize this, and to interpret the output correctly.
}

\seealso{\code{\link{fbvevd}}, \code{\link{fextreme}},
  \code{\link{fgev}}, \code{\link{forder}}}

\examples{
uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
trend <- (-49:50)/100
M1 <- fgev(uvdata, nsloc = trend)
M2 <- fgev(uvdata)
M3 <- fgev(uvdata, shape = 0)
anova(M1, M2, M3)

bvdata <- rbvevd(100, dep = 0.75, model = "log")
M1 <- fbvevd(bvdata, model = "log")
M2 <- fbvevd(bvdata, model = "log", dep = 0.75)
anova(M1, M2)
}

\keyword{models}

