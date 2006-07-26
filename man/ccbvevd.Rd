\name{ccbvevd}
\alias{ccbvevd}

\title{Calculate Conditional Copulas for Parametric Bivariate Extreme
  Value Distributions}
\description{
  Conditional copula functions, conditioning on either margin,
  for nine parametric bivariate extreme value models. 
}
\synopsis{
ccbvevd(x, mar = 2, dep, asy = c(1, 1), alpha, beta, model = c("log", 
    "alog", "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", 
    "amix"), lower.tail = TRUE) 
}
\usage{
ccbvevd(x, mar = 2, dep, asy = c(1, 1), alpha, beta, model = "log",
    lower.tail = TRUE) 
}
\arguments{
  \item{x}{A matrix or data frame, ordinarily with two columns,
    which may contain missing values. A data frame may also
    contain a third column of mode \code{logical}, which
    itself may contain missing values (see \bold{Details}).}
  \item{mar}{One or two; conditions on this margin.}
  \item{dep}{Dependence parameter for the logistic, asymmetric
    logistic, Husler-Reiss, negative logistic and asymmetric
    negative logistic models.}
  \item{asy}{A vector of length two, containing the two asymmetry
    parameters for the asymmetric logistic and asymmetric negative
    logistic models.}
  \item{alpha, beta}{Alpha and beta parameters for the bilogistic,
    negative bilogistic, Coles-Tawn and asymmetric mixed models.}
  \item{model}{The specified model; a character string. Must be
    either \code{"log"} (the default), \code{"alog"}, \code{"hr"},
    \code{"neglog"}, \code{"aneglog"}, \code{"bilog"},
    \code{"negbilog"}, \code{"ct"} or \code{"amix"} (or any unique
    partial match), for the logistic, asymmetric logistic,
    Husler-Reiss, negative logistic, asymmetric negative logistic,
    bilogistic, negative bilogistic, Coles-Tawn and asymmetric
    mixed models respectively.  If parameter arguments are given
    that do not correspond to the specified model those arguments
    are ignored, with a warning.}
  \item{lower.tail}{Logical; if \code{TRUE} (default), the
    conditional distribution function is returned; the conditional
    survivor function is returned otherwise.}
}
\details{
  The function calculates \eqn{P(U1 < x1|U2 = x2)}, where
  \eqn{(U1,U2)} is a random vector with Uniform(0,1) margins
  and with a dependence structure given by the specified
  parametric model. By default, the values of \eqn{x1} and
  \eqn{x2} are given by the first and second columns of the
  argument \code{x}. If \code{mar = 1} then this is reversed.

  If \code{x} has a third column \code{x3} of mode logical, then
  the function returns \eqn{P(U1 < x1|U2 = x2,I = x3)}, according
  to inference proceedures derived by Stephenson and Tawn (2004).
  See \code{\link{fbvevd}}. This requires numerical integration,
  and hence will be slower.
  
  This function is mainly for internal use. It is used by
  \code{\link{plot.bvevd}} to calculate the conditional P-P
  plotting diagnostics.
}
\value{
  A numeric vector of probabilities.
}
\references{
  Stephenson, A. G. and Tawn, J. A. (2004)
  Exploiting Occurence Times in Likelihood Inference for
  Componentwise Maxima.
  \emph{Biometrika} \bold{92}(1), 213--217.
}
\seealso{\code{\link{rbvevd}}, \code{\link{fbvevd}},
  \code{\link{plot.bvevd}}}

\keyword{distribution}

