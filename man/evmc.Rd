\name{evmc}
\alias{evmc}
\title{Simulate Markov Chains With Extreme Value
  Dependence Structures}
\description{
  Simulation of first order Markov chains, such that each pair
  of consecutive values has the dependence structure of one of
  nine parametric bivariate extreme value distributions. 
}
\usage{
evmc(n, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
    margins = c("uniform","rweibull","frechet","gumbel"))
}
\arguments{
  \item{n}{Number of observations.}
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
    bilogistic, negative bilogistic, Coles-Tawn and asymmetric mixed
    models respectively. The definition of each model is given in
    \code{\link{rbvevd}}. If parameter arguments are given that do
    not correspond to the specified model those arguments are
    ignored, with a warning.}
  \item{margins}{The marginal distribution of each value; a
    character string. Must be either \code{"uniform"} (the
    default), \code{"rweibull"}, \code{"frechet"} or
    \code{"gumbel"} (or any unique partial match), for the uniform,
    standard reversed Weibull, standard Gumbel and standard Frechet
    distributions respectively.}
}
\value{
  A numeric vector of length \code{n}.
}

\seealso{\code{\link{marma}}, \code{\link{rbvevd}}}

\examples{
evmc(100, alpha = 0.1, beta = 0.1, model = "bilog")
evmc(100, dep = 10, model = "hr", margins = "gum")
}
\keyword{distribution}




