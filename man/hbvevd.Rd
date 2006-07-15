\name{hbvevd}
\alias{hbvevd}
\title{Parametric Spectral Density Functions of Bivariate Extreme
    Value Models}
\description{
  Calculate or plot the density \eqn{h} of the spectral measure
  \eqn{H} on the interval \eqn{(0,1)}, for nine parametric
  bivariate extreme value models.
}
\usage{
hbvevd(x = 0.5, dep, asy = c(1,1), alpha, beta, model = c("log", "alog",
    "hr", "neglog", "aneglog", "bilog", "negbilog", "ct", "amix"),
     half = FALSE, plot = FALSE, add = FALSE, lty = 1, \dots)
}
\arguments{
  \item{x}{A vector of values at which the function is evaluated
    (ignored if plot or add is \code{TRUE}). \eqn{h(1/2)} is
    returned by default.}
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
    mixed models respectively. The definition of each model is given
    in \code{\link{rbvevd}}. If parameter arguments are given that do
    not correspond to the specified model those arguments are
    ignored, with a warning.}
  \item{half}{Logical; if \code{TRUE} the function is divided by two,
    corresponding to a spectral measure with total mass one rather
    than two.}
  \item{plot}{Logical; if \code{TRUE} the function is plotted. The
    x and y values used to create the plot are returned invisibly.}
  \item{add}{Logical; add to an existing plot?}
  \item{lty}{Line type.}
  \item{\dots}{Other high-level graphics parameters to be passed to
    \code{plot}.}
}
\details{
  Any bivariate extreme value distribution can be written as
  \deqn{G(z_1,z_2) = \exp\left[-\int_0^1 \max\{wy_1, (1-w)y_2\}
    H(dw)\right]}{G(z1,z2) = exp{-\int_0^1 max(w y1, (1-w) y2) H(dw)}}
  for some function \eqn{H(\cdot)}{H()} defined on \eqn{[0,1]},
  satisfying
  \deqn{\int_0^1 w H(dw) = \int_0^1 (1-w) H(dw) = 1}{
    \int_0^1 w H(dw) = int_0^1 (1-w) H(dw) = 1.}
  In particular, the total mass of H is two.
  The functions \eqn{y_1}{y1} and \eqn{y_2}{y2} are as defined in
  \code{\link{abvevd}}.

  \eqn{H(\cdot)}{H()} is called the spectral measure, with density
  \eqn{h} on the interval \eqn{(0,1)}.
}
\section{Point Masses}{
  For differentiable models H may have up to two point masses:
  at zero and one. Assuming that the model parameters are in the
  interior of the parameter space, we have the following. For the
  asymmetric logistic and asymmetric negative logistic models the
  point masses are of size \code{1-asy1} and \code{1-asy2}
  respectively. For the asymmetric mixed model they are of size
  \code{1-alpha-beta} and \code{1-alpha-2*beta} respectively. For
  all other models the point masses are zero.

  At independence, H has point masses of size one at both
  zero and one. At complete dependence [a non-differentiable
  model] H has a single point mass of size two at \eqn{1/2}.
  In either case, \eqn{h} is zero everywhere.
}

\value{
  \code{hbvevd} calculates or plots the spectral density function
  \eqn{h} for one of nine parametric bivariate extreme value models,
  at specified parameter values.
}

\seealso{\code{\link{abvevd}}, \code{\link{fbvevd}},
  \code{\link{rbvevd}}, \code{\link{plot.bvevd}}}

\examples{
hbvevd(dep = 2.7, model = "hr")
hbvevd(seq(0.25,0.5,0.75), dep = 0.3, asy = c(.7,.9), model = "alog")
hbvevd(alpha = 0.3, beta = 1.2, model = "negbi", plot = TRUE)

bvdata <- rbvevd(100, dep = 0.7, model = "log")
M1 <- fitted(fbvevd(bvdata, model = "log"))
hbvevd(dep = M1["dep"], model = "log", plot = TRUE)
}
\keyword{distribution}








