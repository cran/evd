\name{mrlplot}
\alias{mrlplot}
\title{Empirical Mean Residual Life Plot}
\description{
  The empirical mean residual life plot. 
}
\usage{
mrlplot(data, tlim, nt = max(100, length(data)), lty = c(2,1,2),
    col = 1, conf = 0.95, main = "Mean Residual Life Plot", xlab =
    "Threshold", ylab = "Mean Excess", \dots)
}
\arguments{
  \item{data}{A numeric vector.}
  \item{tlim}{A numeric vector of length two, giving the limits for
    the thresholds at which the mean residual life plot is
    evaluated. If \code{tlim} is not given, sensible defaults
    are used.}
  \item{nt}{The number of thresholds at which the mean residual life
    plot is evaluated.}
  \item{lty, col}{Arguments passed to \code{matplot}. The first
    and last elements of \code{lty} correspond to the lower and
    upper confidence limits respectively. Use zero to supress.}
  \item{conf}{The (pointwise) confidence coefficient for the
    plotted confidence intervals.} 
  \item{main}{Plot title.}     
  \item{xlab, ylab}{x and y axis labels.}
  \item{\dots}{Other arguments to be passed to \code{matplot}.}
}
\details{
  The empirical mean residual life plot is the locus of points
  \deqn{\left(u,\frac{1}{n_u} \sum\nolimits_{i=1}^{n_u}
    (x_{(i)} - u) \right)}{{u,1/n_u \sum_{i=1}^{n_u} (x(i) - u)}}
  where \eqn{x_{(1)}, \dots, x_{(n_u)}}{x(1), \dots, x(n_u)} are
  the \eqn{n_u} observations that exceed the threshold \eqn{u}.
  If the exceedances of a threshold \eqn{u_0}{u0}
  are generalized Pareto, the empirical mean residual life plot
  should be approximately linear for \eqn{u > u_0}{u > u0}.
  
  The confidence intervals within the plot are symmetric intervals
  based on the approximate normality of sample means.
}
\value{
  A list with components \code{x} and \code{y} is invisibly returned.
  The components contain those objects that were passed to the formal
  arguments \code{x} and \code{y} of \code{matplot} in order to create
  the mean residual life plot.
}
\author{Stuart Coles and Alec Stephenson}  

\seealso{\code{\link{fpot}}, \code{\link{matplot}},
    \code{\link{tcplot}}}

\examples{
data(portpirie)
mrlplot(portpirie)
}
\keyword{hplot}

