\name{bvtcplot}
\alias{bvtcplot}
\title{Bivariate Threshold Choice Plot}
\description{
  Produces a diagnostic plot to assist with threshold choice for
	bivariate data. 
}
\usage{
bvtcplot(x, spectral = FALSE, xlab, ylab, \dots)
}
\arguments{
  \item{x}{A matrix or data frame, ordinarily with two columns,
    which may contain missing values.}
  \item{spectral}{If \code{TRUE}, an estimate of the spectral measure is
	  plotted instead of the diagnostic plot.}
  \item{ylab, xlab}{Graphics parameters.}
  \item{\dots}{Other arguments to be passed to the plotting
    function.}
}
\details{
  If \code{spectral} is \code{FALSE} (the default), produces a threshold
	choice plot as illustrated in Beirlant et al. (2004). With \eqn{n} 
	non-missing bivariate observations, the integers
	\eqn{k = 1,\dots,n-1}{k = 1,...,n-1} are plotted against the values
	\eqn{(k/n)r_{(n-k)}}{(k/n)r_(n-k)}, where \eqn{r_{(n-k)}}{r_(n-k)}
	is the \eqn{(n-k)}th order statistic of the sum of the margins
  following empirical transformation to standard Frechet.	
	
	A vertical line is drawn at \code{k0}, where \code{k0} is the largest
	integer for which the y-axis is above the value two. If \code{spectral} 
	is \code{FALSE}, the largest \code{k0} data points are used to plot an
	estimate of the spectal measure \eqn{H([0, w])} versus \eqn{w}.
}
\value{
  A list is invisibly returned giving \code{k0} and the values used to
	produce the plot.
}
\references{
  Beirlant, J., Goegebeur, Y., Segers, J. and Teugels, J. L. (2004)
  \emph{Statistics of Extremes: Theory and Applications.},
  Chichester, England: John Wiley and Sons.
}

\seealso{\code{\link{fbvpot}}, \code{\link{tcplot}}}

\examples{
\dontrun{bvtcplot(lossalae)}
\dontrun{bvtcplot(lossalae, spectral = TRUE)}
}
\keyword{hplot}

