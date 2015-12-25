\name{evind.test}
\alias{evind.test}
\title{Perform Hypothesis Test Of Independence}
\description{
  Perform score and likelihood ratio tests of independence for
	bivariate data, assuming a logistic dependence model as the
	alternative.
}
\usage{
evind.test(x, method = c("ratio", "score"), verbose = FALSE)
}
\arguments{
  \item{x}{A matrix or data frame, ordinarily with two columns,
    which may contain missing values.}
  \item{method}{The test methodology; either \code{"ratio"} for
	  the likelihood ratio test or \code{"score"} for the score test.}
  \item{verbose}{If \code{TRUE}, shows estimates of the marginal 
	  parameters in addition to the dependence parameter.}
}
\details{
  This simple function fits a stationary bivariate logistic model to the
	data and performs a hypothesis test of \eqn{\code{dep} = 1} versus
	\eqn{\code{dep} < 1} using the methodology in Tawn (1988). The null 
	distributions for the printed test statistics are chi-squared on one
	df for the likelihood ratio test, and standard normal for the score 
	test.
}
\value{
  An object of class \code{"htest"}.
}
\references{
  Tawn, J. A. (1988)
  Bivariate extreme value theory: models and estimation.
  \emph{Biometrika}, \bold{75}, 397--415.
}

\seealso{\code{\link{fbvevd}}, \code{\link{t.test}}}

\examples{
evind.test(sealevel)
evind.test(sealevel, method = "score")
}
\keyword{htest}

