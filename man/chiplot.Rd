\name{chiplot}
\alias{chiplot}
\title{Dependence Measure Plots}
\description{
  Plots of estimates of the dependence measures chi and chi-bar
  for bivariate data. 
}
\usage{
chiplot(data, nq = 100, qlim = NULL, which = 1:2, conf = 0.95, lty = 1,
    cilty = 2, col = 1, cicol = 1, xlim = c(0,1), ylim1 = NULL,
    ylim2 = c(-1,1), main1 = "Chi Plot", main2 = "Chi Bar Plot", xlab =
    "Quantile", ylab1 = "Chi", ylab2 = "Chi Bar", ask = nb.fig <
    length(which) && dev.interactive(), \dots)
}
\arguments{
  \item{data}{A matrix or data frame with two columns. Rows
    (observations) with missing values are stripped from
    the data before any computations are performed.}
  \item{nq}{The number of quantiles at which the measures
    are evaluated.}
  \item{qlim}{The limits of the quantiles at which the measures
    are evaluated (see \bold{Details}).}
  \item{which}{If only one plot is required, specify \code{1}
    for chi and \code{2} for chi-bar.}
  \item{conf}{The confidence coefficient of the plotted confidence
    intervals.}
  \item{lty, cilty}{Line types for the estimates of the measures and
    for the confidence intervals respectively. Use zero to supress.}
  \item{col, cicol}{Colour types for the estimates of the measures
    and for the confidence intervals respectively.}
  \item{xlim, xlab}{Limits and labels for the x-axis; they apply
    to both plots.}
  \item{ylim1}{Limits for the y-axis of the chi plot. If this
    is \code{NULL} (the default) the upper limit is one, and the
    lower limit is the minimum of zero and the smallest plotted
    value.}
  \item{ylim2}{Limits for the y-axis of the chi-bar plot.}
  \item{main1, main2}{The plot titles for the chi and chi-bar plots
    respectively.}
  \item{ylab1, ylab2}{The y-axis labels for the chi and chi-bar plots
    respectively.} 
  \item{ask}{Logical; if \code{TRUE}, the user is asked before
    each plot.}
  \item{\dots}{Other arguments to be passed to \code{matplot}.}
}

\details{
  These measures are explained in full detail in Coles, Heffernan
  and Tawn (1999). A brief treatment is also given in Section
  8.4 of Coles(2001).
  A short summary is given as follows.
  We assume that the data are \emph{iid} random vectors with common
  bivariate distribution function \eqn{G}, and we define the random
  vector \eqn{(X,Y)} to be distributed according to \eqn{G}.
  
  The chi plot is a plot of \eqn{q} against empirical estimates of
  \deqn{\chi(q) = 2 - \log(\Pr(F_X(X) < q, F_Y(Y) < q)) / \log(q)}{
    chi(q) = 2 - log(Pr(F_X(X) < q, F_Y(Y) < q)) / log(q)}
  where \eqn{F_X} and \eqn{F_Y} are the marginal distribution
  functions, and where \eqn{q} is in the interval (0,1).
  The quantity \eqn{\chi(q)}{chi(q)} is bounded by
  \deqn{2 - \log(2u - 1)/\log(u) \leq \chi(q) \leq 1}{
    2 - log(2u - 1)/log(u) <= chi(q) <= 1}
  where the lower bound is interpreted as \code{-Inf} for
  \eqn{q \leq 1/2}{q <= 1/2} and zero for \eqn{q = 1}.
  These bounds are reflected in the corresponding estimates.

  The chi bar plot is a plot of \eqn{q} against empirical estimates of
  \deqn{\bar{\chi}(q) = 2\log(1-q)/\log(\Pr(F_X(X) > q, F_Y(Y) > q)) - 1}{
    chibar(q) = 2log(1-q)/log(Pr(F_X(X) > q, F_Y(Y) > q)) - 1}
  where \eqn{F_X} and \eqn{F_Y} are the marginal distribution
  functions, and where \eqn{q} is in the interval (0,1).
  The quantity \eqn{\bar{\chi}(q)}{chibar(q)} is bounded by
  \eqn{-1 \leq \bar{\chi}(q) \leq 1}{-1 <= chibar(q) <= 1}
  and these bounds are reflected in the corresponding estimates.

  Note that the empirical estimators for \eqn{\chi(q)}{chi(q)} and
  \eqn{\bar{\chi}(q)}{chibar(q)} are undefined near \eqn{q=0} and \eqn{q=1}. By
  default the function takes the limits of \eqn{q} so that the plots
  depicts all values at which the estimators are defined. This can be
  overridden by the argument \code{qlim}, which must represent a subset
  of the default values (and these can be determined using the
  component \code{quantile} of the invisibly returned list; see
  \bold{Value}).
  
  The confidence intervals within the plot assume that observations are
  independent, and that the marginal distributions are estimated exactly.
  The intervals are constructed using the delta method; this may
  lead to poor interval estimates near \eqn{q=0} and \eqn{q=1}.
  
  The function \eqn{\chi(q)}{chi(q)} can be interpreted as a quantile
  dependent measure of dependence. In particular, the sign of
  \eqn{\chi(q)}{chi(q)} determines whether the variables are positively
  or negatively associated at quantile level \eqn{q}.
  By definition, variables are said to be asymptotically independent
  when \eqn{\chi(1)}{chi(1)} (defined in the limit) is zero.
  For independent variables, \eqn{\chi(q) = 0}{chi(q) = 0} for all
  \eqn{q} in (0,1).
  For perfectly dependent variables, \eqn{\chi(q) = 1}{chi(q) = 1}
  for all \eqn{q} in (0,1).
  For bivariate extreme value distributions, \eqn{\chi(q) =
    2(1-A(1/2))}{chi(q) = 2(1-A(1/2))}
  for all \eqn{q} in (0,1), where \eqn{A} is the dependence function,
  as defined in \code{\link{abvevd}}. If a bivariate threshold model
  is to be fitted (using \code{\link{fbvpot}}), this plot can therefore
  act as a threshold identification plot, since e.g. the use of 95\%
  marginal quantiles as threshold values implies that \eqn{\chi(q)}{chi(q)}
  should be approximately constant above \eqn{q = 0.95}.

  The function \eqn{\bar{\chi}(q)}{chibar(q)} can again be interpreted
  as a quantile dependent measure of dependence; it is most useful
  within the class of asymptotically independent variables.
  For asymptotically dependent variables (i.e. those for which
  \eqn{\chi(1) < 1}{chi(1) < 1}), we have \eqn{\bar{\chi}(1) = 1}{
    chibar(1) = 1}, where
  \eqn{\bar{\chi}(1)}{chibar(1)} is again defined in the limit.
  For asymptotically independent variables, \eqn{\bar{\chi}(1)}{
    chibar(1)} provides a measure that increases with dependence strength.
  For independent variables \eqn{\bar{\chi}(q) = 0}{chibar(q) = 0} for
  all \eqn{q} in (0,1), and hence \eqn{\bar{\chi}(1) = 0}{chibar(1) = 0}.
}
\value{
  A list with components \code{quantile}, \code{chi} (if \code{1} is in
  \code{which}) and \code{chibar} (if \code{2} is in \code{which})
  is invisibly returned.
  The components \code{quantile} and \code{chi} contain those objects
  that were passed to the formal arguments \code{x} and \code{y} of
  \code{matplot} in order to create the chi plot.
  The components \code{quantile} and \code{chibar} contain those objects
  that were passed to the formal arguments \code{x} and \code{y} of
  \code{matplot} in order to create the chi-bar plot.
}
\references{
  Coles, S. G., Heffernan, J. and Tawn, J. A. (1999)
  Dependence measures for extreme value analyses.
  \emph{Extremes}, \bold{2}, 339--365.

  Coles, S. G. (2001)
  \emph{An Introduction to Statistical Modelling of Extreme Values},
  London: Springer--Verlag.  
}
\author{Jan Heffernan and Alec Stephenson}  

\seealso{\code{\link{fbvevd}}, \code{\link{fbvpot}},
  \code{\link{matplot}}}

\examples{
par(mfrow = c(1,2))
smdat1 <- rbvevd(1000, dep = 0.6, model = "log")
smdat2 <- rbvevd(1000, dep = 1, model = "log")
chiplot(smdat1)
chiplot(smdat2)
}
\keyword{hplot}

