\name{clusters}
\alias{clusters}
\title{Identify Clusters of Exceedences}
\description{
  Identify clusters of exceedences. 
}
\usage{
clusters(data, u, r = 1, ulow = -Inf, rlow = 1, cmax = FALSE, keep.names
    = TRUE, plot = FALSE, xdata = seq(along = data), lvals = TRUE, lty =
    1, lwd = 1, pch = par("pch"), col = if(n > 250) NULL else "grey",
    xlab = "Index", ylab = "Data", ...)
}
\arguments{
  \item{data}{A numeric vector, which may contain missing values.}
  \item{u}{A single value giving the threshold, unless a time varying
    threshold is used, in which case \code{u} should be a vector of
    thresholds, typically with the same length as \code{data} (or else
    the usual recycling rules are applied).}
  \item{r}{A postive integer denoting the clustering interval
    length. By default the interval length is one.}
  \item{ulow}{A single value giving the lower threshold, unless a time
    varying lower threshold is used, in which case \code{ulow} should
    be a vector of lower thresholds, typically with the same length as
    \code{data} (or else the usual recycling rules are applied).
    By default there is no lower threshold (or equivalently, the
    lower threshold is \code{-Inf}).}
  \item{rlow}{A postive integer denoting the lower clustering interval
    length. By default the interval length is one.}
  \item{cmax}{Logical; if \code{FALSE} (the default), a list
    containing the clusters of exceedences is returned. If
    \code{TRUE} a numeric vector containing the cluster maxima
    is returned.}
  \item{keep.names}{Logical; if \code{FALSE}, the function makes
    no attempt to retain the names/indices of the observations
    within the returned object. If \code{data} contains a large
    number of observations, this can make the function run much
    faster. The argument is mainly designed for internal use.}
  \item{plot}{Logical; if \code{TRUE} a plot is given that depicts
    the identified clusters, and the clusters (if \code{cmax} is
    \code{FALSE}) or cluster maxima (if \code{cmax} is \code{TRUE})
    are returned invisibly. If \code{FALSE} (the default), the
    following arguments are ignored.}
  \item{xdata}{A numeric vector with the same length as \code{data},
    giving the values to be plotted on the x-axis.}
  \item{lvals}{Logical; should the values below the threshold and
    the line depicting the lower threshold be plotted?}
  \item{lty, lwd}{Line type and width for the lines depicting the
    threshold and the lower threshold.}
  \item{pch}{Plotting character.}
  \item{col}{Strips of colour \code{col} are used to identify the
    clusters. An observation is contained in the cluster if the 
    centre of the corresponding plotting character is contained
    in the coloured strip. If \code{col} is \code{NULL} the strips
    are omitted. By default the strips are coloured \code{"grey"},
    but are omitted whenever \code{data} contains more than 250
    observations.}
  \item{xlab, ylab}{Labels for the x and y axis.}
  \item{\dots}{Other graphics parameters.}
}
\details{
  The clusters of exceedences are identified as follows.
  The first exceedence of the threshold initiates the first cluster.
  The first cluster then remains active until either \code{r}
  consecutive values fall below (or are equal to) the threshold,
  or until \code{rlow} consecutive values fall below (or are equal
  to) the lower threshold.
  The next exceedence of the threshold (if it exists) then initiates
  the second cluster, and so on.
  Missing values are allowed, in which case they are treated as
  falling below (or equal to) the threshold, but falling above the
  lower threshold.
}
\value{
  If \code{cmax} is \code{FALSE} (the default), a list with one
  component for each identified cluster.
  If \code{cmax} is \code{TRUE}, a numeric vector containing the
  cluster maxima.
  In any case, the returned object has an attribute \code{acs},
  giving the average cluster size (where the cluster size is
  defined as the number of exceedences within a cluster), which
  will be \code{NaN} if there are no values above the threshold
  (and hence no clusters).
  
  If \code{plot} is \code{TRUE}, the list of clusters, or vector
  of cluster maxima, is returned invisibly.
}

\seealso{\code{\link{exi}}, \code{\link{exiplot}}}

\examples{
clusters(portpirie, 4.2, 3)
clusters(portpirie, 4.2, 3, cmax = TRUE)
clusters(portpirie, 4.2, 3, 3.8, plot = TRUE)
clusters(portpirie, 4.2, 3, 3.8, plot = TRUE, lvals = FALSE)
tvu <- c(rep(4.2, 20), rep(4.1, 25), rep(4.2, 20))
clusters(portpirie, tvu, 3, plot = TRUE)
}
\keyword{manip}

