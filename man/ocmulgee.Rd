\name{ocmulgee}
\alias{ocmulgee}
\title{Maximum Annual Flood Discharges of the Ocmulgee River}
\usage{data(ocmulgee)}
\description{
  The \code{ocmulgee} data frame has 40 rows and 2 columns.
  The columns contain maximum annual flood discharges, in units
  of 1000 cubed feet per second, from the Ocmulgee River in
  Georgia, USA at Hawkinsville (upstream) and Macon (downstream),
  for the years 1910 to 1949.
  The row names give the years of observation.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{hawk}{A numeric vector containing maximum annual flood
      discharges at Hawkinsville (upstream).}
    \item{macon}{A numeric vector containing maximum annual flood
      discharges at Macon (downstream).}
  }
}
\source{
  Gumbel, E. J. and Goldstein, N. (1964)
  Analysis of empirical bivariate extremal distributions.
  \emph{J. Amer. Statist. Assoc.}, \bold{59}, 794--816.
}
\keyword{datasets}
