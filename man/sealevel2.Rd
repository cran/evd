\name{sealevel2}
\alias{sealevel2}
\title{Annual Sea Level Maxima at Dover and Harwich with Indicator}
\usage{data(sealevel2)}
\description{
  The \code{sealevel2} data frame has 81 rows and 3 columns.
  The first two columns contain annual sea level maxima from 1912
  to 1992 at Dover and Harwich respectively, two sites on the coast
  of Britain.
  The third column is a logical vector denoting whether or not the
  maxima in a given year are assumed to have derived from the
  same storm event; this assumption is made if the times of
  obsevation of the maxima are at most 48 hours apart.
  The row names give the years of observation.
  There are 39 missing data values.
  There are only nine non-missing logical values.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{dover}{A numeric vector containing annual sea level maxima
    at Dover, including 9 missing values.}
    \item{harwich}{A numeric vector containing sea annual level maxima
      at Harwich, including 30 missing values.}
    \item{case}{A logical vector denoting whether or not the maxima
      are assumed to have derived from the same storm event.}
  }
}
\source{
  Coles, S. G. and Tawn, J. A. (1990)
  Statistics of coastal flood prevention.
  \emph{Phil. Trans. R. Soc. Lond., A} \bold{332}, 457--476.
}
\keyword{datasets}

