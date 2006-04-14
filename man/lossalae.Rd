\name{lossalae}
\alias{lossalae}
\title{General Liability Claims}
\usage{data(lossalae)}
\description{
  The \code{lossalae} data frame has 1500 rows and 2 columns.
  The columns contain the indemnity payment (loss), and
  the allocated loss adjustment expense (alae), both in USD.
  The latter is the additional expenses associated with the
  settlement of the claim (e.g. claims investigation expenses
  and legal fees).

  The dataset also has an attribute called \code{capped}, which
  gives the row names of the indemnity payments that were capped
  at their policy limit. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Loss}{A numeric vector containing the indemnity
      payments.}
    \item{ALAE}{A numeric vector containing the allocated
      loss adjustment expenses.}
  }
}
\source{
  Frees, E. W. and Valdez, E. A. (1998)
  Understanding relationships using copulas.
  \emph{North American Actuarial Journal}, \bold{2}, 1--15.
}
\references{
  Klugman, S. A. and Parsa, R. (1999)
  Fitting bivariate loss distributions with copulas.
  \emph{Insurance: Mathematics and Economics}, \bold{24}, 139--148.
  
  Beirlant, J., Goegebeur, Y., Segers, J. and Teugels, J. L. (2004)
  \emph{Statistics of Extremes: Theory and Applications.},
  Chichester, England: John Wiley and Sons.
}
\keyword{datasets}
