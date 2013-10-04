\name{caterpillar}
\alias{caterpillar}
\docType{data}

\title{
  Kingsolver \emph{et al.} caterpillar estimated covariance
}

\description{
  Estimated broad-sense genetic variance-covariance matrix for
  short-term growth rate at different temperatures in fourth-instar
  caterpillars from Kingsolver et al (2004).
  
  Measurements are at temperatures 11, 17, 23, 29, 35, 40.
}

\usage{
data(caterpillar)
}

\format{
  A 6-by-6 estimated genetic covariance matrix. The matrix is obtained
  from the printed version of the paper and, due to rounding, it is not
  positive definite -- one of its eigenvalues is negative. The functions
  in this package, such as \code{simpart} handle this by setting all
  negative eigenvalues to zero and reconstructing the covariance matrix
  before proceeding.
}

\source{
  Kingsolver, J. G., Ragland, G. J. and Shlichta, J. G. (2004).
  Quantitative genetics of continuous reaction norms: thermal
  sensitivity of caterpillar growth rates. \emph{Evolution},
  58:1521--1529.
}

\examples{
library(prinsimp)

data(caterpillar)
cat.sim <- simpart(caterpillar, simpledim = 2,
                   x = c(11, 17, 23, 29, 35, 40), cov = TRUE)
plot(cat.sim)
}

\keyword{datasets}
