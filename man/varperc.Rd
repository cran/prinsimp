\name{varperc}
\alias{varperc}

\title{Varperc method for Simple Partition}

\description{
  Produces the percent of variance explained plot for an object of the
  class \code{simpart}
}

\usage{
varperc(x, ...)
}

\arguments{
  \item{x}{
    an object of class \code{"simpart"}, typically result
    of \code{\link{simpart}}.
  }
  
  \item{...}{
    other parameters passed through to plotting functions
  }
}

\seealso{
  \code{\link{basisplot}}, \code{\link{varsimp}},
  \code{\link{simpart}}
}

\examples{
library(prinsimp)
require(graphics)

## Caterpillar estimated covariance from Kingsolver et al (2004)
## Measurements are at times 11, 17, 23, 29, 35, 40
data(caterpillar)

cat.sim <- simpart(caterpillar, simpledim = 2,
                 x = c(11, 17, 23, 29, 35, 40), cov = TRUE)

varperc(cat.sim)
}

\keyword{multivariate}
