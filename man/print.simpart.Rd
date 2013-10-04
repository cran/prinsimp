\name{print.simpart}
\alias{print.simpart}

\title{Print method for Simple Partition}

\description{
  The \code{\link{print}} method for class \code{"simpart"}. Prints the
  simplicity values for both the model basis vectors and the vectors in
  the simplicity basis of the nearly null space, the percent of variance
  explained by each basis vector, and the cummulative percent of
  variance explained up to each basis vector, relative to the total
  variance explained by the subspace (model or nearly null).
}

\usage{
\method{print}{simpart}(x, \dots)
}

\arguments{
  \item{x}{
    an object of class \code{"simpart"}, typically result
    of \code{simpart()}.
  }
  
  \item{\dots}{
    arguments to be passed to or from other methods.
  }
}

\seealso{
  \code{\link{simpart}}
}

\examples{
library(prinsimp)

## Caterpillar data: estimated covariance from Kingsolver et al (2004)
## Measurements are at temperatures 11, 17, 23, 29, 35, 40
data(caterpillar)

cat.sim <- simpart(caterpillar, simpledim = 2,
                   x = c(11, 17, 23, 29, 35, 40), cov = TRUE)
print(cat.sim)
}

\keyword{multivariate}
