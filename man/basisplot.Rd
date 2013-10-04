\name{basisplot}
\alias{basisplot}
\title{Basisplot method for Simple Partition}

\description{
  Produces the plots for the specified basis vectors of the class
  \code{simpart}
}

\usage{
basisplot(x, display = list(model=TRUE, simple=TRUE), ...)
}

\arguments{
  \item{x}{
    an object of class \code{"simpart"}, typically result
    of \code{\link{simpart}}.
  }
  
  \item{display}{
    a list. Specifies which model and nearly null space simplicity basis
    vectors to print.
  }
  
  \item{...}{
    other parameters passed through to plotting functions
  }
}

\seealso{
  \code{\link{varsimp}}, \code{\link{varperc}},
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

## Plots the 4 model basis vectors and 2 simplicity basis vectors consecutively
basisplot(cat.sim)
    
## Display the 4 model basis vectors on one page
par(mfrow = c(2,2))
basisplot(cat.sim, display = list(model=1:4))
}

\keyword{multivariate}
