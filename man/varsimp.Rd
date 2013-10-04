\name{varsimp}
\alias{varsimp}

\title{Varsimp method for Simple Partition}

\description{
  Produces the variance-simplicity view for the specified basis vectors
  of the d-dimensional sample space.
}

\usage{
varsimp(x, display = list(model=TRUE, simple=TRUE), full.simple = TRUE, ...)
}

\arguments{
  \item{x}{
    an object of class \code{"simpart"}, typically result
    of \code{\link{simpart}}.
  }
  
  \item{display}{
    a list specifying which model basis vectors and which nearly null
    space simplicity basis vectors to print.
  }
  
  \item{full.simple}{
    logical indicating whether to draw the background
    lines at simplicity values of the full space.
  }
  \item{...}{
    other parameters passed through to plotting functions
  }
}

\seealso{
  \code{\link{basisplot}}, \code{\link{varperc}},
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

## Display all 6 basis vectors in the variance-simplicity view
varsimp(cat.sim)
    
## Display only the 4 model basis vectors in the variance-simplicity view
varsimp(cat.sim, display = list(model=1:4))
}

\keyword{multivariate}
