\name{summary.prinsimp}
\alias{summary.prinsimp}
\alias{print.summary.prinsimp}

\title{Summary method for Simplicity basis}

\description{
  The \code{\link{summary}} method for class \code{"prinsimp"}. Prints
  the percent of total variance explained by each basis vector and the
  cumulative percent of total variance explained. Optionally, the basis
  vectors can also be printed.
}

\usage{
\method{summary}{prinsimp}(object, loadings = FALSE, \dots)

\method{print}{summary.prinsimp}(x, digits = 3, loadings = x$print.loadings,
       \dots)
}

\arguments{
  \item{object}{
    an object of class \code{"prinsimp"}, as
    from \code{prinsimp()}.
  }
  
  \item{loadings}{
    logical. If true, all basis vectors are printed.
  }
  
  \item{x}{
    an object of class \code{summary.prinsimp}.
  }
  
  \item{digits}{
    the number of significant digits to be used in listing
    loadings.
  }
  
  \item{\dots}{
    arguments to be passed to or from other methods.
  }
}

\value{
  \code{object} with additional component \code{print.loadings}.
}

\seealso{
  \code{\link{prinsimp}}
}

\examples{
library(prinsimp)

## Caterpillar data: estimated covariance from Kingsolver et al (2004)
## Measurements are at temperatures 11, 17, 23, 29, 35, 40
data(caterpillar)

cat.sim <- prinsimp(caterpillar, x = c(11, 17, 23, 29, 35, 40))
summary(cat.sim)

print(summary(cat.sim, loadings = TRUE), digits = 2)
}

\keyword{multivariate}
