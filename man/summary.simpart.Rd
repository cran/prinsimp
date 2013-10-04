\name{summary.simpart}
\alias{summary.simpart}
\alias{print.summary.simpart}

\title{Summary method for Simple Partition}

\description{
  The \code{\link{summary}} method for class \code{"simpart"}. Prints
  the dimensions of the nearly null space in the \code{simpart} object,
  the percent of total variance explained by each basis vector, the
  cumulative percent of total variance explained (with accumulation
  restarting in the nearly null space), and the simplicity values of
  each basis vector. The model basis vectors are ordered by the percent
  of variance explained in descending order. The vectors in the
  simplicity basis of the nearly null space are ordered by their
  simplicity measures in descending order.
}

\usage{
\method{summary}{simpart}(object, loadings = FALSE, \dots)

\method{print}{summary.simpart}(x, digits = 3, loadings = x$print.loadings,
       \dots)
}

\arguments{
  \item{object}{
    an object of class \code{"simpart"}, as
    from \code{simpart()}.
  }
  
  \item{loadings}{
    logical. If true, all basis vectors are printed.
  }
  
  \item{x}{
    an object of class "summary.simpart".
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
  \code{\link{simpart}}
}

\examples{
library(prinsimp)

## Caterpillar data: estimated covariance from Kingsolver et al (2004)
## Measurements are at temperatures 11, 17, 23, 29, 35, 40
data(caterpillar)

cat.sim <- simpart(caterpillar, simpledim = 2,
                   x = c(11, 17, 23, 29, 35, 40), cov = TRUE)
summary(cat.sim)

print(summary(cat.sim, loadings = TRUE), digits = 2)
}

\keyword{multivariate}
