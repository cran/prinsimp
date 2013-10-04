\name{prinsimp}
\alias{prinsimp}

\title{
  Simplicity basis
}

\description{
  \code{prinsimp} calculates a simplicity basis of a given covariance
  matrix, using a specified quadratic simplicity measure.
}

\usage{
prinsimp(covmat, measure = c('first', 'second'),
         x = seq(d), reverse=FALSE)
}

\arguments{
  \item{covmat}{
    a \eqn{d \times d}{d x d} covariance matrix
  }
  
  \item{measure}{
    a function that calculates a simplicity measure of a vector, based
    on a non-negative definite symmetric matrix Lambda. There are two
    available simplicity measures, specified as 'first' or 'second',
    that correspond to first or second divided difference, respectively.
  }
  
  \item{x}{
    a vector of independent variable values (for functional data),
    length equal to \eqn{d}, the number of columns of \code{covmat}. If
    not supplied, a sequence from 1 to \eqn{d} is used.
  }
  
  \item{reverse}{
    a logical vector of length d. If the i-th element is true, the i-th
    basis vector is "reversed" by multiplication by -1.  Basis vectors
    are arranged with model basis first, then simplicity basis. If
    length of \code{reverse} is less than d, then the remaining entries
    of \code{reverse} are assumed to be false, and the corresponding
    basis vectors remain unchanged.
  }
}

\value{
  \code{prinsimp} returns a list with class \code{"prinsimp"} containing
  the following components:

  \item{loadings}{
    \eqn{d \times d}{d x d} matrix with columns containing the
    simplicity basis of the nearly null space. Basis vectors are
    arranged in descending order of simplicity.
  }
  
  \item{sdev}{
    standard deviations associated with the vectors in the simplicity
    basis of the nearly null space.
  }
  
  \item{scores}{
    eigenvalues of the vectors in the simplicity basis of the nearly
    null space.
  }
  
  \item{call}{
    the matched call
  }
  
  \item{measure}{
    the simplicity measure used: \code{"first"} or \code{"second"}
  }
}

\references{
  T.L. Gaydos, N.E. Heckman, M. Kirkpatrick, J.R. Stinchcombe, J.
  Schmitt, J. Kingsolver, J.S. Marron. (2013). Visualizing genetic
  constraints. \emph{Annals of Applied Statistics} 7: 860-882.
}

\seealso{
  \code{\link{summary.prinsimp}}, \code{\link{simpart}}
}

\examples{
library(prinsimp)

## Caterpillar data: estimated covariance from Kingsolver et al (2004)
## Measurements are at temperatures 11, 17, 23, 29, 35, 40
data(caterpillar)

prinsimp(caterpillar, x = c(11, 17, 23, 29, 35, 40))
}

\keyword{multivariate}
