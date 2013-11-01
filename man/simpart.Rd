\name{simpart}
\alias{simpart}
\alias{simpart.default}
\alias{simpart.formula}

\title{
  Simple Partition
}
\description{
  \code{simpart} partitions a \eqn{d}-dimensional sample space into two
  orthonormal subspaces: a \code{simpledim}-dimensional nearly null
  space and a \eqn{(d-\mathtt{simpledim})}{(d-simpledim)}-dimensional
  model space. It provides an orthonormal basis for each subspace. The
  nearly null space basis is defined in terms of a simplicity measure
  and is ordered from most simple to least simple. The model space basis
  is made up of leading eigenvectors of the covariance matrix and is
  ordered by proportion of variance explained.

  Returns the result as an object of class \code{simpart}.
}

\usage{
simpart(y, simpledim, ...)

\method{simpart}{formula}(formula, simpledim, data = NULL, \dots)

\method{simpart}{default}(y, simpledim, measure = c('first', 'second', 'periodic'),
        x = seq(d), cov=FALSE, reverse=rep(FALSE, d), na.action, ...)
}

\arguments{
  \item{formula}{
    a formula with no response variable, referring only to numeric
    variables.
  }
  
  \item{y}{
    a matrix or data frame that specifies the data, or a covariance
    matrix. Data matrix has d columns, covariance matrix is \eqn{d
    \times d}{d x d}.
  }
  
  \item{simpledim}{
    the dimension of the  nearly null space of the covariance matrix. It
    is equal to \eqn{d} minus the dimension of the model space.
  }

  \item{measure}{
    a function that calculates a simplicity measure of a vector, based
    on a non-negative definite symmetric matrix Lambda. There are three
    built in simplicity measures, specified by 'first', 'second', or
    'periodic' that correspond to first divided difference, second
    divided difference and periodic simplicity respectively. The
    argument \code{measure} can take a user specified function.
  }
  
  \item{data}{
    an optional data frame (or similar: see \code{\link{model.frame}})
    containing the variables in the formula \code{formula}. By default
    the variables are taken from \code{environment(formula)}.
  }
  
  \item{x}{
    a vector of independent variable values (for functional data),
    length equal to \eqn{d}, the number of columns of \code{y}. If not
    supplied, a sequence from 1 to \eqn{d} is used.
  }
  
  \item{cov}{
    a logical value. If true, then \code{y} is assumed to be a \eqn{d
    \times d}{d x d} covariance matrix. If false, y is assumed to be an
    \eqn{n \times d}{n x d} data matrix which \code{simpart} uses to
    calculate a \eqn{d \times d}{d x d} covariance matrix.
  }
  
  \item{reverse}{
    a logical vector of length d. If the i-th element is true, the i-th
    basis vector is "reversed" by multiplication by -1.  Basis vectors
    are arranged with model basis first, then simplicity basis. If
    length of \code{reverse} is less than d, then the remaining entries
    of \code{reverse} are assumed to be false, and the corresponding
    basis vectors remain unchanged.
  }
  
  \item{na.action}{
    specify how missing data should be treated.
  }
  
  \item{\dots}{
    arguments passed to or from other methods. If \code{x} is a formula
    one might specify \code{cov} or \code{reverse}. If \code{"periodic"}
    is chosen as the \code{measure}, period is specified as a numeric.
    If \code{measure} is user specified, its arguments are passed here.
  }
}

\details{
  \code{simpart} is a generic function with \code{"formula"} and
  \code{"default"} methods.

  \code{simpart} implements a method described in Gaydos et al (2013).

  When \code{cov=FALSE}, the covariance matrix is calculated using the
  data matrix \code{y}. The calculation uses divisor \eqn{n}, the number
  of rows of \code{y}.
}

\value{
  \code{simpart} returns a list with class \code{"simpart"} containing
  the following components:

  \item{model}{
    a \eqn{d \times (d-\mathtt{simpledim})}{d x (d-simpledim)} matrix
    with columns containing the basis of the model space, that is,
    containing the first \eqn{(d-\mathtt{simpledim})}{(d-simpledim)}
    eigenvectors of the covariance matrix. Basis vectors are arranged in
    descending order of eigenvalue, that is, in descending order of the
    proportion of variance explained.
  }

  \item{simple}{
    \eqn{d \times \code{simpledim}}{d x simpledim} matrix with columns
    containing the simplicity basis of the nearly null space. Basis
    vectors are arranged in descending order of simplicity.
  }
  
  \item{variance}{
    list of three components:
    \describe{
      \item{model}{variances associated with the vectors in the model
        basis.}

      \item{simple}{variances associated with the vectors in the
        simplicity basis of the nearly null space.}

      \item{full}{variances associated with eigenvectors of the
        covariance matrix, that is, its eigenvalues.}
    }
  }
  
  \item{simplicity}{
    list of three components:
    \describe{
      \item{model}{simplicity values of the vectors in the model basis.}

      \item{simple}{eigenvalues of the vectors in the simplicity basis
        of the nearly null space.}

      \item{full}{simplicity values of the simplicity basis when
        \code{simpledim=d}.}
    }
  }
  
  \item{call}{
    the matched call
  }
  
  \item{measure}{
    the simplicity measure used: \code{"first"}, \code{"second"},
    \code{"periodic"} or an user specified measure function
  }
  
  \item{varperc}{
    the percent of variance explained by the corresponding basis vector,
    as a list of two components:
    \describe{
      \item{model}{percent of variance explained by the vectors in the
        model basis.}

      \item{simple}{percent of variance explained by the vectors in the
        simplicity basis of the nearly null space.}
    }
  }
  
  \item{scores}{
    if \code{y} is the data matrix, the scores on the basis vector
    loadings.
  }
}

\references{
  T.L. Gaydos, N.E. Heckman, M. Kirkpatrick, J.R. Stinchcombe, J.
  Schmitt, J. Kingsolver, J.S. Marron. (2013). Visualizing genetic
  constraints. \emph{Annals of Applied Statistics} 7: 860-882.
}

\note{  
  The simplicity values of the simplicity basis when \code{simpledim=d}
  are equal to the eigenvalues of the non-negative definite matrix,
  Lambda, that defines the simplicity measure.
}

\seealso{
  \code{\link{summary.simpart}}, \code{\link{plot.simpart}}
}

\examples{
library(prinsimp)
require(graphics)

## Caterpillar data: estimated covariance from Kingsolver et al (2004)
## Measurements are at temperatures 11, 17, 23, 29, 35, 40
data(caterpillar)

## Analyze 5 dimensional model space, 1 dimensional nearly null space
## First divided difference simplicity measure
simpart(caterpillar, simpledim=1, cov=TRUE)  # Need to specify x

simpart(caterpillar, simpledim=1,
        x=c(11, 17, 23, 29, 35, 40), cov=TRUE)

## Second divided difference simplicity measure and 3-dimensional model space
simpart(caterpillar, simpledim=3, measure="second",
        x=c(11, 17, 23, 29, 35, 40), cov=TRUE)
}

\keyword{multivariate}
