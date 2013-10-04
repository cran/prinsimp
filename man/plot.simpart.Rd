\name{plot.simpart}
\alias{plot.simpart}

\title{Plot method for Simple Partition}

\description{
  Produces plots for the model basis vectors and the vectors forming the
  simplicity basis for the nearly null space, along with the
  variance-simplicity view and the percent of variance explained panel.
  Only six basis vectors can be displayed at a time.
}

\usage{
\method{plot}{simpart}(x, display = list(model=TRUE, simple=TRUE), layout, ...)
}

\arguments{
  \item{x}{
    an object of class \code{"simpart"}, typically result
    of \code{\link{simpart}}.
  }
  \item{display}{
    a list. Specifies which model and simplicity basis vectors to plot.
    Must be consistent with \code{x}.
  }
  \item{layout}{
    a matrix of locations for the component plots, following the format
    used in the \code{\link[graphics]{layout}} function. The default is
    to draw basis plots in the left half of the figure area, arranged in
    a U-curve with model basis starting at the top of the left arm and
    simple basis at the top right. The variance-simplicity view and
    percent of variance explained panel take up the right half of the
    figure.

    The subplots are drawn in the following order: first model then
    simple basis, followed by the variance-simplicity plot, and
    variance-explained last.
  }
  
  \item{...}{
    other parameters passed through to plotting functions
  }
}

\details{
  The variance-simplicity view and the percent of variance explained
  panel are always produced everytime \code{plot.simpart} is called. The
  variance-simplicity view plots the percent of variance explained and
  the simplicity values of the basis vectors specified in
  \code{display}. If the number of basis vectors is greater than 6 and
  the user specifies more than 6 in \code{display}, an error message is
  generated. In \code{display}, non-existent basis numbers are ignored.
}

\seealso{
  \code{\link{basisplot}}, \code{\link{varsimp}}, \code{\link{varperc}},
  \code{\link{simpart}}
}

\examples{
library(prinsimp)
require(graphics)

## Caterpillar estimated covariance matrix from Kingsolver et al (2004)
## Measurements are at temperatures 11, 17, 23, 29, 35, 40
data(caterpillar)

## Default plot method for a 2-dimensional nearly null space
cat.sim <- simpart(caterpillar, simpledim = 2,
                   x = c(11, 17, 23, 29, 35, 40), cov = TRUE)
plot(cat.sim)

## Plot only the two simplicity basis vectors in the 2-dimensional
## nearly null space
plot(cat.sim, display = list(simple=1:2))
    
## Try, by mistake to plot first four simple vectors when simpledim=2.
## Will plot only 2.
plot(cat.sim, display = list(simple=1:4))

## Only plot the first two model basis vectors, the variance simplicity
## view, and the percent variance explained panel using the 'layout'
## argument
plot(cat.sim, display = list(model=1:2), layout = matrix(1:4, nrow=2, ncol=2))
}

\keyword{multivariate}
