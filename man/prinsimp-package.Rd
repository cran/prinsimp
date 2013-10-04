\docType{package}
\name{prinsimp-package}
\alias{plyr-package}

\title{plyr: the split-apply-combine paradigm for R.}

\description{
  The plyr package is a set of clean and consistent tools
  that implement the split-apply-combine pattern in R. This
  is an extremely common pattern in data analysis: you
  solve a complex problem by breaking it down into small
  pieces, doing something to each piece and then combining
  the results back together again.
}

\details{
  Principal Component Analysis (PCA) is simply an eigenanalysis of a
  covariance matrix, \eqn{G}. Its eigenvalues \eqn{\lambda_j} can be
  interpreted as the variance of \eqn{G} in the direction of the
  eigenvector \eqn{v_j}, and \eqn{\lambda_j / \sum \lambda_k} as the
  proportion of variance explained by \eqn{v_j}. Often, \eqn{G} is
  well-approximated using the first \eqn{M} eigenvectors and
  eigenvalues, which we call the \emph{model space}.

  The orthogonal complement of the model space (that is, the space
  spanned by the remaining eigenvectors), we call the \emph{nearly null}
  space. The nearly null space is interesting as a space of low
  variability, which may be particularly important in, for instance,
  evolutionary biology. This package provides functions for analyzing
  the nearly null space and finding interesting structures of low
  variability, as defined by a quadratic simplicity measure. It is an
  expanded reimplementation in R of the method described by Gaydos et
  al (2013).
}

\references{
  T.L. Gaydos, N.E. Heckman, M. Kirkpatrick, J.R. Stinchcombe, J.
  Schmitt, J. Kingsolver, J.S. Marron. (2013). Visualizing genetic
  constraints. \emph{Annals of Applied Statistics} 7: 860-882.
}
