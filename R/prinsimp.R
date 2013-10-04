### Performs a simplicity analysis on the principal components of the
### given covariance matrix and transforms them into a simple basis.
###
### covmat  - covariance matrix
### measure - simplicity measure to use ('first', 'second')
### x       - optional vector containing ordered values of the
###           environment levels (default: evenly spaced)
### reverse - (logical) if TRUE, changes direction of the
###           corresponding simple basis vector(s)
###
### Returns a list with class `prinsimp` containing the following
### components:
###   scores - simplicity scores of the basis vectors
###   loadings - the matrix whose columns contain the basis vectors
###   sdev -  standard deviations of the principal components
prinsimp <- function(covmat, measure = c('first', 'second'),
                     x=seq(d), reverse=FALSE) {
    cl <- match.call()
    cl[[1]] <- as.name('prinsimp')
    d <- nrow(covmat)
    if (!identical(d, ncol(covmat)))
        stop('Covariance matrix is not square')

    eig <- eigen(covmat, symmetric = TRUE)
    Gnnd <- eig$vectors %*% diag(ifelse(eig$values < 0,
                                        0, eig$values)) %*% t(eig$vectors)

    measure <- match.arg(measure)
    lambda_fun <- switch(measure,
                         first = lambda_first,
                         second = lambda_second)
    
    simple_basis <- simplify(eig$vectors, lambda_fun, x)

    simple_eigval <- diag(t(simple_basis) %*%
                          Gnnd %*% simple_basis)

    structure(list(scores=attr(simple_basis, 'simplicity'),
                   loadings=apply(simple_basis, 2, fix_direction) %*%
                            diag(ifelse(rep(reverse, length.out=d),
                                        -1, 1)),
                   sdev=sqrt(simple_eigval),
                   call=cl,
                   measure = measure),
              class='prinsimp')
}

## Choose the direction of the eigenvector so that the largest element
## is positive
fix_direction <- function(x) {
    if (length(x) > 0 && (max(x) != max(abs(x)))) -x else x
}


### Print methods for the prinsimp class
print.prinsimp <- function(x, ...) {
    cat(paste0("Simplicity scores (",
               if (!is.null(x$measure)) {
                   paste(x$measure, 'divided differences')
               } else 'unknown',
               "):\n"))
    print(x$scores, ...)
    cat("\nStandard deviations:\n")
    print(x$sdev, ...)
    cat("\n", nrow(x$loadings), " variables.\n")
}


summary.prinsimp <- function (object, loadings = FALSE, ...) {
    object$print.loadings <- loadings
    class(object) <- "summary.prinsimp"
    object
}


print.summary.prinsimp <- function (x, digits = 3, loadings = x$print.loadings, ...) {
    vars <- x$sdev^2
    vars <- vars/sum(vars)
    cat(paste0("Simplicity of basis (",
               if (!is.null(x$measure)) {
                   paste(x$measure, 'divided differences')
               } else 'unknown',
               "):\n"))
    print(rbind(`Standard deviation` = x$sdev,
                `Proportion of Variance` = vars,
                `Cumulative Proportion` = cumsum(vars)))
    if (loadings) {
        cat("\nLoadings:\n")
        cx <- format(round(x$loadings, digits = digits))
        print(cx, quote = FALSE, ...)
    }
    invisible(x)
}
