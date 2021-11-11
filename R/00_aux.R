
#' Check numeric scalar
#'
#' Check whether an object is bounded and coercible to a numeric value.
#'
#' @param x Numeric scalar.
#' @param min Numeric scalar. Minimum value of \emph{x}.
#' @param max Numeric scalar. Maximum value of \emph{x}.
#' @param fun Function to apply to \emph{x} before returning.
#' @param msg String fed to \code{\link[base]{stop}} if an error occurs.
#'
#' @return Returns \code{fun(x)}.
#'
#' @noRd
num_check <- function(
  x, min = 0, max = Inf,
  msg = "Please check the numeric parameters.",
  fun = as.numeric) {

  if(!is.numeric(x) || length(x) != 1 || x < min || x > max) {stop(msg)}

  return(fun(x))
}


#' @noRd
int_check <- function(
  x, min = 0L, max = Inf,
  msg = "Please check the integer parameters.") {

  num_check(x, min, max, msg, fun = as.integer)
}

#' @noRd
num_default <- function(x, default, min, max, msg) {
  if(missing(x)) {return(default)}
  num_check(x, min, max, msg)
}


#' Multivariate Normal
#'
#' Draw from a multivariate Normal using the precision instead of variance.
#'
#' @param n Integer scalar. Number of draws.
#' @param mu Numeric vector.
#' @param precision Numeric matrix.
#'
#' @return Returns a matrix with \emph{n} rows of draws.
#'
#' @importFrom stats rnorm
#'
#' @noRd
rmvn <- function(n, mu, precision) {

  # Spectral  ---
  # ev <- eigen(precision, symmetric = TRUE)
  # m <- length(ev[["values"]])
  # R <- t(ev[["vectors"]] %*% (t(ev[["vectors"]]) * sqrt(1 / pmax(ev[["values"]], 0))))
  # out <- matrix(rnorm(n * m), nrow = n, ncol = m, byrow = TRUE) %*% R

  # Cholesky ---
  m <- ncol(precision)
  R <- chol(precision)
  out <- t(backsolve(R, matrix(rnorm(n * m), nrow = m, ncol = n, byrow = TRUE)))

  if(!missing(mu)) {out <- sweep(out, 2, mu, "+")}

  return(out)
}


#' Check Sparsity
#'
#' @param x Matrix.
#'
#' @return Returns a logical scalar indicating sparsity.
#'
#' @noRd
is_sparse <- function(x) {
  if(inherits(x, "function")) {return(isTRUE(attr(x, "sparse")))} # Allow setting via attribute
  isTRUE(inherits(x, "dgCMatrix"))
}

#' Check Symmetry
#'
#' @param x A numeric matrix.
#'
#' @return Returns a logical scalar indicating symmetry.
#'
#' @noRd
is_symmetric <- function(x) {
  if(inherits(x, "function")) {return(isTRUE(attr(x, "symmetric")))} # Allow setting via attribute
  isSymmetric(x)
}


#' Construct integer interval from vector
#'
#' @param x A numeric vector with the start, end, and length of the interval.
#'
#' @return Returns an integer interval.
#'
#' @noRd
i_seq <- function(x) {
  seq.int(x[1], x[2], length.out = x[3])
}


#' Sum of squares
#'
#' @param ... Numeric vectors or scalars that are squared and summed.
#'
#' @return Returns an numeric scalar.
#'
#' @noRd
sq_sum <- function(...) {
  sum((...)^2)
}


#' Check whether a package is installed
#'
#' @param package Character scalar.
#'
#' @noRd
has_package <- function(package) {

  if(!requireNamespace(package, quietly = TRUE)) {
    stop("Package \'", package, "\' required for this method.", call. = FALSE)
  }

  return(NULL)
}
