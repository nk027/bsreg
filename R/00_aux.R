
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


#' Density of a Beta Binomial distribution
#'
#' @param x Numeric scalar with the value.
#' @param size Integer scalar. Number of trials.
#' @param alpha,beta Numeric scalar. Shape parameters.
#' @param log Logical scalar. Whether to log probabilities.
#'
#' @noRd
dbbinom <- function(x, size, alpha, beta, log = FALSE) {
  out <- lchoose(size, x) + lbeta(x + alpha, size - x + beta) - lbeta(alpha, beta)
  if(isFALSE(log)) {out <- exp(out)}
  return(out)
}


#' Density of a Beta Negative Binomial distribution
#'
#' @param x Numeric scalar with the value.
#' @param size Integer scalar. Number of trials.
#' @param alpha,beta Numeric scalar. Shape parameters.
#' @param log Logical scalar. Whether to log probabilities.
#'
#' @noRd
dbnbinom <- function(x, size, alpha, beta, log = FALSE) {
  out <- lgamma(size + x) - lgamma(x + 1) - lgamma(size) + lbeta(alpha + size, beta + x) - lbeta(alpha, beta)
  if(isFALSE(log)) {out <- exp(out)}
  return(out)
}


#' p(tau | lambda) \propto (lambda - lambda^2)^tau tau^(alpha - 1) exp(-beta tau)
#'
#' @param x Numeric scalar with the value.
#' @param lambda Numeric scalar in (0, 1). Beta value to condition on.
#' @param alpha,beta Numeric scalar. Shape parameters.
#' @param log Logical scalar. Whether to log probabilities.
#'
#' @noRd
dtau <- function(x, lambda, alpha, beta, log = FALSE) {
  out <- x * log(lambda - lambda^2) + (alpha - 1) * log(x) - beta * x +
    alpha * log(beta) - lgamma(alpha) # Convenient for sampling
  if(isFALSE(log)) {out <- exp(out)}
  return(out)
}


#' Draw from p(tau | lambda)
#'
#' @param n Integer scalar with the number of draws.
#' @param lambda Numeric scalar in (0, 1). Beta value to condition on.
#' @param alpha,beta Numeric scalar. Shape parameters.
#'
#' @noRd
rtau <- \(n, lambda, alpha, beta, verbose = TRUE) {
  out <- numeric(n)
  i <- j <- 1
  # The Gamma density we draw from should be encompassing already
  m <- 1
  # exp(dtau(1e-32, lambda, alpha, beta, log = TRUE) - dgamma(1e-32, alpha, beta, log = TRUE))
  while (i <= n) {
    x <- rgamma(1L, alpha, beta)
    p_acc <- exp(dtau(x, lambda, alpha, beta, log = TRUE) - dgamma(x, alpha, beta, log = TRUE) - log(m))
    if(runif(1L) < p_acc) {
      out[i] <- x
      i <- i + 1
    }
    j <- j + 1
  }
  if(isTRUE(verbose)) cat("Tries: ", j, " (", round(j / n, 2), " draws per sample)", sep = "")
  return(out)
}
