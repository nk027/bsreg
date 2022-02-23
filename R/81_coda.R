
#' Methods for \pkg{coda} Markov chain Monte Carlo objects
#'
#' Methods to convert parameter and/or coefficient draws to \pkg{coda}'s
#' \code{\link[coda]{mcmc}} format for further processing.
#'
#' @name coda
#'
#' @param x A \code{bm} object, obtained from \code{\link{bm}}.
#' @param ... Other parameters for \code{\link[coda]{as.mcmc}}.
#'
#' @return Returns a \pkg{coda} \code{\link[coda]{mcmc}} object.
NULL


#' @rdname coda
as.mcmc.bm <- function( # Dynamic export (zzz.R)
  x, ...) {

  # Checks ---

  if(!inherits(x, "bm")) {
    stop("Please provide a `bm` object.")
  }
  has_package("coda")

  out <- coda::as.mcmc(x[["draws"]], ...)

  return(out)
}
