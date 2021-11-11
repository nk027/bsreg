
#' @export
print.bm <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Draws (total): ", NROW(x$draw), " (", x$model$get_meta, ")\n", sep = "")
  coefs <- colMeans(x$draw)
  cat("Coefficients (mean):\n")
  print.default(format(coefs, digits = 3L), print.gap = 2L, quote = FALSE)
  cat("\n")
  # invisible(x)
}


#' @export
summary.bm <- function(x, ...) {
  summary(x$draws)
}


#' @export
plot.bm <- function(x, ...) {
  plot.ts(x$draws)
}
