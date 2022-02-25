
#' Fit a Bayesian model
#'
#' @param x Formula or \code{bm} object to sample with.
#' @param data A \code{\link{data.frame}} containing the variables in the model.
#' @param n_save,n_burn Integer scalar. Number of draws for the burn-in period and to store for inference.
#' @param type Character scalar with the type of prior setup.
#' @param options Settings for the prior setup. See \code{\link{set_options}}.
#' @param mh Settings to tune the Metropolis-Hastings step. See \code{\link{set_mh}}.
#' @param verbose Logical scalar. Whether to print status updates.
#' @param W Numeric matrix (or function to construct one) with the spatial connectivities.
#' @param X_SLX Numeric matrix with explanatory variables that should be lagged spatially.
#' @param type Character scalar used to specify the desired model.
#' @param ... Not used.
#'
#' @return Returns a list with draws from the specified Bayesian model and an object to obtain further samples.
#' @export
#'
#' @examples
#' N <- 100L
#' beta <- 1:5
#' X <- matrix(rnorm(N * 5), N, 5)
#' y <- X %*% beta + rnorm(N)
#'
#' bm(y ~ X, n_burn = 100, n_draw = 100)
#'
#' \donttest{
#' # Reproduce the linear model in Kuschnig (2022)
#' blm(log(sales) ~ log(price / cpi) + log(ndi / cpi) +
#'   factor(name) + factor(year), data = cigarettes)
#'}
bm <- function(x, ...) {UseMethod("bm", x)}


#' @export
#' @importFrom stats model.frame model.response model.matrix
#' @rdname bm
bm.formula <- function(x, data = NULL,
  n_save = 1000L, n_burn = 500L,
  options = set_options(), mh = set_mh(), verbose = TRUE,
  W, X_SLX,
  type = c("lm", "slx", "sar", "sem", "sdm", "sdem", "sv"),
  ...) {

  # Check inputs ---
  call <- match.call()
  type <- match.arg(type)
  getter <- switch(type, lm = get_blm, slx = get_bslx, sar = get_bsar, sem = get_bsem,
    sdm = get_bsdm, sdem = get_bsdem, sv = get_bsv)

  # Prepare data ---
  mf <- model.frame(x, data = data)
  y <- model.response(mf, "numeric")
  X <- model.matrix(attr(mf, "terms"), mf, contrasts = NULL)

  if(all(X[, 1] == X[, 2])) {X <- X[, -1]} # Drop double intercept
  if(type %in% c("slx", "sdm", "sdem") && missing(X_SLX)) {X_SLX <- X[, -1]} # Use all regressors except the intercept

  # Get model and estimate ---
  mdl <- getter(y = y, X = X, options = options, Psi = W, X_SLX = X_SLX, ...)

  draws <- sample(mdl, n_save = n_save, n_burn = n_burn, mh = mh, verbose = verbose)

  # Done ---
  return(structure(list("draws" = draws, "model" = mdl, "call" = call), class = "bm"))
}


#' @export
#' @rdname bm
bm.bm <- function(x, n_save = 1000L, n_burn = 0L, verbose = TRUE, ...) {

  draws <- rbind(x$draws, sample(x$model, n_save = n_save, n_burn = n_burn, verbose = verbose))

  # Done ---
  return(structure(list("draws" = draws, "model" = x$model, "call" = x$call), class = "bm"))

}


#' @export
#' @rdname bm
blm <- function(...) {bm(..., type = "lm")}
#' @export
#' @rdname bm
bslx <- function(...) {bm(..., type = "slx")}
#' @export
#' @rdname bm
bsar <- function(...) {bm(..., type = "sar")}
#' @export
#' @rdname bm
bsem <- function(...) {bm(..., type = "sem")}
#' @export
#' @rdname bm
bsdm <- function(...) {bm(..., type = "sdm")}
#' @export
#' @rdname bm
bsdem <- function(...) {bm(..., type = "sdem")}
#' @export
#' @rdname bm
bsv <- function(...) {bm(..., type = "sv")}
