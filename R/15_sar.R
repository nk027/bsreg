
#' Bayesian spatial autoregressive model
#'
#' @docType class
#'
#' @param parent \code{\link{R6Class}} object to inherit from.
#'
#' @importFrom R6 R6Class
#' @importFrom stats splinefun
#'
#' @noRd
get_sar_class <- function(parent = NormalGamma) {

  SpatialAR <- R6Class("SpatialAR", inherit = parent,

  public = list(

    initialize_SAR = function(priors, ...) {

      # Store prior settings ---

      if(missing(priors) || is.null(priors$SAR)) {priors <- list(SAR = set_SAR())}
      private$SAR$priors <- priors$SAR

      # Build SAR object ---

      # Function to update the latent, i.e. z(lambda) = y - lambda W y, with new parameters
      private$SAR$set_latent <- function(lambda = private$SAR$lambda, Wy = self$Wy) {
        private$SAR$z <- super$y - lambda * Wy
      }

      # Function to set a new lambda and update the latent
      private$SAR$set_lambda <- function(lambda) {
        private$SAR$lambda <- lambda
        # Update the latent
        private$SAR$set_latent(lambda = lambda)
      }

      # Initialise MH ---

      self$MH$SAR_lambda <- MH_SAR_lambda$new(value = priors$SAR$lambda, scale = priors$SAR$lambda_scale,
        shape_a = priors$SAR$lambda_a, shape_b = priors$SAR$lambda_b,
        lower = priors$SAR$lambda_min, upper = priors$SAR$lambda_max)
    },

    # High priority 9 for latent
    setup_9SAR = function(Psi_SAR = NULL, ldet_SAR = list(grid = FALSE, i_lambda = c(-1, 1 - 1e-12, 100L), reps = 1L),
      ...) {

      # Work out connectivity ---

      if(is.null(Psi_SAR)) {stop("Please provide a connectivity matrix 'Psi_SAR'.")}

      private$SAR$Psi_fixed <- TRUE
      private$SAR$W <- private$SAR$Psi <- Psi_SAR
      # # Set cache
      private$SAR$Wy <- private$SAR$W %*% super$y
      # Set lambda and obtain the latent
      private$SAR$set_lambda(private$SAR$priors$lambda)

      # Set up MH ---

      self$MH$SAR_lambda$setup(N = private$cache$N, M = private$cache$M)

      # Set up the log determinant ---

      # Initialise object with options
      private$SAR$ldet <- list(size = private$cache$N / ldet_SAR$reps, reps = ldet_SAR$reps, # Kronecker settings
        grid = isTRUE(ldet_SAR$grid), i_lambda = ldet_SAR$i_lambda) # Grid settings

      # If W is repeated via Kronecker product we only need one submatrix and scale up using `reps` later
      private$SAR$ldet$get_W <- if(private$SAR$ldet$reps == 1) {
        function(W = private$SAR$W) {W}
      } else { # Just retrieve the sub-matrix
        function(W = private$SAR$W) {W[seq(private$SAR$ldet$size), seq(private$SAR$ldet$size)]}
      }

      # Provide a function for the log-determinant
      if(private$SAR$ldet$grid) { # If a grid is requested we fit a spline to a grid over lambda

        pars <- i_seq(private$SAR$ldet$i_lambda)
        ldets <- vapply(pars, function(x) {
          determinant(diag(private$SAR$ldet$size) - x * private$SAR$ldet$get_W(),
          logarithm = TRUE)$modulus * private$SAR$ldet$reps
        }, numeric(1L))
        private$SAR$ldet$splinefun <- splinefun(pars, y = ldets)
        private$SAR$ldet$get_ldet <- function(lambda, ...) {private$SAR$ldet$splinefun(x = lambda)}

      } else { # Otherwise we use a spectral decomposition

        private$SAR$ldet$ev <- eigen(private$SAR$ldet$get_W(),
          symmetric = is_symmetric(private$SAR$Psi), only.values = TRUE)$values
        # The log-determinant of I - lambda W is just the sum of log(1 - lambda * omega)
        private$SAR$ldet$get_ldet <- function(lambda, ...) {
          Re(sum(log(1 - lambda * private$SAR$ldet$ev))) * private$SAR$ldet$reps
        }
      }
    },

    # Sample spatial autoregressive parameter lambda ---
    sample_latent = function() {

      # Prepare RSS and the log-determinant as functions of lambda
      get_rss <- {function() {
        prec_ch <- chol(private$NG$prec0 + self$XX / self$sigma)
        b0 <- backsolve(prec_ch, forwardsolve(prec_ch, (private$NG$prec0 %*% private$NG$mu0 + crossprod(self$X, super$y) / self$sigma),
          upper.tri = TRUE, transpose = TRUE))
        b1 <- backsolve(prec_ch, forwardsolve(prec_ch, (private$NG$prec0 %*% private$NG$mu0 + self$XWy / self$sigma),
          upper.tri = TRUE, transpose = TRUE))
        e0 <- super$y - self$X %*% b0
        e1 <- self$Wy - self$X %*% b1
        e0e0 <- sum(e0^2)
        e1e0 <- sum(e1 * e0)
        e1e1 <- sum(e1^2)
        return(function(value) {(e0e0) - (2 * value * e1e0) + (value^2 * e1e1)})
      }}()
      get_ldet <- function(value) {private$SAR$ldet$get_ldet(lambda = value)}

      # Metropolis-Hastings step for lambda
      self$MH$SAR_lambda$propose()
      self$MH$SAR_lambda$acceptance(get_rss = get_rss, get_ldet = get_ldet)
      self$MH$SAR_lambda$finalize()
      lambda <- self$MH$SAR_lambda$get_value # Assign and recompute
      if(abs(lambda - private$SAR$lambda) > 1e-12) {private$SAR$set_lambda(lambda)}
    }

  ),

  active = list(

    # Variables that are adapted using the latent ---
    y = function() {private$SAR$z},
    Xy = function() {crossprod(self$X, private$SAR$z)},
    Wy = function() {private$SAR$Wy},
    XWy = function() {crossprod(self$X, private$SAR$Wy)},

    # Access functions ---
    get_parameters = function() {
      pars <- super$get_parameters
      pars$lambda_SAR <- private$SAR$lambda
      return(pars)
    },
    get_effects = function() { # To-do: use eigendecomposition or provide alternative methods to be more efficient
      total <- as.numeric(self$beta / (1 - private$SAR$lambda))
      direct <- as.numeric(sum(diag(solve(diag(private$cache$N) - private$SAR$lambda * private$SAR$W))) /
        private$cache$N * self$beta)

      list("total" = total, "direct" = direct, "indirect" = total - direct)
    },
    get_SAR = function() {private$SAR}
  ),
  private = list(
    SAR = NULL
  )

  )
}
