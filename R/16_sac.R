
#' Bayesian spatial autoregressive combined model
#'
#' @docType class
#'
#' @param parent \code{\link{R6Class}} object to inherit from.
#'
#' @importFrom R6 R6Class
#' @importFrom stats splinefun
#'
#' @noRd
get_sac_class <- function(parent = NormalGamma) {

  SpatialAC <- R6Class("SpatialAC", inherit = parent,

  public = list(

    initialize_SEM = function(priors, ...) {

      # Store prior settings ---

      if(missing(priors) || is.null(priors$SEM)) {priors <- list(SEM = set_SEM())}
      private$SEM$priors <- priors$SEM

      # Build SEM object ---

      # Function to update the normalizer, i.e. s(lambda) = I - lambda W, with new parameters
      private$SEM$set_normalizer <- function(lambda = private$SEM$lambda, W = private$SEM$W) {
        private$SEM$normalizer <- diag(private$cache$N) - lambda * W
        # Apply to y and X
        private$SEM$y <- private$SEM$normalizer %*% self$z
        private$SEM$X <- private$SEM$normalizer %*% super$X
        # Update the cache
        private$SEM$Xy <- crossprod(private$SEM$X, private$SEM$y)
        private$SEM$XX <- crossprod(private$SEM$X)
      }

      # Function to set a new lambda and update the normalizer
      private$SEM$set_lambda <- function(lambda) {
        private$SEM$lambda <- lambda
        # Update the normalizer
        private$SEM$set_normalizer(lambda = lambda)
      }

      # Initialise MH ---

      self$MH$SEM_lambda <- MH_SEM_lambda$new(value = priors$SEM$lambda, scale = priors$SEM$lambda_scale,
        shape_a = priors$SEM$lambda_a, shape_b = priors$SEM$lambda_b,
        lower = priors$SEM$lambda_min, upper = priors$SEM$lambda_max)
    },

    initialize_SAR = function(priors, ...) {

      # Store prior settings ---

      if(missing(priors) || is.null(priors$SAR)) {priors <- list(SAR = set_SAR())}
      private$SAR$priors <- priors$SAR

      # Build SAR object ---

      # Function to update the latent, i.e. z(lambda) = y - lambda W y, with new parameters
      private$SAR$set_latent <- function(lambda = private$SAR$lambda, Wy = self$Wy_s) {
        private$SAR$z <- self$y_s - lambda * Wy
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


    # Medium priority 8 for volatility filter
    setup_8SEM = function(Psi_SEM = NULL, ldet_SEM = list(grid = FALSE, i_lambda = c(-1, 1 - 1e-12, 100L), reps = 1L),
      ...) {

      # Work out connectivity ---

      if(is.null(Psi_SEM)) {stop("Please provide a connecitivity matrix 'Psi_SEM'.")}

      private$SEM$Psi_fixed <- TRUE
      private$SEM$W <- private$SEM$Psi <- Psi_SEM
      # Set cache
      # private$SEM$Wy <- private$SEM$W %*% super$y
      private$SEM$WX <- private$SEM$W %*% super$X
      # Set lambda and obtain the normalizer
      private$SAR$z <- super$y
      private$SEM$set_lambda(private$SEM$priors$lambda)

      # Set up MH ---

      self$MH$SEM_lambda$setup(N = private$cache$N, M = private$cache$M)

      if(!is.null(self$MH$SEM_delta)) {
        if(private$SEM$Psi_fixed) {stop("Connectivity function 'Psi' (SEM) required to sample 'delta'.")}
        self$MH$SEM_delta$setup(N = private$cache$N, M = private$cache$M)
      }

      # Set up the log determinant ---

      # Initialise object with options
      private$SEM$ldet <- list(size = private$cache$N / ldet_SEM$reps, reps = ldet_SEM$reps, # Kronecker settings
        grid = isTRUE(ldet_SEM$grid), i_lambda = ldet_SEM$i_lambda) # Grid settings

      # If W is repeated via Kronecker product we only need one submatrix and scale up using `reps` later
      private$SEM$ldet$get_W <- if(private$SEM$ldet$reps == 1) {
        function(W = private$SEM$W) {W}
      } else { # Just retrieve the sub-matrix
        function(W = private$SEM$W) {W[seq(private$SEM$ldet$size), seq(private$SEM$ldet$size)]}
      }

      # Provide a function for the log-determinant
      if(private$SEM$ldet$grid) { # If a grid is requested we fit a spline to a grid over lambda

        pars <- i_seq(private$SEM$ldet$i_lambda)
        ldets <- vapply(pars, function(x) {
          determinant(diag(private$SEM$ldet$size) - x * private$SEM$ldet$get_W(),
            logarithm = TRUE)$modulus * private$SEM$ldet$reps
        }, numeric(1L))
        private$SEM$ldet$splinefun <- splinefun(pars, y = ldets)
        private$SEM$ldet$get_ldet <- function(lambda, ...) {private$SEM$ldet$splinefun(x = lambda)}

      } else {

        private$SEM$ldet$ev <- eigen(private$SEM$ldet$get_W(),
          symmetric = is_symmetric(private$SEM$Psi), only.values = TRUE)$values
        # The log-determinant of I - lambda W is just the sum of log(1 - lambda * omega)
        private$SEM$ldet$get_ldet <- \(lambda, ...) {
          Re(sum(log(1 - lambda * private$SEM$ldet$ev))) * private$SEM$ldet$reps
        }
      }
    },

    # High priority 9 for latent
    setup_9SAR = function(Psi_SAR = NULL, ldet_SAR = list(grid = FALSE, i_lambda = c(-1, 1 - 1e-12, 100L), reps = 1L),
      ...) {

      # Work out connectivity ---

      if(is.null(Psi_SAR)) {stop("Please provide a connectivity matrix 'Psi_SAR'.")}

      private$SAR$Psi_fixed <- TRUE
      private$SAR$W <- private$SAR$Psi <- Psi_SAR
      # # Set cache
      # private$SAR$Wy <- private$SAR$W %*% super$y
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
        b0 <- backsolve(prec_ch, forwardsolve(prec_ch, (private$NG$prec0 %*% private$NG$mu0 + self$Xy_s / self$sigma),
          upper.tri = TRUE, transpose = TRUE))
        b1 <- backsolve(prec_ch, forwardsolve(prec_ch, (private$NG$prec0 %*% private$NG$mu0 + self$XWy_s / self$sigma),
          upper.tri = TRUE, transpose = TRUE))
        e0 <- self$y_s - self$X %*% b0
        e1 <- self$Wy_s - self$X %*% b1
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
    },

    # Sample spatial error parameter lambda ---
    sample_volatility = function() {

      # Prepare RSS and the log-determinant as functions of lambda
      get_rss <- function(value) {
        y_s <- self$z - value * private$SEM$W %*% self$z
        X_s <- super$X - value * private$SEM$W %*% super$X
        sq_sum(y_s) - sq_sum(crossprod(qr.Q(qr(X_s)), y_s)) # See e.g. Bivand and Piras, 2015
      }
      get_ldet <- function(value) {private$SEM$ldet$get_ldet(lambda = value)}

      # Metropolis-Hastings step for lambda
      self$MH$SEM_lambda$propose()
      self$MH$SEM_lambda$acceptance(get_rss = get_rss, get_ldet = get_ldet)
      self$MH$SEM_lambda$finalize()
      lambda <- self$MH$SEM_lambda$get_value # Assign and recompute
      if(abs(lambda - private$SEM$lambda) > 1e-12) {private$SEM$set_lambda(lambda)}
    }

  ),

  active = list(

    # Variables that are adapted using the normalizer ---
    y = function() {private$SAR$z - private$SEM$lambda * private$SEM$W %*% private$SAR$z},
    Xy = function() {crossprod(private$SEM$X, self$y)},
    X = function() {private$SEM$X},
    XX = function() {private$SEM$XX},
    z = function() {private$SAR$z},
    y_s = function() {private$SEM$y},
    Xy_s = function() {private$SEM$Xy},
    Wy_s = function() {private$SAR$W %*% private$SEM$y},
    XWy_s = function() {crossprod(private$SEM$X, private$SAR$W %*% private$SEM$y)}
  ),

  private = list(SAR = NULL, SEM = NULL)

  )
}
