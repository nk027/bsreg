
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

      # Function to update the latent, i.e. z(lambda, delta) = y - lambda W(delta) y, with new parameters
      private$SAR$set_latent <- function(lambda = private$SAR$lambda, Wy = self$Wy) {
        private$SAR$z <- super$y - lambda * Wy
      }

      # Function to set a new lambda and update the latent
      private$SAR$set_lambda <- function(lambda) {
        private$SAR$lambda <- lambda
        # Update the latent
        private$SAR$set_latent(lambda = lambda)
      }

      # Function to set delta and update W(delta), cache, and latent
      private$SAR$set_delta <- function(delta) {
        if(private$SAR$Psi_fixed) {stop("Connectivity function 'Psi' (SAR) required to update delta.")}
        private$SAR$delta <- delta
        private$SAR$W <- private$SAR$Psi(delta)
        # Update the cache
        private$SAR$Wy <- private$SAR$W %*% super$y
        # private$SAR$XWy <- crossprod(super$X, private$SAR$Wy)
        # Update the latent
        private$SAR$set_latent(Wy = private$SAR$Wy)
      }

      # Initialise MH ---

      self$MH$SAR_lambda <- MH_SAR_lambda$new(value = priors$SAR$lambda, scale = priors$SAR$lambda_scale,
        prior = priors$SAR$lambda_prior, shape_a = priors$SAR$lambda_a, shape_b = priors$SAR$lambda_b,
        lower = priors$SAR$lambda_min, upper = priors$SAR$lambda_max)

      # Only sample delta if a proposal scale is set
      if(priors$SAR$delta_scale > 0) {
        self$MH$SAR_delta <- MH_SAR_delta$new(value = priors$SAR$delta, scale = priors$SAR$delta_scale,
          prior = priors$SAR$delta_prior, shape_a = priors$SAR$delta_a, shape_b = priors$SAR$delta_b,
          lower = priors$SAR$delta_min, upper = priors$SAR$delta_max)
      }

    },

    # High priority 9 for latent
    setup_9SAR = function(Psi_SAR = NULL,
      ldet_SAR = list(grid = FALSE, reps = 1L, i_lambda = c(-1 + 1e-12, 1 - 1e-12, 100L), i_delta = c(1e-12, 10, 20)),
      ...) {

      # Work out connectivity ---

      if(is.null(Psi_SAR)) {stop("Please provide a connectivity function or matrix 'Psi_SAR'.")}

      if(is.matrix(Psi_SAR)) {
        private$SAR$Psi_fixed <- TRUE
        private$SAR$W <- private$SAR$Psi <- Psi_SAR
        # # Set cache
        private$SAR$Wy <- private$SAR$W %*% super$y
        # Set lambda and obtain the latent
        private$SAR$set_lambda(private$SAR$priors$lambda)
      } else {
        private$SAR$Psi_fixed <- FALSE
        private$SAR$Psi <- Psi_SAR
        # Manually set lambda and then set delta to obtain W(delta), the cache, and the latent
        private$SAR$lambda <- private$SAR$priors$lambda
        private$SAR$set_delta(private$SAR$priors$delta)
      }

      # Set up MH ---

      self$MH$SAR_lambda$setup(N = private$cache$N, M = private$cache$M)

      if(!is.null(self$MH$SAR_delta)) {
        if(private$SAR$Psi_fixed) {stop("Connectivity function 'Psi' (SAR) required to sample 'delta'.")}
        self$MH$SAR_delta$setup(N = private$cache$N, M = private$cache$M)
      }

      # Set up the log determinant ---

      # Initialise object with options
      private$SAR$ldet <- list(size = private$cache$N / ldet_SAR$reps, reps = ldet_SAR$reps, # Kronecker settings
        grid = isTRUE(ldet_SAR$grid), i_lambda = ldet_SAR$i_lambda, i_delta = ldet_SAR$i_delta) # Grid settings

      # If W is repeated via Kronecker product we only need one submatrix and scale up using `reps` later
      private$SAR$ldet$get_W <- if(private$SAR$ldet$reps == 1) {
        function(W = private$SAR$W) {W}
      } else { # Just retrieve the sub-matrix
        function(W = private$SAR$W) {W[seq(private$SAR$ldet$size), seq(private$SAR$ldet$size)]}
      }

      # Provide a function for the log-determinant

      if(isTRUE(private$SAR$Psi_fixed)) { # If Psi is fixed we only have the lambda dimension

        if(private$SAR$ldet$grid) { # If a grid is requested we fit a spline to a grid over lambda

          pars <- i_seq(private$SAR$ldet$i_lambda)
          ldets <- vapply(pars, function(x) {
            determinant(diag(private$SAR$ldet$size) - x * private$SAR$ldet$get_W(),
            logarithm = TRUE)$modulus * private$SAR$ldet$reps
          }, numeric(1L))
          private$SAR$ldet$splinefun <- splinefun(pars, y = ldets)
          private$SAR$ldet$get_ldet <- function(lambda = private$SAR$lambda, ...) {
            private$SAR$ldet$splinefun(x = lambda)
          }

        } else { # Otherwise we use a spectral decomposition

          private$SAR$ldet$ev <- eigen(private$SAR$ldet$get_W(),
            symmetric = is_symmetric(private$SAR$Psi), only.values = TRUE)$values
          # The log-determinant of I - lambda W is just the sum of log(1 - lambda * omega)
          private$SAR$ldet$get_ldet <- function(lambda = private$SAR$lambda, ...) {
            Re(sum(log(1 - lambda * private$SAR$ldet$ev))) * private$SAR$ldet$reps
          }
        }

      } else { # If Psi is free we need to account for both dimensions

        if(private$SAR$ldet$grid) { # If a grid is requested we approximate using a Gaussian process

          # We train on a grid over lambda and delta
          ldets <- cbind("ldet" = NA_real_, as.matrix(expand.grid("lambda" = i_seq(private$SAR$ldet$i_lambda),
            "delta" = i_seq(private$SAR$ldet$i_delta))))
          for(i in seq(private$SAR$ldet$i_delta[3L])) {
            ev <- eigen(private$SAR$ldet$get_W(
              private$SAR$Psi(ldets[(i - 1) * private$SAR$ldet$i_lambda[3L] + 1, "delta"])),
              symmetric = is_symmetric(private$SAR$Psi), only.values = TRUE)$values
            ldets[, "ldet"] <- vapply(i_seq(private$SAR$ldet$i_lambda), function(lambda) {
              Re(sum(log(1 - lambda * ev))) * private$SAR$ldet$reps
            }, numeric(1L))
          }
          # ldets <- apply(pars, 1, function(x) {
          #   determinant(diag(private$SAR$ldet$size) - x[1] * private$SAR$ldet$get_W(private$SAR$Psi(x[2])),
          #     logarithm = TRUE)$modulus * private$SAR$ldet$reps
          # })
          private$SAR$ldet$gp <- GauPro::GauPro(ldets[, c("lambda", "delta")], ldets[, "ldet"], parallel = FALSE)
          # The log-determinant is predicted using the Gaussian process
          private$SAR$ldet$get_ldet <- function(lambda = private$SAR$lambda, delta = private$SAR$delta, ...) {
            private$SAR$ldet$gp$predict(cbind(lambda, delta))[[1]]
          }

        } else { # Without a grid we brute-force every time via a LU decomposition

          private$SAR$ldet$get_ldet <- function(lambda = private$SAR$lambda, delta = private$SAR$delta, ...) {
            determinant(diag(private$SAR$ldet$size) - lambda * private$SAR$ldet$get_W(private$SAR$Psi(delta)),
              logarithm = TRUE)$modulus * private$SAR$ldet$reps
          }
        }

      }
    },

    # Sample spatial autoregressive parameter lambda ---
    sample_latent = function() {

      # Prepare RSS and the log-determinant as functions of lambda
      get_rss <- {function() {
        prec_ch <- chol(self$prior_precision + self$XX / self$sigma)
        b0 <- backsolve(prec_ch, forwardsolve(prec_ch, (self$prior_precision %*% private$NG$mu0 + crossprod(self$X, super$y) / self$sigma),
          upper.tri = TRUE, transpose = TRUE))
        b1 <- backsolve(prec_ch, forwardsolve(prec_ch, (self$prior_precision %*% private$NG$mu0 + self$XWy / self$sigma),
          upper.tri = TRUE, transpose = TRUE))
        e0 <- super$y - self$X %*% b0
        e1 <- self$Wy - self$X %*% b1
        e0e0 <- sum(e0^2)
        e1e0 <- sum(e1 * e0)
        e1e1 <- sum(e1^2)
        return(function(value) {(e0e0) - (2 * value * e1e0) + (value^2 * e1e1)})
      }}()
      get_ldet <- function(value) {private$SAR$ldet$get_ldet(lambda = value, delta = private$SAR$delta)}

      # Metropolis-Hastings step for lambda
      self$MH$SAR_lambda$propose()
      self$MH$SAR_lambda$acceptance(get_rss = get_rss, get_ldet = get_ldet)
      self$MH$SAR_lambda$finalize()
      lambda <- self$MH$SAR_lambda$get_value # Assign and recompute
      if(abs(lambda - private$SAR$lambda) > 1e-12) {private$SAR$set_lambda(lambda)}
    },

    # Sample connectivity parameter delta ---
    sample_extra1 = function() {

      if(!is.null(self$MH$SAR_delta)) {
        # Prepare RSS and log-determinant as functions of delta
        res <- super$y - self$X %*% self$beta
        get_rss <- function(value) {sq_sum(res - private$SAR$lambda * private$SAR$Psi(value) %*% super$y)}
        get_ldet <- function(value) {private$SAR$ldet$get_ldet(lambda = private$SAR$lambda, delta = value)}

        # Metropolis-Hastings step for delta
        self$MH$SAR_delta$propose()
        self$MH$SAR_delta$acceptance(get_rss = get_rss, get_ldet = get_ldet)
        self$MH$SAR_delta$finalize()
        delta <- self$MH$SAR_delta$get_value # Assign and recompute
        if(abs(delta - private$SAR$delta) > 1e-12) {private$SAR$set_delta(delta)}
      }
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
      if(!private$SAR$Psi_fixed) {pars$delta_SAR <- private$SAR$delta}
      if(!private$SAR$priors$lambda_prior == "bgamma") {
        pars$tau_SAR <- private$MH$SAR_lambda$get_tau
      }
      return(pars)
    },
    get_effects = function() { # To-do: use eigendecomposition or provide alternative methods to be more efficient
      total <- as.numeric(self$beta / (1 - private$SAR$lambda))
      direct <- as.numeric(sum(diag(solve(diag(private$cache$N) - private$SAR$lambda * private$SAR$W))) /
        private$cache$N * self$beta)
      return(list("total" = total, "direct" = direct, "indirect" = total - direct))
    },
    get_SAR = function() {private$SAR}
  ),
  private = list(
    SAR = NULL
  )

  )
}
