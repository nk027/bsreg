
#' Bayesian spatial error model
#'
#' @docType class
#'
#' @param parent \code{\link{R6Class}} object to inherit from.
#'
#' @importFrom R6 R6Class
#' @importFrom stats splinefun
#'
#' @noRd
get_sem_class <- function(parent = NormalGamma) {

  SpatialEM <- R6Class("SpatialEM", inherit = parent,

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
        private$SEM$y <- private$SEM$normalizer %*% super$y
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

    # Medium priority 8 for volatility filter
    setup_8SEM = function(Psi_SEM = NULL, ldet_SEM = list(grid = FALSE, i_lambda = c(-1, 1 - 1e-12, 100L), reps = 1L),
      ...) {

      # Work out connectivity ---

      if(is.null(Psi_SEM)) {stop("Please provide a connecitivity matrix 'Psi_SEM'.")}

      private$SEM$Psi_fixed <- TRUE
      private$SEM$W <- private$SEM$Psi <- Psi_SEM
      # Set cache
      private$SEM$Wy <- private$SEM$W %*% super$y
      private$SEM$WX <- private$SEM$W %*% super$X
      # Set lambda and obtain the normalizer
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
        private$SEM$ldet$get_ldet <- function(lambda, ...) {
          Re(sum(log(1 - lambda * private$SEM$ldet$ev))) * private$SEM$ldet$reps
        }
      }
    },

    # Sample spatial error parameter lambda ---
    sample_volatility = function() {

      # Prepare RSS and the log-determinant as functions of lambda
      get_rss <- function(value) {
        y_s <- super$y - value * private$SEM$Wy
        X_s <- super$X - value * private$SEM$WX
        sq_sum(y_s) - sq_sum(crossprod(qr.Q(qr(X_s)), y_s)) # See e.g. Bivand and Piras, 2015
      }
      get_ldet <- function(value) {private$SEM$ldet$get_ldet(lambda = value)}

      # Metropolis-Hastings step for lambda
      self$MH$SEM_lambda$propose()
      self$MH$SEM_lambda$acceptance(get_rss = get_rss, get_ldet = get_ldet)
      self$MH$SEM_lambda$finalize()
      lambda <- self$MH$SEM_lambda$get_value # Assign and recompute
      private$SEM$set_lambda(lambda) # Do it always to accommodate changing X
    }

  ),

  active = list(

    # Variables that are adapted using the normalizer ---
    y = function() {private$SEM$y},
    X = function() {private$SEM$X},
    XX = function() {private$SEM$XX},
    Xy = function() {private$SEM$Xy},

    # Access functions ---
    get_parameters = function() {
      pars <- super$get_parameters
      pars$lambda_SEM <- private$SEM$lambda
      return(pars)
    },
    get_SEM = function() {private$SEM}
  ),

  private = list(
    SEM = NULL
  )

  )
}
