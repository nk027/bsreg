
#' Bayesian spatially lagged explanatories model
#'
#' @docType class
#'
#' @param parent \code{\link{R6Class}} object to inherit from.
#'
#' @importFrom R6 R6Class
#'
#' @noRd
get_slx_class <- function(parent = NormalGamma) {

  SpatialLX <- R6Class("SpatialLX", inherit = parent,

  public = list(

    initialize_SLX = function(priors, ...) {

      # Store prior settings ---

      if(missing(priors) || is.null(priors$SLX)) {priors <- list(SLX = set_SLX())}
      private$SLX$priors <- priors$SLX

      # Build SLX object ---

      # Function to set a new delta and update W(delta) and cache
      private$SLX$set_delta <- function(delta) {
        if(private$SLX$Psi_fixed) {stop("Connectivity function 'Psi' (SLX) required to update delta.")}
        private$SLX$delta <- delta
        private$SLX$W <- private$SLX$Psi(delta)
        # Update the cache
        private$SLX$WX <- private$SLX$W %*% private$SLX$X_SLX
        private$SLX$X <- cbind(super$X, private$SLX$WX)
        private$SLX$XX <- crossprod(self$X)
        private$SLX$Xy <- crossprod(self$X, self$y)
      }

      # Initialise MH ---

      # Only sample delta if a proposal scale is set
      if(priors$SLX$delta_scale > 0) {
        self$MH$SLX_delta <- MH_SLX_delta$new(value = priors$SLX$delta, scale = priors$SLX$delta_scale,
          shape_a = priors$SLX$delta_a, shape_b = priors$SLX$delta_b,
          lower = priors$SLX$delta_min, upper = priors$SLX$delta_max)
      }

      return(NULL)
    },

    # Priority 2 since it updates 'M' (the number of columns in the design matrix)
    setup_2SLX = function(X_SLX, Psi_SLX = NULL, ...) {

      # Add new data and adapt the cache ---
      private$SLX$X_SLX <- X_SLX
      private$cache$M <- private$cache$M + NCOL(X_SLX) # Add columns of spatial lag

      # Work out connectivity ---

      if(is.null(Psi_SLX)) {stop("Please provide a connecitivity function or matrix 'Psi_SLX'.")}

      if(is.matrix(Psi_SLX)) {
        private$SLX$Psi_fixed <- TRUE
        private$SLX$W <- Psi_SLX
        # Set cache manually
        private$SLX$WX <- private$SLX$W %*% private$SLX$X_SLX
        private$SLX$X <- cbind(super$X, private$SLX$WX)
        private$SLX$XX <- crossprod(private$SLX$X)
        private$SLX$Xy <- crossprod(private$SLX$X, self$y)
      } else {
        private$SLX$Psi_fixed <- FALSE
        private$SLX$Psi <- Psi_SLX
        # Set delta to obtain spatially lagged X
        private$SLX$set_delta(private$SLX$priors$delta)
      }

      # Set up MH ---

      if(!is.null(self$MH$SLX_delta)) {
        if(private$SLX$Psi_fixed) {stop("Connectivity function 'Psi' (SLX) required to sample 'delta'.")}
        self$MH$SLX_delta$setup(N = private$cache$N, M = private$cache$M)
      }

      return(NULL)
    },

    # Sample connectivity parameter delta ---
    sample_extra3 = function() {

      if(!is.null(self$MH$SLX_delta)) {
        # Prepare RSS as a function of delta
        get_rss <- function(value) {
          X <- cbind(super$X, private$SLX$Psi(value) %*% private$SLX$X_SLX)
          beta <- solve(private$NG$prec0 + crossprod(X) / self$sigma,
            (crossprod(X, self$y) / self$sigma + self$prior_precision %*% private$NG$mu0))
          sq_sum(self$y - X %*% beta)
          # sq_sum(self$y - cbind(super$X, private$SLX$Psi(value) %*% private$SLX$X_SLX) %*% self$beta)
        }

        # Metropolis-Hastings step for delta
        self$MH$SLX_delta$propose()
        self$MH$SLX_delta$acceptance(get_rss = get_rss)
        self$MH$SLX_delta$finalize()
        delta <- self$MH$SLX_delta$get_value # Assign and recompute
        if(abs(delta - private$SLX$delta) > 1e-12) {private$SLX$set_delta(delta)}
      }
    }
  ),

  active = list(

    # Variables that are updated using the connectivity ---
    X = function() {private$SLX$X},
    XX = function() {private$SLX$XX},
    Xy = function() {crossprod(private$SLX$X, self$y)},

    # Acessor functions ---
    get_parameters = function() {
      pars <- super$get_parameters
      pars$delta_SLX <- private$SLX$delta
      return(pars)
    },
    get_SLX = function() {private$SLX}
  ),

  private = list(
    SLX = NULL
  )

  )
}