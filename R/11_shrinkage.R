
#' Bayesian model with Normal-Gamma shrinkage prior (Polson and Scott, 2010)
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @noRd
ShrinkageNormalGamma <- R6Class("ShrinkageNormalGamma", inherit = NormalGamma,

  public = list(

    initialize_SNG = function(priors, ...) {

      # Store prior settings ---

      if(missing(priors) || is.null(priors$SNG)) {priors <- list(SNG = set_SNG())}
      private$SNG$priors <- priors$SNG

      # Update meta info ---

      private$meta$priortype <- "Normal-Gamma shrinkage"

      # Initialise MH ---

      # Only sample theta if a proposal scale is set
      if(priors$SNG$theta_scale > 0) {
        self$MH$SNG_theta <- MH_SNG_theta$new(value = priors$SNG$theta, scale = priors$SNG$theta_scale,
          rate = priors$SNG$theta_a)
      }

      return(NULL)
    },

    setup_SNG = function(...) {

      # Set up quantities ---

      private$SNG$tau <- rep(private$SNG$priors$tau, times = private$cache$M)
      private$SNG$lambda <- private$SNG$priors$lambda
      private$SNG$theta <- private$SNG$priors$theta
      # Overwrite prior precision
      self$prior_precision <- diag(1 / private$SNG$tau, nrow = private$cache$M)

      # Set up MH ---
      if(!is.null(self$MH$SNG_theta)) {self$MH$SNG_theta$setup()}

      return(NULL)
    },

    sample_shrinkage = function() {

      # Sample Normal-Gamma shrinkage ---

      # Lambda from Gamma
      private$SNG$lambda <- rgamma(1L, shape = private$SNG$priors$lambda_a + private$cache$M * private$SNG$theta,
        rate = private$SNG$priors$lambda_b + (private$SNG$theta * sum(private$SNG$tau)) / 2)
      # Tau from GIG
      for(i in seq_along(private$SNG$tau)) {
        private$SNG$tau[i] <- max(GIGrvg::rgig(1L, lambda = private$SNG$theta - 0.5,
          chi = (self$beta[i] - self$prior_mean[i])^2, psi = private$SNG$lambda * private$SNG$theta), 1e-12)
      }

      # Update prior precision
      self$prior_precision <- diag(1 / private$SNG$tau, nrow = private$cache$M)

      # Metropolis-Hastings step for theta
      if(!is.null(self$MH$SNG_theta)) {
        self$MH$SNG_theta$sample(tau = private$SNG$tau, lambda = private$SNG$lambda)
        private$SNG$theta <- self$MH$SNG_theta$get_value
      }
    }
  ),

  active = list(
    # Access functions ---
    get_SNG = function() {private$SNG}
  ),

  private = list(
    SNG = NULL
  )

)


#' Bayesian model with Horseshoe shrinkage prior (Makalic and Schmidt, 2015)
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @noRd
Horseshoe <- R6Class("Horseshoe", inherit = NormalGamma,

  public = list(

    initialize_HS = function(priors, ...) {

      # Store prior settings ---

      if(missing(priors) || is.null(priors$HS)) {priors <- list(HS = set_HS())}
      private$HS$priors <- priors$HS

      # Update meta info ---

      private$meta$priortype <- "Horseshoe"

      return(NULL)
    },

    setup_HS = function(...) {

      # Set up quantities ---

      private$HS$lambda <- rep(private$HS$priors$lambda, times = private$cache$M)
      private$HS$nu <- rep(private$HS$priors$nu, times = private$cache$M)
      private$HS$tau <- private$HS$priors$tau
      private$HS$zeta <- private$HS$priors$zeta
      # Overwrite prior precision
      self$prior_precision <- diag(1 / c(private$HS$lambda * private$HS$tau), nrow = private$cache$M)
    },

    sample_shrinkage = function() {

      # Sample Horseshoe shrinkage ---

      # Lambda from inverse Gamma
      private$HS$lambda <- 1 / rgamma(private$cache$M, shape = 1,
        rate = 1 / private$HS$nu + self$beta^2 / (2 * private$HS$tau * self$sigma))
      # Tau from inverse Gamma
      private$HS$tau <- 1 / rgamma(1L, shape = (private$cache$M + 1L) / 2,
        rate = 1 / private$HS$zeta + sum(self$beta^2 / private$HS$lambda) / (2 * self$sigma))
      # Tau from inverse Gamma
      private$HS$nu <- 1 / rgamma(private$cache$M, shape = 1, rate = 1 + 1 / private$HS$lambda)
      # Zeta from inverse Gamma
      private$HS$zeta <- 1 / rgamma(1L, shape = 1, rate = 1 + 1 / private$HS$tau)

      # Update prior precision
      self$prior_precision <- diag(1 / pmax((private$HS$lambda * private$HS$tau), 1e-12), nrow = private$cache$M)
    }
  ),

  active = list(
    # Access functions ---
    get_HS = function() {private$HS}
  ),

  private = list(
    HS = NULL
  )

)
