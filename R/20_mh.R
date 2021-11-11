
#' Metropolis-Hastings step
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @noRd
MetropolisHastings <- R6Class("MetropolisHastings",

  public = list(

    print = function() {

      cat("Metropolis-Hastings object for ", private$name ,".\n", sep = "")
      cat(private$accepted, " / ", private$proposed, " accepted / proposed values (rate ",
        round(self$get_acceptance, 2), ").\n", private$accepted_tune, " / ", private$proposed_tune,
        " total accepted / proposed values.\n", sep = "")
      cat("Current value / scale: ", private$value, " / ", private$scale, ".\n", sep = "")

      return(invisible(self))
    },

    initialize = function(value, scale = 0.1, ...) {

      private$value <- value
      private$scale <- scale

      # Keep track of totals and since scale adjustment
      private$proposed <- private$accepted <- private$proposed_tune <- private$accepted_tune <- 0L

      # Initialize children
      for(l in ls(self)) {if(grepl("^initialize_[0-9a-zA-Z]+$", l)) {self[[l]](...)}}

      return(NULL)
    },

    setup = function(...) {

      # Update children
      for(l in ls(self)) {if(grepl("^setup_[0-9a-zA-Z]+$", l)) {self[[l]](...)}}

      return(NULL)
    },

    sample = function(...) {
      self$propose(...)
      self$acceptance(...)
      self$finalize(...)
    },

    propose = function(location = private$value, scale = private$scale, ...) {
      private$proposal <- rnorm(1L, location, scale)
    },

    acceptance = function(proposal = private$proposal, current = private$value, ...) {
      private$probability <- exp(self$posterior(proposal, ...) -
        self$posterior(current, ...) + self$adjustment(proposal, ...))
    },

    adjustment = function(proposal = private$proposal, current = private$value, ...) {
      0
    },

    posterior = function(value = private$value, ...) {
      0 + prior(value, ...)
    },

    prior = function(value = private$value, ...) {
      0
    },

    finalize = function(...) {
      if(isTRUE(runif(1L) < private$probability)) {
        private$value <- private$proposal
        private$accepted <- private$accepted + 1L
        private$accepted_tune <- private$accepted_tune + 1L
      }
      private$proposed <- private$proposed + 1L
      private$proposed_tune <- private$proposed_tune + 1L
    }

  ),

  active = list(
    get_value = function() {private$value},
    get_scale = function() {private$scale},
    set_scale = function(value) {
      if(missing(value)) return(private$scale)
      private$scale <- value
      private$accepted <- 0L
      private$proposed <- 0L
    },
    get_accepted = function() {private$accepted},
    get_proposed = function() {private$proposed},
    get_acceptance = function(value) {private$accepted / max(private$proposed, 1)},
    get_tuning = function(value) {private$accepted_tune / max(private$proposed_tune, 1)}
  ),

  private = list(
    value = NULL, proposal = NULL, scale = NULL,
    probability = NULL,
    proposed = NULL, accepted = NULL, proposed_tune = NULL, accepted_tune = NULL,
    name = NULL
  )
)


#' Metropolis-Hastings step for theta in the Normal-Gamma shrinkage setup
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @noRd
MH_SNG_theta <- R6Class("MH_SNG_theta", inherit = MetropolisHastings,

  public = list(

    initialize_theta = function(rate = 1, ...) {
      private$rate <- rate
      private$name <- "Normal-Gamma theta"
    },

    propose = function(location = private$value, scale = private$scale, ...) {
      private$proposal <- exp(rnorm(1L, 0, scale)) * location
    },

    adjustment = function(proposal = private$proposal, current = private$value, ...) {
      log(proposal) - log(current)
    },

    posterior = function(value = private$value, ...) {
      dots <- list(...)
      sum(dgamma(dots$tau, value, (value * dots$lambda / 2), log = TRUE)) + self$prior(value)
    },

    prior = function(value = private$value, ...) {
      dexp(value, rate = private$rate, log = TRUE)
    }
  ),

  private = list(
    rate = NULL
  )
)


#' Metropolis-Hastings step for lambda in the spatial autoregressive model
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @noRd
MH_SAR_lambda <- R6Class("MH_SAR_lambda", inherit = MetropolisHastings,

  public = list(

    initialize_lambda = function(shape_a = 1.01, shape_b = 1.01, lower = -Inf, upper = Inf, ...) { # To-do: priors
      private$shape_a <- shape_a
      private$shape_b <- shape_b
      private$lower <- lower
      private$upper <- upper
      private$name <- "Spatial lambda"
    },

    setup_lambda = function(N, M) {
      private$N <- N
      private$M <- M
    },

    propose = function(location = private$value, scale = private$scale, ...) { # To-do: interweaving
      while(TRUE) {
        private$proposal <- rnorm(1L, location, scale)
        if(private$proposal < private$upper && private$proposal > private$lower) {break}
      }
      # private$proposal2 <- rnorm(1L, private$value2, private$scale)
      # private$proposal <- self$untransform(private$proposal2)
    },

    acceptance = function(proposal = private$proposal, current = private$value, ...) {
      private$probability <- exp(self$posterior(proposal, ...) -
        self$posterior(current, ...) + self$adjustment(proposal, current))
    },

    adjustment = function(proposal = private$proposal, current = private$value) {
      0
      # (1 - proposal^2) - (1 - current^2) # Proposed z transformed variable
    },

    posterior = function(value = private$value, ...) {
      dots <- list(...)
      dots$get_ldet(value) - (private$N - private$M) / 2 * log(dots$get_rss(value)) + self$prior(value)
    },

    prior = function(value = private$value, ...) {
      dbeta((value + 1) / 2, shape1 = private$shape_a, shape2 = private$shape_b, log = TRUE) - log(2)
    }
  ),

  private = list(
    shape_a = NULL, shape_b = NULL,
    lower = NULL, upper = NULL,
    N = NULL, M = NULL
  )
)



#' Metropolis-Hastings step for lambda in the spatial error model
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @noRd
MH_SEM_lambda <- R6Class("MH_SEM_lambda", inherit = MH_SAR_lambda)



#' Metropolis-Hastings step for delta in the spatially lagged explanatories model
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @noRd
MH_SLX_delta <- R6Class("MH_SLX_delta", inherit = MetropolisHastings,

  public = list(

    initialize_delta = function(shape_a = 2, shape_b = 1, lower = 1e-12, upper = Inf, ...) {
      private$shape_a <- shape_a
      private$shape_b <- shape_b
      private$lower <- lower
      private$upper <- upper
      private$name <- "Connectivity delta"
    },

    setup_delta = function(N, M) {
      private$N <- N
      private$M <- M
    },

    propose = function(location = private$value, scale = private$scale, ...) {
      while(TRUE) {
        private$proposal <- rnorm(1L, location, scale)
        if(private$proposal < private$upper && private$proposal > private$lower) {break}
      }
    },

    acceptance = function(proposal = private$proposal,
      current = private$value, ...) {
      private$probability <- exp(self$posterior(proposal, ...) -
        self$posterior(current, ...) + self$adjustment(proposal, current))
    },

    adjustment = function(proposal = private$proposal, current = private$value) {
      0
    },

    posterior = function(value = private$value, ...) {
      dots <- list(...)
      -(private$N - private$M) / 2 * log(dots$get_rss(value)) + self$prior(value)
    },

    prior = function(value = private$value, ...) {
      dgamma(1 / value, shape = private$shape_a, rate = private$shape_b, log = TRUE)
    }
  ),

  private = list(
    shape_a = NULL, shape_b = NULL,
    lower = NULL, upper = NULL,
    N = NULL, M = NULL
  )
)
