
#' Bayesian stochastic volatility model
#'
#' @docType class
#'
#' @param parent \code{\link{R6Class}} object to inherit from.
#'
#' @importFrom R6 R6Class
#'
#' @noRd
get_sv_class <- function(parent = NormalGamma) {

  has_package("stochvol")
  StochasticVolatility <- R6Class("StochasticVolatility", inherit = parent,

  public = list(

    initialize_SV = function(priors, ...) {

      # Store prior settings ---

      if(missing(priors) || is.null(priors$SV)) {priors <- list(SV = set_SV())}
      private$SV$priors <- priors$SV$priors
      private$SV$parameters <- priors$SV$parameters

      # Build SV object ---

      # Function to update the normalizer and cache with a new latent
      private$SV$set_latent <- function(latent) {
        private$SV$latent <- latent
        normalizer <- exp(-latent / 2)
        private$SV$y <- super$y * normalizer
        private$SV$X <- super$X * normalizer
        # Cache
        private$SV$XX <- crossprod(private$SV$X)
        private$SV$Xy <- crossprod(private$SV$X, private$SV$y)
      }

      return(NULL)
    },

    setup_SV = function(...) {

      # Set up the latent ---

      private$SV$set_latent(rep(private$SV$parameters$latent0, private$cache$N))

      return(NULL)
    },

    # Sample latent quantity
    sample_volatility = function() {

      # Use stochvol's one step sampler
      sv <- stochvol::svsample_fast_cpp(self$residuals,
        startpara = private$SV$parameters, startlatent = private$SV$latent, priorspec = private$SV$priors)

      private$SV$parameters <- sv$para
      private$SV$set_latent(drop(sv$latent))
    }
  ),

  active = list(

    # Variables that are adapted using the normalizer ---
    y = function() {private$SV$y},
    X = function() {private$SV$X},
    XX = function() {private$SV$XX},
    Xy = function() {private$SV$Xy},

    # Access functions ---
    get_SV = function() {private$SV}
  ),

  private = list(
    SV = NULL
  )

  )

}