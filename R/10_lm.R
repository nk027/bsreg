
#' Base class with common functionality
#'
#' This class provides basic functionality to build a hierarchical Bayesian model.
#' The three public functions are (1) 'initialize' to provide settings at construction time, (2) 'setup' to provide
#' additional settings and data, (3) 'finalize' to update the status and meta information after an iteration. The first
#' two functions also call 'initialize' and 'setup' methods of descendants as well as the 'starting' methods to further
#' allow layering setup processes.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @noRd
Base <- R6Class("Base",

  public = list(

    initialize = function(...) {

      # Meta ---
      private$meta <- list("iterations" = 0L)

      # Initialize children
      lapply(ls(self, pattern = "^initialize_[0-9a-zA-Z]+$"), function(n) self[[n]](...))

      invisible(self)
    },

    setup = function(y, X, ...) {

      # Data ---
      private$data <- list("y" = y, "X" = X)

      # Cache ---
      private$cache <- list("N" = NROW(X), "M" = NCOL(X), "XX" = crossprod(X), "Xy" = crossprod(X, y))

      # Set up children ---
      lapply(ls(self, pattern = "^setup_[0-9a-zA-Z]+$"), function(n) self[[n]](...))
      lapply(ls(self, pattern = "^starting_[0-9a-zA-Z]+$"), function(n) self[[n]](...))

      # Done ---
      invisible(self)
    },

    sample = function(...) { # To-do: build a list of sampling functions once at initialization

      # Sample expected quantities ---
      self$sample_sigma()
      self$sample_beta()
      self$sample_shrinkage() # Normal-Gamma or Horseshoe
      self$sample_latent() # Spatial autoregressive or limited dependent
      self$sample_volatility() # Spatial error or stochastic volatility
      # Sample potential extra quantities
      self$sample_extra1()
      self$sample_extra2()
      self$sample_extra3()

      # Update status ---
      self$finalize()

      return(invisible(NULL))
    },
    sample_beta = function() {NULL},
    sample_sigma = function() {NULL},
    sample_shrinkage = function() {NULL},
    sample_latent = function() {NULL},
    sample_volatility = function() {NULL},
    sample_extra1 = function() {NULL},
    sample_extra2 = function() {NULL},
    sample_extra3 = function() {NULL},

    finalize = function() {

      private$meta$iterations <- private$meta$iterations + 1L

    },

    # Slot for Metropolis-Hastings steps
    MH = list()
  ),
  active = list(

    # Access functions ---
    y = function() {private$data$y},
    X = function() {private$data$X},
    XX = function() {private$cache$XX},
    Xy = function() {private$cache$Xy},

    get_meta = function() {c(iterations = private$meta$iterations)}
  ),

  private = list(
    data = NULL, cache = NULL, meta = NULL
  )

)



#' Bayesian model with independent Normal-Gamma prior
#'
#' This class serves as the base for most practical models. It extends the proper 'Base' class with basic functionality
#' to estimate a linear model. A print function is available, a slot for Metropolis-Hastings objects is provided, and
#' several functions to access parts of the model are made available.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @noRd
NormalGamma <- R6Class("NormalGamma", inherit = Base,

  public = list(

    print = function() {

      cat("Bayesian model with a", private$meta$priortype, "prior setup.\n")
      cat(private$meta$iterations, "total samples have been drawn so far.\n")
      cat("The sampler contains", length(self$MH), "Metropolis-Hastings steps.\n")
      for(i in seq_along(self$MH)) {
        cat("\t", names(self$MH)[i], " at ", self$MH[[i]]$get_accepted, " accepted draws out of ",
          self$MH[[i]]$get_proposed, " proposals (rate ", round(self$MH[[i]]$get_acceptance, 2), ").\n", sep = "")
      }

      return(invisible(self))
    },

    initialize_NG = function(priors, ...) {

      # Store prior settings ---

      if(missing(priors) || is.null(priors$NG)) {priors <- list(NG = set_NG())}
      private$NG$priors <- priors$NG

      # Update meta info ---

      private$meta$priortype <- "Normal-Gamma"
    },

    # Priority 4 since it relies on 'M'
    setup_4NG = function(...) {

      # Update prior settings to fit the data ---

      self$prior_mean <- matrix(private$NG$priors$mu, nrow = private$cache$M)
      self$prior_precision <- diag(private$NG$priors$precision, nrow = private$cache$M)
      private$NG$shape0 <- private$NG$priors$shape
      private$NG$rate0 <- private$NG$priors$rate

      # Calculate known posterior quantities
      private$NG$shape1 <- private$NG$shape0 + private$cache$N / 2
    },

    starting_NG = function(...) {

      # Set sensible values of beta and sigma using LS
      self$beta <- if(is.null(private$NG$priors$beta)) {
        qr.solve(self$XX, self$Xy)
      } else {matrix(private$NG$priors$beta, nrow = private$cache$M)}
      self$sigma <- if(is.null(private$NG$priors$sigma)) {
        sq_sum(self$residuals) / (private$cache$N - private$cache$M)
      } else {private$NG$priors$sigma}
    },

    sample_beta = function() {

      private$NG$prec1 <- self$XX / self$sigma + self$prior_precision
      private$NG$mu1 <- solve(private$NG$prec1, (self$Xy / self$sigma + self$prior_precision %*% private$NG$mu0))
      # Draw from multivariate Normal
      self$beta <- t(rmvn(1L, mu = private$NG$mu1, precision = private$NG$prec1))
    },

    sample_sigma = function() {

      private$NG$rate1 <- private$NG$rate0 + sq_sum(self$residuals) / 2
      # Draw from inverse Gamma
      self$sigma <- 1 / rgamma(1L, shape = private$NG$shape1, rate = private$NG$rate1)
    }
  ),

  active = list(

    # We reserve beta and sigma
    beta = function(value) {if(missing(value)) {private$NG$beta} else {private$NG$beta <- value}},
    sigma = function(value) {if(missing(value)) {private$NG$sigma} else {private$NG$sigma <- value}},
    # Shrinkage priors may change the prior precision
    prior_precision = function(value) {if(missing(value)) {private$NG$prec0} else {private$NG$prec0 <- value}},
    prior_mean = function(value) {if(missing(value)) {private$NG$mu0} else {private$NG$mu0 <- value}},

    # Access functions ---
    residuals = function() {self$y - self$X %*% self$beta},
    get_parameters = function() {list("beta" = self$beta, "sigma" = sqrt(self$sigma))},
    get_NG = function() {private$NG}
  ),

  private = list(
    NG = NULL
  )

)


#' Bayesian model with a conjugate Normal-Gamma prior
#'
#' The conjugate Normal-Gamma prior is a simple adaptation of the independent Normal-Gamma prior. An additional
#' 'starting' method is provided to compute all known posterior quantities. Sampling steps for 'beta' and 'sigma'
#' are adapted to use the known posteriors.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @noRd
ConjugateNormalGamma <- R6Class("ConjugateNormalGamma", inherit = NormalGamma,

  public = list(
    starting_NG = function(...) { # No need for LS estimates to initialise
      self$sample_beta()
      self$sample_sigma()
    },
    sample_beta = function() {
      private$NG$prec1 <- self$XX + self$prior_precision
      private$NG$mu1 <- solve(private$NG$prec1, (self$Xy + self$prior_precision %*% private$NG$mu0))
      self$beta <- t(rmvn(1L, mu = private$NG$mu1, precision = private$NG$prec1))
    },
    sample_sigma = function() {
      residuals <- self$y - self$X %*% private$NG$mu1
      private$NG$rate1 <- private$NG$rate0 + sq_sum(residuals) / 2
      self$sigma <- 1 / rgamma(1L, shape = private$NG$shape1, rate = private$NG$rate1)
    }
  )

)
