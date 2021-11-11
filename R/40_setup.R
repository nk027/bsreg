
#' Build a Bayesian linear model
#'
#' @param y Numeric vector with the dependent variable.
#' @param X Numeric matrix with the explanatory variables.
#' @param options List with settings and prior information.
#'
#' @return Returns an object with the desired Bayesian model.
#'
#' @noRd
get_blm <- function(y, X, options = set_options(), ...) {

  class <- switch(options$type, Independent = NormalGamma,
    Conjugate = ConjugateNormalGamma, Shrinkage = ShrinkageNormalGamma, Horseshoe = Horseshoe)

  mdl <- class$new(priors = options$priors)
  mdl$setup(y = y, X = X, ...)

  return(mdl)
}


#' Build a Bayesian spatially lagged explanatories model
#'
#' @inheritParams get_blm
#' @param Psi Numeric matrix (or function to construct one) with the spatial connectivities.
#' @param X_SLX Numeric matrix with explanatory variables that should be lagged spatially.
#'
#' @return Returns an object with the desired Bayesian model.
#'
#' @noRd
get_bslx <- function(y, X, options = set_options(), Psi, X_SLX, ...) {

  class <- get_slx_class(parent = switch(options$type, Independent = NormalGamma,
    Conjugate = ConjugateNormalGamma, Shrinkage = ShrinkageNormalGamma, Horseshoe = Horseshoe))

  mdl <- class$new(priors = options$priors)
  mdl$setup(y = y, X = X, Psi_SLX = Psi, X_SLX = X_SLX, ...)

  return(mdl)
}


#' Build a Bayesian spatial autoregressive model
#'
#' @inheritParams get_blm
#' @param Psi Numeric matrix (or function to construct one) with the spatial connectivities.
#'
#' @return Returns an object with the desired Bayesian model.
#'
#' @noRd
get_bsar <- function(y, X, options = set_options(), Psi, ...) {

  class <- get_sar_class(parent = switch(options$type, Independent = NormalGamma,
    Conjugate = ConjugateNormalGamma, Shrinkage = ShrinkageNormalGamma, Horseshoe = Horseshoe))

  mdl <- class$new(priors = options$priors)
  mdl$setup(X = X, y = y, Psi_SAR = Psi, ...)

  return(mdl)
}


#' Build a Bayesian spatial error model
#'
#' @inheritParams get_blm
#' @param Psi Numeric matrix (or function to construct one) with the spatial connectivities.
#'
#' @return Returns an object with the desired Bayesian model.
#'
#' @noRd
get_bsem <- function(y, X, options = set_options(), Psi, ...) {

  class <- get_sem_class(parent = switch(options$type, Independent = NormalGamma,
    Conjugate = ConjugateNormalGamma, Shrinkage = ShrinkageNormalGamma, Horseshoe = Horseshoe))

  mdl <- class$new(priors = options$priors)
  mdl$setup(X = X, y = y, Psi_SEM = Psi, ...)

  return(mdl)
}


#' Build a Bayesian spatial Durbin model
#'
#' @inheritParams get_blm
#' @inheritParams get_bslx
#' @param Psi,Psi_SLX Numeric matrix (or function to construct one) with the spatial connectivities.
#'
#' @return Returns an object with the desired Bayesian model.
#'
#' @noRd
get_bsdm <- function(y, X, options = set_options(), X_SLX, Psi, Psi_SLX, ...) {

  class <- get_sar_class(parent = get_slx_class(parent = switch(options$type, Independent = NormalGamma,
    Conjugate = ConjugateNormalGamma, Shrinkage = ShrinkageNormalGamma, Horseshoe = Horseshoe)))

  mdl <- class$new(priors = options$priors)
  mdl$setup(X = X, y = y, X_SLX = X_SLX, Psi_SAR = Psi, Psi_SLX = if(missing(Psi_SLX)) Psi else Psi_SLX, ...)

  return(mdl)
}


#' Build a Bayesian spatial Durbin error model
#'
#' @inheritParams get_blm
#' @inheritParams get_bslx
#' @param Psi,Psi_SLX Numeric matrix (or function to construct one) with the spatial connectivities.
#'
#' @return Returns an object with the desired Bayesian model.
#'
#' @noRd
get_bsdem <- function(y, X, options = set_options(), X_SLX, Psi, Psi_SLX, ...) {

  class <- get_sem_class(parent = get_slx_class(parent = switch(options$type, Independent = NormalGamma,
    Conjugate = ConjugateNormalGamma, Shrinkage = ShrinkageNormalGamma, Horseshoe = Horseshoe)))

  mdl <- class$new(priors = options$priors)
  mdl$setup(X = X, y = y, X_SLX = X_SLX, Psi_SEM = Psi, Psi_SLX = if(missing(Psi_SLX)) Psi else Psi_SLX, ...)

  return(mdl)
}


#' Build a Bayesian stochastic volatility model
#'
#' @inheritParams get_blm
#'
#' @return Returns an object with the desired Bayesian model.
#'
#' @noRd
get_bsv <- function(y, X, options = set_options(), ...) {

  class <- get_sv_class(parent = switch(options$type, Independent = NormalGamma,
    Conjugate = ConjugateNormalGamma, Shrinkage = ShrinkageNormalGamma, Horseshoe = Horseshoe))

  mdl <- class$new(priors = options$priors)
  mdl$setup(y = y, X = X)

  return(mdl)
}
