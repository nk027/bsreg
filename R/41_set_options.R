
#' Set up Bayesian model priors and settings
#'
#' @param type Character scalar with the prior type for the nested linear model.
#' @param NG Settings for the Normal-Gamma prior (independent or conjugate). See \code{\link{set_NG}}.
#' @param SNG Settings for the Normal-Gamma shrinkage prior (Polson and Scott, 2010). See \code{\link{set_NG}}.
#' @param HS Settings for the Horseshoe shrinkage prior (Makalic and Schmidt, 2015). See \code{\link{set_NG}}.
#' @param SAR Settings for the spatial autoregressive setup. See \code{\link{set_SAR}}.
#' @param SLX Settings for the spatially lagged explanatory setup. See \code{\link{set_SAR}}. Note that settings for
#' the spatial term 'theta' are provided to \emph{NG} instead.
#' @param SEM Settings for the spatial error setup. See \code{\link{set_SAR}}.
#' @param SV Settings for the stochastic volatility setup. See \code{\link{set_SV}}.
#' @param ... Used to provide custom prior elements.
#'
#' @return Returns a list with priors and settings for a Bayesian model.
#' @export
#'
#' @examples
#' set_options("Shrinkage", SNG = set_SNG(lambda_a = 1, lambda_b = 1))
set_options <- function(
  type = c("Independent", "Conjugate", "Shrinkage", "Horseshoe"),
  NG = set_NG(), SNG = set_SNG(), HS = set_HS(),
  SAR = set_SAR(), SLX = set_SLX(), SEM = set_SEM(),
  SV = set_SV(), ...
  ) {

  type <- match.arg(type)

  priors <- list("NG" = NG, # Normal-Gamma (independent or conjugate)
    "SNG" = SNG, "HS" = HS, # Shrinkage Normal-Gamma or Horseshoe
    "SAR" = SAR, "SLX" = SLX, "SEM" = SEM, "SV" = SV, ...)

  structure(list(
    "type" = type, "priors" = priors
  ), class = "priors")
}


#' Set up a Normal-Gamma prior
#'
#' @param mu Numeric scalar or vector with the prior mean of 'beta'.
#' @param precision Numeric scalar or matrix with the prior precision of 'beta'. Not used for shrinkage priors.
#' @param shape,rate Numeric scalars with the prior shape and rate of 'sigma'.
#' @param lambda_a,lambda_b Numeric scalars with the prior shape and rate of 'lambda'.
#' @param theta_scale Numeric scalar with the proposal scale of 'theta'. Defaults to zero for a fixed value.
#' @param theta_a Numeric scalar with the prior rate of 'theta'.
#' @param lambda,tau,theta,zeta,nu,beta,sigma Numerics with starting values for the respective parameter.
#'
#' @return Returns a list with priors and settings.
#' @export
set_NG <- function(mu = 0, precision = 1e-8, shape = 0.01, rate = 0.01, beta = NULL, sigma = NULL) {

  structure(list(
    mu = vapply(mu, num_check, numeric(1L), min = -Inf, max = Inf,
      msg = "Please provide a valid value for the prior mean via 'mu'."),
    precision = vapply(precision, num_check, numeric(1L), min = 0, max = Inf,
      msg = "Please provide a valid prior variance via 'precision'."),
    shape = num_check(shape, min = 1e-12, max = Inf, msg = "Please provide a valid prior shape via 'shape'."),
    rate = num_check(rate, min = 1e-12, max = Inf, msg = "Please provide a valid prior rate via 'rate'."),
    beta = beta, sigma = sigma
  ), class = "prior_NG")
}

#' @export
#' @rdname set_NG
set_SNG <- function(lambda_a = 0.01, lambda_b = 0.01, theta_scale = 0, theta_a = 1,
  lambda = 1, tau = 10, theta = 0.1) {

  structure(list(
    lambda_a = num_check(lambda_a, min = 1e-12, max = Inf,
      msg = "Please provide a valid shape for lambda (Normal-Gamma) via 'lambda_a'."),
    lambda_b = num_check(lambda_b, min = 1e-12, max = Inf,
      msg = "Please provide a valid rate for lambda (Normal-Gamma) via 'lambda_b'."),
    theta_scale = num_check(theta_scale, min = 0, max = Inf,
      msg = "Please provide a valid proposal scale for theta (Normal-Gamma) via 'theta_scale'."),
    theta_a = num_check(theta_a, min = 0, max = Inf,
      msg = "Please provide a valid rate for theta (Normal-Gamma) via 'theta_a'."),
    lambda = num_check(lambda, min = 1e-12, max = Inf,
      msg = "Please provide a valid starting value for 'lambda' (Normal-Gamma)."),
    tau = vapply(tau, num_check, numeric(1L), min = 1e-12, max = Inf,
      msg = "Please provide a valid starting value for 'tau' (Normal-Gamma)."),
    theta = num_check(theta, min = 1e-12, max = Inf,
      msg = "Please provide a valid starting value for 'theta' (Normal-Gamma).")
  ), class = "prior_SNG")
}

#' @export
#' @rdname set_NG
set_HS <- function(lambda = 1, tau = 1, zeta = 1, nu = 1) {

  structure(list(
    lambda = vapply(lambda, num_check, numeric(1L), min = 1e-12, max = Inf,
      msg = "Please provide a valid starting value for 'lambda' (Horseshoe)."),
    tau = num_check(tau, min = 1e-12, max = Inf,
      msg = "Please provide a valid starting value for 'tau' (Horseshoe)."),
    zeta = num_check(zeta, min = 1e-12, max = Inf,
      msg = "Please provide a valid starting value for 'zeta' (Horseshoe)."),
    nu = vapply(nu, num_check, numeric(1L), min = 1e-12, max = Inf,
      msg = "Please provide a valid starting value for 'nu' (Horseshoe).")
  ), class = "prior_HS")
}


#' Set up a spatial prior
#'
#' @param lambda_a,lambda_b Numeric scalars with the prior shapes of the connectivity strength 'lambda'.
#' @param lambda_scale Numeric scalar with the proposal scale of 'lambda'.
#' @param lambda_min,lambda_max Numeric scalars with upper and lower bounds for 'lambda'.
#' @param delta_a,delta_b Numeric scalars with the prior shapes of the connectivity parameter 'delta'.
#' @param delta_scale Numeric scalar with the proposal scale of 'delta'. Defaults to zero for a fixed value.
#' @param delta_min,delta_max Numeric scalars with upper and lower bounds for 'delta'.
#' @param lambda,delta Numerics with starting values for the respective parameter.
#'
#' @return Returns a list with priors and settings.
#' @export
set_SAR <- function(
  lambda_a = 1.01, lambda_b = 1.01, lambda = 0, lambda_scale = 0.1, lambda_min = -1, lambda_max = 1 - 1e-12,
  delta_a = 1.01, delta_b = 1.01, delta = 1, delta_scale = 0, delta_min = 1e-12, delta_max = Inf) {

  structure(list(
    lambda_a = num_check(lambda_a, min = 1e-12, max = Inf,
      msg = "Please provide a valid shape for lambda (spatial) via 'lambda_a'."),
    lambda_b = num_check(lambda_b, min = 1e-12, max = Inf,
      msg = "Please provide a valid rate for lambda (spatial) via 'lambda_b'."),
    lambda_min = num_check(lambda_min, min = -Inf, max = Inf,
      msg = "Please provide a valid lower bound for lambda (spatial) via 'lambda_min'."),
    lambda_max = num_check(lambda_max, min = -Inf, max = Inf,
      msg = "Please provide a valid upper bound for lambda (spatial) via 'lambda_max'."),
    lambda = num_check(lambda, min = lambda_min, max = lambda_max,
      msg = "Please provide a valid starting value for lambda (spatial) via 'lambda'."),
    lambda_scale = num_check(lambda_scale, min = 1e-12, max = Inf,
      msg = "Please provide a valid proposal scale for lambda (spatial) via 'lambda_scale'."),
    delta_a = num_check(delta_a, min = 1e-12, max = Inf,
      msg = "Please provide a valid shape for delta (spatial) via 'delta_a'."),
    delta_b = num_check(delta_b, min = 1e-12, max = Inf,
      msg = "Please provide a valid rate for delta (spatial) via 'delta_b'."),
    delta_min = num_check(delta_min, min = -Inf, max = Inf,
      msg = "Please provide a valid lower bound for delta (spatial) via 'delta_min'."),
    delta_max = num_check(delta_max, min = -Inf, max = Inf,
      msg = "Please provide a valid upper bound for delta (spatial) via 'delta_max'."),
    delta = num_check(delta, min = delta_min, max = delta_max,
      msg = "Please provide a valid starting value for delta (spatial) via 'delta'."),
    delta_scale = num_check(delta_scale, min = 0, max = Inf,
      msg = "Please provide a valid proposal scale for delta (spatial) via 'delta_scale'.")
  ), class = "prior_SAR")
}
#' @export
#' @rdname set_SAR
set_SLX <- set_SAR
#' @export
#' @rdname set_SAR
set_SEM <- set_SAR


#' Set up a volatility prior
#'
#' @param priors Prior settings from \code{\link[stochvol]{specify_priors}}.
#' @param mu,phi,sigma,nu,rho,beta,latent0 Numerics with starting values for the respective parameter.
#'
#' @return Returns a list with priors and settings.
#' @export
set_SV <- function(
  priors, mu = 0, phi = 0.5, sigma = 1, nu = Inf, rho = 0, beta = 0, latent0 = 0) {

  if(missing(priors)) {
    has_package("stochvol")
    priors <- stochvol::specify_priors()
  }

  structure(list(
    priors = priors, parameters = list(
      mu = num_check(mu, min = -Inf, max = Inf, msg = "Please provide a valid starting value for 'mu' (SV)."),
      phi = num_check(phi, min = -1, max = 1, msg = "Please provide a valid starting value for 'phi' (SV)."),
      sigma = num_check(sigma, min = 0, max = Inf, msg = "Please provide a valid starting value for 'sigma' (SV)."),
      nu = num_check(nu, min = 2, max = Inf, msg = "Please provide a valid starting value for 'nu' (SV)."),
      rho = num_check(rho, min = -1, max = 1, msg = "Please provide a valid starting value for 'rho' (SV)."),
      beta = vapply(beta, num_check, numeric(1L), min = -Inf, max = Inf,
        msg = "Please provide a valid starting value for 'beta' (SV)."),
      latent0 = num_check(latent0, min = -Inf, max = Inf,
        msg = "Please provide a valid starting value for the latent variable (SV) via 'latent0'.")
  )), class = "prior_SV")
}
