
#' Settings to tune a Metropolis-Hastings step
#' 
#' @param adjust_burn Numeric scalar with the percentage of burn-in that should be used to tune the MH step.
#' @param acc_target Numeric vector with the lower and upper bound of the target acceptance rate for the MH step.
#' @param acc_change Numeric scalar with the percentage adjustment to the proposal scale for tuning. 
#' 
#' @return Returns a list with settings to tune the Metropolis-Hastings step of a Bayesian model.
#' @export
#' 
#' @examples
#' set_mh(0.5, c(0.1, 0.5), .05)
set_mh <- function(
  adjust_burn = 0.8,
  acc_target = c(0.20, 0.45),
  acc_change = 0.01) {

  structure(list(
    adjust_burn = num_check(adjust_burn, min = 0, max = 1,
      msg = "Please provide a valid length for the MH tuning period (in percent of burn-in) via 'adjust_burn'."),
    acc_target = vapply(acc_target, num_check, numeric(1L), min = 0, max = 1,
      msg = "Please provide a valid target range for the MH tuning (in percent of acceptance) via 'acc_target'."),
    acc_change = num_check(adjust_burn, min = 0, max = 1e6,
      msg = "Please provide a valid scale adjustment factor for the MH tuning period (in percent) via 'acc_change'.")
  ), class = "mh_settings")
}
