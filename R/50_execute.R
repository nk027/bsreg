
#' Obtain draws from a Bayesian model sampler
#'
#' @param x Bayesian model
#' @param n_save,n_burn Integer scalar with number of draws to save / burn.
#' @param mh Settings to tune the Metropolis-Hastings step. See \code{\link{set_mh}}.
#' @param verbose Logical scalar. Whether to print status updates.
#'
#' @return Returns a numeric matrix with stored draws. The Bayesian model is modified in place.
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
sample <- function(x, n_save = 1000L, n_burn = 0L, mh = set_mh(), verbose = TRUE) {

  if(n_burn > 0) {tune(x, n_burn = n_burn, mh = mh, verbose = verbose)}

  if(n_save > 0) {

    if(verbose) {
      timer <- Sys.time()
      cat("Starting sampler with", n_save, "draws.\n")
      pb <- txtProgressBar(min = 0, max = n_save, style = 3)
    }

    draw <- unlist(x$get_parameters)
    storage <- matrix(NA_real_, n_save, length(draw), dimnames = list(NULL, names(draw)))

    for(i in seq.int(n_save)) {
      x$sample()
      storage[i, ] <- unlist(x$get_parameters)

      if(verbose) {setTxtProgressBar(pb, i)}
    }

    if(verbose) {
      cat("\nFinished sampling after ", format(round(Sys.time() - timer, 2)), ".\n", sep = "")
      close(pb)
    }

  }

  return(storage)
}


#' Burn-in and tune a Bayesian model sampler
#'
#' @inheritParams sample
#'
#' @return Modifies the Bayesian model in place and returns it invisibly.
tune <- function(x, n_burn = 1000L, mh = set_mh(), verbose = TRUE) {

  if(verbose) {
    timer <- Sys.time()
    cat("Starting burn-in with", n_burn, "draws.\n")
    pb <- txtProgressBar(min = 0, max = n_burn, style = 3)
  }

  for(i in seq.int(n_burn)) {

    x$sample()

    if(i %% 10 == 0 && i <= mh$adjust_burn * n_burn) { # Every tenth step we consider tuning the MH step

      for(obj in x$MH) { # Loop over each Metropolis-Hastings object
        acc_rate <- obj$get_tuning
        if(acc_rate < mh$acc_target[1L]) { # Loosen
          obj$set_scale <- max(obj$get_scale * (1 - mh$acc_change), 1e-12)
        } else if(acc_rate > mh$acc_target[2L]) { # Tighten
          obj$set_scale <- obj$get_scale * (1 + mh$acc_change)
        }
      }
    }

    if(verbose) {setTxtProgressBar(pb, i)}
  }

  if(verbose) {
    cat("\nFinished burn-in after ", format(round(Sys.time() - timer, 2)), ".\n", sep = "")
    close(pb)
  }

  return(invisible(x))
}


#' @rdname tune
burn <- function(x, n_burn = 1000L, verbose = TRUE) {

  return(tune(x, n_burn = n_burn, mh = set_mh(adjust_burn = 0), verbose = verbose))
}
