
# Load packages ---
library("bsreg")

# Settings ---
n_save <- 50000L
n_burn <- 10000L

# Estimate ---
out_slxd2 <- bslx(y ~ X, W = Psi(2), X_SLX = X_lag,
  n_save = n_save, n_burn = n_burn)

out_slxd3 <- bslx(y ~ X, W = Psi(3), X_SLX = X_lag,
  n_save = n_save, n_burn = n_burn)

out_slxd4 <- bslx(y ~ X, W = Psi(4), X_SLX = X_lag,
  n_save = n_save, n_burn = n_burn)

out_slxdx <- bslx(y ~ X, W = Psi, X_SLX = X_lag,
  n_save = n_save, n_burn = n_burn, options = set_options(
    SLX = set_SLX(delta = 3, delta_scale = 0.05, delta_a = 2, delta_b = 2)))

# Store results
save.image("paper_jse/cigarettes_dist.Rda")
