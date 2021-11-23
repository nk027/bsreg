
# Load packages ---
devtools::load_all()
library("spatialreg")

# Settings ---
n_save <- 25000L
n_burn <- 5000L
# Prepare the full matrix with lagged explanatories
X_SLX <- cbind(X, X_cont)
listw <- spdep::mat2listw(W)

# Estimate ---
out_blm <- blm(y ~ X - 1,
  n_save = n_save, n_burn = n_burn)
out_lm <- lm(y ~ X - 1)

out_bslx <- bslx(y ~ X - 1, W = W, X_SLX = X_lag,
  n_save = n_save, n_burn = n_burn)
out_slx <- lm(y ~  - 1)

out_bsar <- bsar(y ~ X - 1, W = W,
  n_save = n_save, n_burn = n_burn, ldet_SAR = list(reps = n_time))
out_sar <- lagsarlm(y ~ X - 1, listw = listw)

out_bsem <- bsem(y ~ X, W = W,
  n_save = n_save, n_burn = n_burn, ldet_SEM = list(reps = n_time))
out_sem <- errorsarlm(y ~ X - 1, listw = listw)

out_bsdm <- bsdm(y ~ X, W = W, X_SLX = X_lag,
  n_save = n_save, n_burn = n_burn, ldet_SAR = list(reps = n_time))
out_sdm <- lagsarlm(y ~ X_LX - 1, listw = listw)

out_bsdem <- bsdem(y ~ X, W = W, X_SLX = X_lag,
  n_save = n_save, n_burn = n_burn, ldet_SEM = list(reps = n_time))
out_sdem <- errorsarlm(y ~ X_LX - 1, listw = listw)

# Store results
save.image("paper/cigarettes_contig.Rda")
