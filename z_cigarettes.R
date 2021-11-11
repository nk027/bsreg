
devtools::load_all()

df <- readxl::read_excel("~/repos/bse/data/cigarette+2var.xls")
contig_data <- readxl::read_excel("~/repos/bse/data/Spat-Sym-US.xls", col_names = FALSE)
xy_data <- readxl::read_excel("~/repos/bse/data/cigar_states.xls", col_names = TRUE)

# Utah has a wrong longitude
xy_data_fixed <- xy_data
xy_data_fixed[40, 3] <- 112

x <- as.matrix(df[, c("logc", "logp", "logy")])
y <- x[, 1]
X <- cbind(1, x[, -1])
# Add TFE
n_time <- length(unique(df$year))
X <- cbind(X, kronecker(diag(n_time), matrix(1, nrow(x) / n_time))[, -1])
# Add IFE
X <- cbind(X, kronecker(matrix(1, n_time), diag(nrow(x) / n_time))[, -1])

# Use Frisch-Waugh-Lovell theorem to get rid of fixed effects
# variables <- 1:3
# Q_fwl <- qr.Q(qr(X[, -variables, drop = FALSE]))
# y <- y - Q_fwl %*% crossprod(Q_fwl, y)
# X <- X[, variables, drop = FALSE] - Q_fwl %*% # We don't care about the intercept
#   crossprod(Q_fwl, X[, variables, drop = FALSE])

# Connectivity
W <- as.matrix(contig_data)
W <- kronecker(diag(n_time), W / rowSums(W))
# W <- kronecker(diag(n_time), W / max(Re(eigen(W)$values)))
dist <- as.matrix(dist(xy_data))
dist <- as.matrix(dist(xy_data_fixed))
diag(dist) <- Inf
Psi <- function(x) {
  W <- 1 / dist ^ x
  W <- W / max(Re(eigen(W[1:46, 1:46])$values))
  kronecker(diag(n_time), W)
}
# W <- Psi(3)

# Reproduce -----

# Row-stochastic binary contiguity ---

X_LX <- cbind(X, W %*% X[, 2:3]) # Easier for lm and spatialreg

out_lm <- blm(y ~ X - 1, n_save = 10000)
print(out_lm)
print(lm(y ~ X - 1))

out_slx <- bslx(y ~ X - 1, W = W, X_SLX = X[, 2:3], n_save = 10000)
print(out_slx)
print(lm(y ~ X_LX - 1))

out_sar <- bsar(y ~ X - 1, W = W, n_save = 10000,
  ldet_SAR = list(grid = FALSE, reps = n_time))
print(out_sar)
spatialreg::lagsarlm(y ~ X - 1, listw = spdep::mat2listw(W))

out_sem <- bsem(y ~ X, W = W, n_save = 10000,
  ldet_SEM = list(grid = FALSE, reps = n_time))
print(out_sem)
spatialreg::errorsarlm(y ~ X - 1, listw = spdep::mat2listw(W))

out_sdm <- bsdm(y ~ X, W = W, X_SLX = X[, 2:3], n_save = 10000,
  ldet_SAR = list(grid = FALSE, reps = n_time))
print(out_sdm)
spatialreg::lagsarlm(y ~ X_LX - 1, listw = spdep::mat2listw(W))

out_sdem <- bsdem(y ~ X, W = W, X_SLX = X[, 2:3], n_save = 10000,
  ldet_SEM = list(grid = FALSE, reps = n_time))
print(out_sdem)
spatialreg::errorsarlm(y ~ X_LX - 1, listw = spdep::mat2listw(W))


# Inverse-distance decay ---

out_slxd2 <- bslx(y ~ X, W = Psi(2), X_SLX = X[, 2:3], n_save = 10000)
print(out_slxd2)

out_slxd3 <- bslx(y ~ X, W = Psi(3), X_SLX = X[, 2:3], n_save = 10000)
print(out_slxd3)
summary(lm(y ~ cbind(X, W %*% X[, 2:3])))

out_slxdx <- bslx(y ~ X, W = Psi, X_SLX = X[, 2:3], n_save = 10000, n_burn = 5000,
  options = set_options(SLX = set_SLX(delta = 3, delta_scale = 0.05)))
print(out_slxdx)

save.image()
