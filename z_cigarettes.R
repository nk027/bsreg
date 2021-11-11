
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
W <- kronecker(diag(n_time), W / max(Re(eigen(W)$values)))
dist <- as.matrix(dist(xy_data))
dist <- as.matrix(dist(xy_data_fixed))
diag(dist) <- Inf
Psi <- function(x) {
  W <- 1 / dist ^ x
  W <- W / max(Re(eigen(W[1:46, 1:46])$values))
  kronecker(diag(n_time), W)
}
W <- Psi(3)

# Reproduce -----

# Row-stochastic binary contiguity ---

out_lm <- blm(y ~ X - 1)
summary(out_lm)
summary(lm(y ~ X))

out_lx <- bslx(y ~ X - 1, W = W, SLX = X[, 2:3])
summary(out_lx)
summary(lm(y ~ cbind(X, W %*% X[, 2:3])))

out_ar <- bsar(y ~ X - 1, W = W, ldet_SAR = list(grid = FALSE, reps = n_time), n_save = 10000)
print(out_ar)
summary(out_ar)
summary(lm(y ~ cbind(X, W %*% X[, 2:3]) - 1))

out1 <- sar_mh(cbind(y, X), W, f_ldet = function(rho, W) {
  log(det(diag(nrow(W[1:46, 1:46])) - rho * W[1:46, 1:46])) * n_time
}, n_draw = 10000)
summary(out1$rho)

out_em <- bsem(y ~ X, W = W, ldet_SAR = list(grid = FALSE, reps = n_time), n_save = 10000)
summary(out_em$draws[, c(1:3, 79)])
summary(lm(y ~ cbind(X, W %*% X[, 2:3])))


# Inverse-distance decay ---

out_lx <- bslx(y ~ X, W = Psi(2), X_SLX = X[, 2:3])
summary(out_lx$draws[, 1:3])
summary(lm(y ~ cbind(X, W %*% X[, 2:3])))

out_lx <- bslx(y ~ X, W = Psi, X_SLX = X[, 2:3], n_save = 5000,
  options = set_options(SLX = set_SLX(delta = 3, delta_scale = 0.05)))
summary(out_lx$draws[, 1:3])
summary(lm(y ~ cbind(X, W %*% X[, 2:3])))



# Benchmark ---
lw <- spdep::mat2listw(W)

mb({
  out_sr <- spatialreg::spBreg_lag(y ~ X, listw = lw,
    control = list(ndraw = 10000L, nomit = 1000L))
}, {
  out_ar <- bsar(y ~ X, W = W, n_save = 10000L, n_burn = 1000L)
}, times = 5)
