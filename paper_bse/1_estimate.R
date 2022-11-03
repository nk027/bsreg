
devtools::load_all()

# Simulate data ---
N <- 50

# Simulate connectivity
xy <- cbind(runif(N), runif(N))
# xy <- expand.grid(seq(0, 10), seq(0, 10))
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
Psi <- function(x) {
  W <- 1 / dist ^ x
  W / max(eigen(W)$values)
}

# Simulate a SAR
X <- rmatrix(N, 5)

lambda <- -2
y <- solve(diag(N) - lambda * Psi(1), X %*% 1:5 + rnorm(N, 0, .1))

# Simulate SLX
y <- X %*% 1:5 + Psi(1) %*% X[, 1:2] %*% c(10, 0) + rnorm(N, 0, .1)
y <- X %*% 1:5 + (Psi(1) / (1 + c(0, 1, rep(1e6, N - 2)))) %*% X[, 1:2] %*% c(10, 0) + rnorm(N, 0, .1)

# Estimation ---

class <- get_slx_class()
m <- class$new(priors = set_options(SLX = set_SLX(xi_prior = "gamma",
  xi_a = .1, xi_b = .01, xi_scale = 0.1))$priors)
m$setup(y = y, X = X, X_SLX = X[, 1:2], Psi = Psi(1))
# burn(m)
x <- sample(m, n_save = 5000)

plot.ts(x[, 1:10])
plot.ts(x[, 9:18])

# Estimation ---

class <- get_sar_class()
n_save <- 5000

# With delta
# m <- class$new(priors = set_options(SAR = set_SAR("bgamma",
#   lambda_a = 1, lambda_b = .01, delta_scale = 0.1))$priors)
# m$setup(y = y, X = X, Psi = Psi, ldet_SAR = list(grid = TRUE, reps = 1L,
#   i_lambda = c(-1 + 1e-12, 1 - 1e-12, 100L), i_delta = c(1e-12, 10, 20)))

# Without delta
m <- class$new(priors = set_options(SAR = set_SAR("beta",
  lambda_a = 1, lambda_b = 1, delta_scale = 0, lambda_min = -1.5))$priors)
m$setup(y = y, X = X, Psi = Psi(1))
burn(m, 500L)

s1 <- matrix(NA_real_, n_save, length(unlist(m$get_parameters)))
for(i in seq(n_save)) {
  m$sample()
  s1[i, ] <- unlist(m$get_parameters)
  if(i %% 100 == 0) cat(i / n_save * 100, "%\t", sep = "")
}
# plot.ts(s1)
summary(s1)
assign(paste0("s1_", N), s1)

m <- class$new(priors = set_options(SAR = set_SAR("bgamma",
  lambda_a = .5, lambda_b = .0005, delta_scale = 0))$priors)
m$setup(y = y, X = X, Psi = Psi(2))
burn(m, 500L)

s2 <- matrix(NA_real_, n_save, length(unlist(m$get_parameters)))
for(i in seq(n_save)) {
  m$sample()
  s2[i, ] <- unlist(m$get_parameters)
  if(i %% 100 == 0) cat(i / n_save * 100, "%\t", sep = "")
}
# plot.ts(s2)
summary(s2)
assign(paste0("s2_", N), s2)

m <- class$new(priors = set_options(SAR = set_SAR("bgamma",
  lambda_a = 1, lambda_b = .001, delta_scale = 0))$priors)
m$setup(y = y, X = X, Psi = Psi(2))
burn(m, 500L)

s3 <- matrix(NA_real_, n_save, length(unlist(m$get_parameters)))
for(i in seq(n_save)) {
  m$sample()
  s3[i, ] <- unlist(m$get_parameters)
  if(i %% 100 == 0) cat(i / n_save * 100, "%\t", sep = "")
}
# plot.ts(s3)
summary(s3)
assign(paste0("s3_", N), s3)

op <- par(mfrow = c(2, 2))
plot(s2[, 7])
points(s3[, 7], col = 3)
points(s1[, 7], col = "gray")
abline(h = lambda)

plot(density(s2[, 7]))
lines(density(s3[, 7]), col = 3)
lines(density(s1[, 7]), col = "gray")
abline(v = c(mean(s2[, 7]), mean(s3[, 7]), mean(s1[, 7])),
  col = c(1, 3, "gray"), lty = 2:4)
abline(v = lambda)

plot(s2[, 8])
points(s3[, 8], col = 3)
points(s1[, 8], col = "gray")

plot(density(s2[, 8]))
lines(density(s3[, 8]), col = 3)
lines(density(s1[, 8]), col = "gray")
par(op)


op <- par(mfrow = c(2, 2))
plot(density(s1_50[, 7]), main = "Beta")
lines(density(s2_50[, 7]), col = "gray")
abline(v = lambda)

plot(density(s1_100[, 7]))
lines(density(s2_100[, 7]), col = "gray")
abline(v = lambda)

plot(density(s1_200[, 7]))
lines(density(s2_200[, 7]), col = "gray")
abline(v = lambda)

plot(density(s1_500[, 7]))
lines(density(s2_500[, 7]), col = "gray")
abline(v = lambda)


set.seed(42)
devtools::load_all()

# Construct an inverse distance-decay matrix ---
xy <- cigarettes[cigarettes[["year"]] == 1980, c("longitude", "latitude")]
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
Psi <- function(x) {
  W <- dist ^ -x
  W <- W / max(eigen(W, symmetric = TRUE)$values)
  kronecker(diag(n_time), W)
}
delta <- 3 # Decay parameter
W_decay <- dist ^ -delta
W_scaled <- W_decay / max(eigen(W_decay, symmetric = TRUE)[["values"]])
n_time <- length(unique(cigarettes[["year"]]))
W <- kronecker(diag(n_time), W_scaled) # Repeated for every year

cigarettes$p <- cigarettes$price / cigarettes$cpi
cigarettes$i <- cigarettes$ndi / cigarettes$cpi

# Spatial Durbin model (Uniform lambda)
x1 <- bslx(log(sales) ~ log(p) + log(i) + as.factor(name) + as.factor(year),
  X_SLX = log(cigarettes[, c("p", "i")]) |> as.matrix(), W = W,
  data = cigarettes,
  n_save = 1000L, n_burn = 000L)

# Spatial Durbin model (Uniform lambda)
x2 <- bsdm(log(sales) ~ log(p) + log(i) + as.factor(name) + as.factor(year),
  X_SLX = log(cigarettes[, c("p", "i")]) |> as.matrix(), Psi_SLX = W,
  data = cigarettes,
  W = W, options = set_options(
    SAR = set_SAR(lambda_a = 1, lambda_b = 1)),
  n_save = 10000L, n_burn = 000L)

x2 <- bsar(log(sales) ~ log(p) + log(i) + as.factor(name) + as.factor(year),
  # X_SLX = log(cigarettes[, c("p", "i")]) |> as.matrix(), Psi_SLX = W,
  data = cigarettes,
  W = W, options = set_options(
    SAR = set_SAR(lambda_a = 1, lambda_b = 1)),
  n_save = 10000L, n_burn = 000L)

# Spatial Durbin model (Uniform lambda)
x3 <- bsdm(log(sales) ~ log(p) + log(i) + as.factor(name) + as.factor(year),
  X_SLX = log(cigarettes[, c("p", "i")]) |> as.matrix(), Psi_SLX = W,
  data = cigarettes, W = Psi,
  options = set_options(
    SAR = set_SAR(lambda_a = 1, lambda_b = 1, delta_scale = .2)),
  ldet_SAR = list(grid = FALSE, reps = 1L, i_lambda = c(-1, 1 - 1e-12, 100L), i_delta = c(1e-12, 10, 20)),
  n_save = 1000L, n_burn = 000L)

x3 <- bsar(log(sales) ~ log(p) + log(i) + as.factor(name) + as.factor(year),
  # X_SLX = log(cigarettes[, c("p", "i")]) |> as.matrix(), Psi_SLX = W,
  data = cigarettes, W = Psi,
  options = set_options(
    SAR = set_SAR(lambda_a = 1, lambda_b = 1, delta_scale = .2, delta = 1.3)),
  ldet_SAR = list(grid = FALSE, reps = 1L, i_lambda = c(-1, 1 - 1e-12, 100L), i_delta = c(1e-12, 10, 20)),
  n_save = 1000L, n_burn = 000L)

# Outputs ---
x3 <- bm(x3)

# Convergence
plot(x)
plot(as.mcmc(x))
coda::geweke.diag(as.mcmc(x))

# Analysis
print(x)
apply(x[[1]], 2, quantile, c(0.025, 0.5, 0.975))
coda::HPDinterval(as.mcmc(x))
plot(density(x[[1]][, "lambda_SAR"]))



set.seed(42)
devtools::load_all()
library("sf")

data <- readRDS("paper/data.rds") |> dplyr::filter(date %in% 2006:2017)
W_qu <- readRDS("~/Documents/70_work/projects/stats_sustain/weights_queen.rds")
n_time <- length(unique(data$date))
xy <- data |> dplyr::filter(date == 2010) |> st_centroid() |> st_coordinates()
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
Psi <- function(x) {
  W <- dist ^ -x
  W <- W / max(eigen(W, symmetric = TRUE)$values)
  kronecker(diag(n_time), W)
}
delta <- 2 # Decay parameter
W <- Psi(3)

X_SLX <- as.matrix(st_drop_geometry(data)[, c("forest_px_km2_lag",
  "pasture_px_km2_lag", "crop_px_km2_lag", "cattle_dens_lag_log",
  "soy_filled_lag", "pop_km2_lag_log", "spei_dry")])

# Spatial Durbin model (Uniform lambda)
x0 <- bslx(forest_ch_km2 ~ forest_px_km2_lag + pasture_px_km2_lag +
  crop_px_km2_lag + cattle_dens_lag_log + soy_filled_lag +
  pop_km2_lag_log + spei_dry + as.factor(date) + as.factor(name),
  X_SLX = X_SLX,
  W = kronecker(diag(n_time), W_qu),
  data = data,
  n_save = 5000L, n_burn = 1000L)

n <- c(colnames(x0$model$.__enclos_env__$self$X), "sigma", "delta")

# Spatial Durbin model (Uniform lambda)
x1 <- bslx(forest_ch_km2 ~ forest_px_km2_lag + pasture_px_km2_lag +
  crop_px_km2_lag + cattle_dens_lag_log + soy_filled_lag +
  pop_km2_lag_log + spei_dry + as.factor(date) + as.factor(name), W = W,
  X_SLX = X_SLX,
  data = data,
  n_save = 5000L, n_burn = 1000L)

# Spatial Durbin model (Uniform lambda)
x2 <- bslx(forest_ch_km2 ~ forest_px_km2_lag + pasture_px_km2_lag +
  crop_px_km2_lag + cattle_dens_lag_log + soy_filled_lag +
  pop_km2_lag_log + spei_dry + as.factor(date) + as.factor(name), W = Psi,
  X_SLX = X_SLX, data = data,
  options = set_options(SLX = set_SLX(delta_scale = .1)),
  n_save = 1000L, n_burn = 500L)

cbind(n, c(colMeans(x1$draws), delta), colMeans(x2$draw))[-9:-159, ]

n_draw <- 1000

x1a <- bsar(forest_ch_km2 ~ forest_px_km2_lag + pasture_px_km2_lag +
  crop_px_km2_lag + cattle_dens_lag_log + soy_filled_lag +
  pop_km2_lag_log + spei_dry + as.factor(date) + as.factor(name), W = Psi(delta),
    options = set_options(
    ldet_SAR = list(grid = FALSE, reps = n_time)),
  data = data, n_save = n_draw, n_burn = 1000L)

x1b <- bsar(forest_ch_km2 ~ forest_px_km2_lag + pasture_px_km2_lag +
  crop_px_km2_lag + cattle_dens_lag_log + soy_filled_lag +
  pop_km2_lag_log + spei_dry + as.factor(date) + as.factor(name), W = Psi(1.4),
  options = set_options(
    ldet_SAR = list(grid = FALSE, reps = n_time)),
  data = data, n_save = n_draw, n_burn = 1000L)

# Spatial Durbin model (Uniform lambda)
x2 <- bsar(forest_ch_km2 ~ forest_px_km2_lag + pasture_px_km2_lag +
  crop_px_km2_lag + cattle_dens_lag_log + soy_filled_lag +
  pop_km2_lag_log + spei_dry + as.factor(date) + as.factor(name), W = Psi,
  data = data, options = set_options(
    ldet_SAR = list(grid = TRUE, reps = n_time, i_lambda = c(-1, 1 - 1e-12, 100L), i_delta = c(1e-12, 5, 20)),
    SAR = set_SAR(delta_scale = .1)),
  n_save = n_draw, n_burn = 100L)

te_sar1a <- te_sar1b <- te_sar2 <- numeric(n_draw)

n_reg <- nrow(data) / n_time

for(i in seq(n_draw)) {
  Si <- qr.solve(diag(n_reg) - x1a$draws[i, "lambda_SAR"] * Psi(delta)[seq(n_reg), seq(n_reg)])
  te_sar1a[i] <- sum(Si) / n_reg * x1a$draws[i, "beta4"]
}
for(i in seq(n_draw)) {
  Si <- qr.solve(diag(n_reg) - x1b$draws[i, "lambda_SAR"] * Psi(1.4)[seq(n_reg), seq(n_reg)])
  te_sar1b[i] <- sum(Si) / n_reg * x1b$draws[i, "beta4"]
}
for(i in seq(n_draw)) {
  Si <- qr.solve(diag(n_reg) - x2$draws[i, "lambda_SAR"] * Psi(x2$draws[i, "delta_SAR"])[seq(n_reg), seq(n_reg)])
  te_sar2[i] <- sum(Si) / n_reg * x2$draws[i, "beta4"]
}

plot(density(te_sar2))
lines(density(te_sar1a), lty = 2)
lines(density(te_sar1b), lty = 3)

plot.ts(x1a$draws[, 1:10])
plot.ts(x1a$draws[, 160:161])
plot.ts(x1b$draws[, 1:10])
plot.ts(x1b$draws[, 160:161])
plot.ts(x2$draws[, 1:10])
plot.ts(x2$draws[, 160:162])

sum(qr.solve(diag(n_reg) - mean(x2$draws[, "lambda_SAR"]) * Psi(mean(x2$draws[i, "delta_SAR"]))[seq(n_reg), seq(n_reg)])) / n_reg *
  mean(x2$draws[, "beta4"])
sum(qr.solve(diag(n_reg) - mean(x1a$draws[, "lambda_SAR"]) * Psi(2)[seq(n_reg), seq(n_reg)])) / n_reg *
  mean(x1a$draws[, "beta4"])
sum(qr.solve(diag(n_reg) - mean(x1a$draws[, "lambda_SAR"]) * Psi(1.4)[seq(n_reg), seq(n_reg)])) / n_reg *
  mean(x1b$draws[, "beta4"])

(x2$draws[, "beta4"] / x2$draws[, "lambda_SAR"]) |> density() |> plot()
(x1b$draws[, "beta4"] / x1b$draws[, "lambda_SAR"]) |> density() |> lines(lty = 2)

a <- b <- numeric(1000)

p1 <- p2 <- matrix(NA_real_, 1000, 4)

for(i in seq(1000)) {
  j <- base::sample(20000, 1)
  p1[i, 4] <- sum(qr.solve(diag(n_reg) - x1b$draws[j, "lambda_SAR"] * Psi(1.4)[seq(n_reg), seq(n_reg)])) / n_reg
  p1[i, 2] <- x1b$draws[j, "beta4"]
  p1[i, 3] <- x1b$draws[j, "lambda_SAR"]
  p1[i, 1] <- p1[i, 4] * x1b$draws[j, "beta4"]
}
for(i in seq(1000)) {
  j <- base::sample(20000, 1)
  p2[i, 4] <- sum(qr.solve(diag(n_reg) - x2$draws[j, "lambda_SAR"] * Psi(x2$draws[i, "delta_SAR"])[seq(n_reg), seq(n_reg)])) / n_reg
  p2[i, 2] <- x2$draws[j, "beta4"]
  p2[i, 3] <- x2$draws[j, "lambda_SAR"]
  p2[i, 1] <- p2[i, 4] * x2$draws[j, "beta4"]
}

plot.ts(cbind(p1, p2))



f <- \(x, y) x * dist^(-y)

f(.5, 2) |> det()

xs <- seq(-1, 1, length = 100)
ys <- seq(1e-12, 2, length = 100)

ns <- matrix(NA_real_, 100, 100)

for(i in seq_along(xs))
  for(j in seq_along(ys))
    ns[i, j] <- norm(f(xs[i], ys[j]))
