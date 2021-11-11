
# Setup ---

set.seed(42)

# Parameters
N <- 500
beta <- c(2, -1)
M <- length(beta)
theta <- c(1)
sigma <- .5
lambda <- .4
# Connectivity
xy <- cbind(runif(N), runif(N))
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
Psi <- function(x) {
  W <- 1 / dist ^ x
  W / rowSums(W)
}
W <- Psi(2)
S <- diag(N) - lambda * W
D <- diag(N) - -lambda * W
# Data
X <- cbind(1, matrix(rnorm(N * (M - 1)), N, M - 1))
X_SLX <- X[, -1]
y_lm <- X %*% beta + rnorm(N, 0, sigma)
y_slx <- X %*% beta + W %*% X_SLX %*% theta + rnorm(N, 0, sigma)
y_sar <- solve(S, X %*% beta + rnorm(N, 0, sigma))
y_sem <- X %*% beta + solve(D, rnorm(N, 0, sigma))
y_sdm <- solve(S, X %*% beta + W %*% X_SLX %*% theta + rnorm(N, 0, sigma))
y_sdem <- X %*% beta + W %*% X_SLX %*% theta + solve(D, rnorm(N, 0, sigma))
y_sac <- solve(S, X %*% beta + solve(D, rnorm(N, 0, sigma)))

# Estimate ---

y <- y_sar

out_lm <- blm(y ~ X, n_save = 10000, n_burn = 1000)
print(out_lm)
c(beta, sigma^2)
summary(out_lm)
plot(out_lm)

lm(y_slx ~ X + W %*% X_SLX - 1)
out_lx <- bslx(y ~ X, W = W, n_save = 10000, n_burn = 1000)
print(out_lx)
c(beta, theta, sigma^2)
summary(out_lx)
plot(out_lx)

out_ar <- bsar(y ~ X, W = W, n_save = 10000, n_burn = 1000)
print(out_ar)
c(beta, sigma^2, lambda)
summary(out_ar)
plot(out_ar)

out_em <- bsem(y ~ X, W = W, n_save = 10000, n_burn = 1000)
print(out_em)
c(beta, sigma^2, lambda)
summary(out_em)
plot(out_em)

out_dm <- bsdm(y ~ X, W = W, n_save = 10000, n_burn = 1000)
print(out_dm)
c(beta, theta, sigma^2, lambda)
summary(out_dm)
plot(out_dm)

out_de <- bsdem(y ~ X, W = W, n_save = 10000, n_burn = 1000)
print(out_de)
c(beta, theta, sigma^2, lambda)
summary(out_de)
plot(out_de)


library("dplyr")
library("ggplot2")
library("ggdist")

# SAR is a hassle
# eff_dir[i, ] <- sum(diag(B)) / N * betas[i, index %in% c("alpha", "beta")] +
#   if(LX) {
#     c(0, sum(diag(B %*% W)) / N * betas[i, index == "theta"])
#   } else {0}
# eff_tot <- sum(B) / N * betas[i, index %in% c("alpha", "beta")] +
#   if(LX) {
#     c(0, sum(B %*% W) / N * betas[i, index == "theta"])
#   } else {0}

d_ar <- out_ar$draws %>% as_tibble()
direct <- vector("numeric", 5000)
for(i in seq_along(direct)) {
  S <- solve(diag(N) - d_ar$lambda_SAR[i] * W)
  direct[i] <- as.numeric(sum(diag(S)) / N * d_ar$beta2[i])
}
total <- d_ar$beta2[seq(1, 10000, 2)] / (1 - d_ar$lambda_SAR[seq(1, 10000, 2)])
d_ar$direct <- rep(direct, 2)
d_ar$indirect <- rep(total, 2) - d_ar$direct

d_ar <- d_ar %>% transmute(type = "SAR", total = direct + indirect, direct, indirect, sigma)
d_lm <- out_lm$draws %>% as_tibble() %>%
  transmute(type = "LM", total = beta2, direct = beta2, indirect = NA, sigma)
d_lx <- out_lx$draws %>% as_tibble() %>%
  transmute(type = "SLX", total = beta2 + beta3, direct = beta2, indirect = beta3, sigma)
d_em <- out_em$draws %>% as_tibble() %>%
  transmute(type = "SEM", total = beta2, direct = beta2, indirect = NA, sigma)

x <- rbind(d_lm, d_lx, d_ar, d_em)

# True effects
S <- solve(diag(N) - lambda * W)
total_true <- beta[2]
total_true <- beta[2] + theta[2]
total_true <- beta[2] / (1 - lambda)

direct_true <- beta[2]
direct_true <- beta[2]
direct_true <- beta[2] * sum(diag(S)) / N

ggplot(x) +
  stat_halfeye(aes(x = type, y = total)) +
  geom_hline(yintercept = total_true)
ggplot(x) +
  stat_halfeye(aes(x = type, y = direct)) +
  geom_hline(yintercept = direct_true)
ggplot(x) +
  stat_halfeye(aes(x = type, y = indirect)) +
  geom_hline(yintercept = total_true - direct_true)
