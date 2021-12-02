
devtools::load_all()

# Setup ---

set.seed(42)

# Parameters
N <- 500
beta <- c(2, -1)
M <- length(beta)
theta <- beta[-1] * .5
sigma <- .1
lambda <- .5
delta <- 1
# Connectivity
xy <- cbind(runif(N), runif(N))
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
Psi <- function(x) {
  W <- 1 / dist ^ x
  W / max(eigen(W)$values)
}
W <- Psi(delta)
S <- diag(N) - lambda * W
D <- diag(N) - -lambda * W
# Data
X <- cbind(1, matrix(rnorm(N * (M - 1)), N, M - 1))
X_SLX <- X[, -1]
# Regressors depend
y_lm <- X %*% beta + rnorm(N, 0, sigma)
y_slx <- X %*% beta + W %*% X_SLX %*% theta + rnorm(N, 0, sigma)
y_sar <- solve(S, X %*% beta + rnorm(N, 0, sigma))
y_sem <- X %*% beta + solve(D, rnorm(N, 0, sigma))
y_sdm <- solve(S, X %*% beta + W %*% X_SLX %*% theta + rnorm(N, 0, sigma))
y_sdem <- X %*% beta + W %*% X_SLX %*% theta + solve(D, rnorm(N, 0, sigma))
y_sac <- solve(S, X %*% beta + solve(D, rnorm(N, 0, sigma)))

# Estimate rows ---

y <- y_sar

x <- get_bsar(y, X, Psi = Psi,
  options = set_options(SAR = set_SAR(
    lambda = 0, lambda_a = 1.01, lambda_b = 1.01, lambda_scale = 0.1, lambda_min = -1, lambda_max = 1,
    delta = 2, delta_a = 2, delta_b = 1, delta_scale = 0.1, delta_min = 0.01, delta_max = 10)),
  ldet_SAR = list(grid = FALSE, reps = 1L, i_lambda = c(-1, 1 - 1e-12, 100L), i_delta = c(0.1, 5, 10)))

x$sample()
x$get_parameters

s <- sample(x, 1000)
for(i in seq(1000)) s <- rbind(s, sample(x, 100))
plot.ts(s)
summary(s[-1:-100, ])

pars <- colMeans(s[-1:-100, ])
e <- solve(diag(N) - pars["lambda_SAR"] * Psi(pars["delta_SAR"]), diag(N) * pars["beta2"])
e_tot <- sum(e) / N
(e_dir <- sum(diag(e)) / N)
(e_ind <- e_tot - e_dir)

e_true <- solve(diag(N) - lambda * Psi(delta), diag(N) * beta[2])
sum(diag(e_true)) / N
(sum(e_true) / N) - (sum(diag(e_true)) / N)
