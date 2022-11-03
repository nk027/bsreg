
# Basic slice samplers ---

slice <- function(n, init_theta, target, A) {
  u <- theta <- rep(NA, n)
  theta[1] <- init_theta
  u[1] <- runif(1, 0, target(theta[1])) # This never actually gets used
  for (i in 2:n) {
    u[i] <- runif(1, 0, target(theta[i - 1]))
    endpoints <- A(u[i], theta[i - 1]) # The second argument is used in the second example
    theta[i] <- runif(1, endpoints[1], endpoints[2])
  }
  return(list(theta = theta, u = u))
}

set.seed(6)
A <- function(u, theta = NA) c(0, -log(u))
res <- slice(10, 0.1, dexp, A)

n = 5
y = rnorm(n,.2)
f = Vectorize(function(theta, y.=y) exp(sum(dnorm(y., theta, log=TRUE)) + dexp(abs(theta), log=TRUE)))
# Function to numerically find endpoints
A = function(u, xx, f.=f) {
  left_endpoint = uniroot(function(x) f.(x) - u, c(-10^10, xx))
  right_endpoint = uniroot(function(x) f.(x) - u, c( 10^10, xx))
c(left_endpoint$root, right_endpoint$root)
}

slice2 <- function(n, init_theta = .5, like<, qprior) {
  u <- theta <- rep(NA, n)
  theta[1] <- init_theta
  u[1] <- runif(1, 0, like(theta[1]))
  for (i in 2:n) {
    u[i] <- runif(1, 0, like(theta[i - 1]))
    success <- FALSE
    endpoints <- 0:1
    while (!success) {
      # Inverse CDF
      up <- runif(1, endpoints[1], endpoints[2])
      theta[i] <- qprior(up)
      if (u[i] < like(theta[i])) {
        success <- TRUE
      } else {
        # Updated endpoints when proposed value is rejected
        if (theta[i] > theta[i - 1]) {
          endpoints[2] <- up
        }
        if (theta[i] < theta[i - 1]) {
          endpoints[1] <- up
        }
      }
    }
  }
  return(list(theta = theta, u = u))
}


# Target function ---

N <- 100

xy <- cbind(runif(N), runif(N))
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
Psi <- function(x) {
  W <- 1 / dist ^ x
  W / max(eigen(W)$values)
}

X <- rmatrix(N, 5)
y <- solve(diag(N) - .5 * Psi(3), X %*% 1:5 + rnorm(N))

rss <- \(x) {
    prec_ch <- chol(crossprod(X) / 1)
    b0 <- backsolve(prec_ch, forwardsolve(prec_ch, crossprod(X, y) / 1, upper.tri = TRUE, transpose = TRUE))
    b1 <- backsolve(prec_ch, forwardsolve(prec_ch, crossprod(X, Psi(x) %*% y / 1), upper.tri = TRUE, transpose = TRUE))
    e0 <- y - X %*% b0
    e1 <- Psi(x) %*% y - X %*% b1
    e0e0 <- sum(e0^2)
    e1e0 <- sum(e1 * e0)
    e1e1 <- sum(e1^2)
    rss <- ((e0e0) - (2 * .5 * e1e0) + (.5^2 * e1e1))
    return((NROW(X) - NCOL(X)) / 2 * log(rss))
}
ldet <- \(x) {
  determinant(diag(N) - Psi(x), log = TRUE)[[1]]
}
prior <- \(x) {
  dinvgamma(x, .5, 4, log = TRUE)
}

p1 <- matrix(NA_real_, 500, 5,
  dimnames = list(NULL, c("delta", "rss", "ldet", "priors", "posterior")))
p1[, 1] <- seq(.01, 15, len = nrow(p1))

for(i in seq(nrow(p1))) {
  p1[i, 2] <- rss(p1[i, 1])
  p1[i, 3] <- ldet(p1[i, 1])
  p1[i, 4] <- prior(p1[i, 1])
  p1[i, 5] <- p1[i, 3] - p1[i, 2] + p1[i, 4]
}
plot.ts(p1[, -1], main = "Posterior components")
