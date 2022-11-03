
# Simulate data ---
N <- 5

# Simulate connectivity
xy <- cbind(runif(N), runif(N))
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
Psi <- function(x) {
  W <- 1 / dist ^ x
  # W / max(eigen(W)$values)
  W
}
d <- diag(N)

v <- eigen(Psi(2))$values
r <- rowSums(Psi(2))

eigen(d - Psi(2))
# If we scale by the maximum
eigen(Psi(2) / v[1]) # Maximum is 1
eigen(d - Psi(2) / v[1]) # one will be degenerate (1 - 1)
eigen(d - .9 * Psi(2) / v[1]) # one will be degenerate (1 - 1)
# If we scale by the minimum
eigen(Psi(2) / v[5]) # Maximum is 1
eigen(d - Psi(2) / v[5]) # 1 will be degenerate (1 - 1)
eigen(d - .9 * Psi(2) / v[5]) # 1 will be degenerate (1 - 1)
eigen(d - .9 * Psi(2) / v[5]) # 1 will be degenerate (1 - 1)
eigen(d - .9 * Psi(2) / v[5]) # 1 will be degenerate (1 - 1)

eigen(Psi(2) / r)
eigen(d - Psi(2) / r)

# we rather scale by the spectral radius?

W <- Psi(2)
v <- eigen(W)$values

(s <- seq(1 / v[1], 1 / v[5], len = 200))
sapply(s, \(x) eigen(d - W * x)$values)
sapply(s, \(x) determinant(d - W * x)[[1]])
sapply(s, \(x) eigen(W * x)$values)

eigen(W * s[length(s)])$values
eigen(W * -s[1])$values

A <- B <- W / v[5]
ev1 <- ev2 <- matrix(NA_real_, 50, N)

for(i in seq(50)) {
  ev1[i, ] <- eigen(0.9 * B)$values
  ev2[i, ] <- eigen(d - 0.9 * B)$values
  B <- B %*% A
}

plot.ts(ev1)
plot.ts(ev2)

solve(d - .9 * A)