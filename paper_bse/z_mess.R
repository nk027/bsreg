
n <- 10
xy <- cbind(runif(n), runif(n))
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
Psi <- function(x) {
  W <- 1 / dist ^ x
  return(W)
}
rs <- \(x, xi = rep(0, nrow(x))) x / (rowSums(x)) / (1 + xi)
ev <- \(x, xi = rep(0, nrow(x))) x / (max(eigen(x, only.values = TRUE)$values)) / (1 + xi)

Psi(2) |> eigen(only.values = TRUE)
Psi(10) |> eigen(only.values = TRUE)
Psi(.1) |> eigen(only.values = TRUE)

Psi(2) |> rs() |> eigen(only.values = TRUE)
Psi(10) |> rs() |> eigen(only.values = TRUE)
Psi(.1) |> rs() |> eigen(only.values = TRUE)

Psi(2) |> rs(xi = rgamma(n, 1, 1)) |> eigen(only.values = TRUE)
Psi(10) |> rs(xi = rgamma(n, 1, 1)) |> eigen(only.values = TRUE)
Psi(.1) |> rs(xi = rgamma(n, 1, 1)) |> eigen(only.values = TRUE)

Psi(2) |> ev() |> eigen(only.values = TRUE)
Psi(10) |> ev() |> eigen(only.values = TRUE)
Psi(.1) |> ev() |> eigen(only.values = TRUE)

Psi(2) |> ev(xi = rgamma(n, 1, 1)) |> eigen(only.values = TRUE)
Psi(10) |> ev(xi = rgamma(n, 1, 1)) |> eigen(only.values = TRUE)
Psi(.1) |> ev(xi = rgamma(n, 1, 1)) |> eigen(only.values = TRUE)

k <- 3
X <- rmatrix(n, k)
y <- rmatrix(n, 1)
xi <- c(10, rep(0, n - 1))

Psi(2) %*% X
(Psi(2) |> ev()) %*% X
(Psi(2) |> ev(xi = xi)) %*% X

# Extra step for SLX

# Update xi_i -->
# X <- (X, Psi * X) # Only (i, (K+1):(2K)) changes
# b <- inv(prec0 + XX / sigma, Xy / sigma + prec0 * mu0) # Updating formula for XX and Xy
# sum( (y - X * b)^2 )

XLX <- cbind(X, ev(Psi(2)) %*% X)
XLX_a <- cbind(X, ev(Psi(2), xi = xi) %*% X)

(XX <- crossprod(cbind(X, ev(Psi(2)) %*% X)))
(XX_a <- crossprod(cbind(X, ev(Psi(2), xi = xi) %*% X)))

(XL <- ev(Psi(2), xi = 0) %*% X)
(XL_a <- ev(Psi(2), xi = xi) %*% X)
ev(Psi(2), xi = xi)[1, ] %*% X

X_rm <- matrix(c(rep(0, k), XL[1, ] - XL_a[1, ]), 1) # The change to induce

-2 * X_rm[2] * XLX[1, 2] + X_rm[2]^2 # w1*² = (w1 - psi1)²
X[1, 1] * X_rm[2] # x1w1* = (x1 w1 - x1 psi1)

X[1, 1] %*% X_rm[4:6]
X[1, 2] %*% X_rm[4:6]
ur <- t(X[1, , drop = FALSE]) %*% X_rm[4:6] # upper right / lower left
bd <- 2 * XLX[1, 4:6] %*% X_rm[1, 4:6, drop = FALSE] - crossprod(X_rm[1, 4:6, drop = FALSE]) # block-diagonal

XLX[1, ] %*% X_rm[1, , drop = FALSE] - crossprod(X_rm[1, , drop = FALSE]) # bd still needs * 2

norm(XX[4:6, 4:6] - bd - XX_a[4:6, 4:6]) # Change in X'X
norm(XX[1:3, 4:6] - ur - XX_a[1:3, 4:6])

crossprod(XLX, y) - crossprod(XLX_a, y)

crossprod(X_rm[, 4:6, drop = FALSE], y[1]) # Change in X'y

# Change in beta
solve(XX, crossprod(XLX, y)) - solve(XX_a, crossprod(XLX_a, y))

# times change in x



XL[1, ] - (XL_a[1, ] + X_rm[seq(-k)])

XX - crossprod(t(c(rep(0, k), XL_a[1, ] - XL[1, ])),
  c(rep(0, k), XL_a[1, ] - XL[1, ]))
update_cp <- function(XY, X_rm, Y_rm = X_rm) {
  XY - crossprod(X_rm, Y_rm)
}
update_cp(XX, X_rm)
update_cp(XX, X_rm) - XX
update_cp(XX, X_rm) - XX_a
(update_cp(XX, XLX[1, , drop = FALSE]) - crossprod(XLX[-1, ])) |> norm()
update_cp(XX, cbind(0, XLX[1, -1, drop = FALSE])) - crossprod(XLX[-1, ])
