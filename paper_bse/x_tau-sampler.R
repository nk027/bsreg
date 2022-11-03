
# Functions for lambda ~ Beta(1 + t, 1 + t), t ~ Gamma(alpha, beta)

f <- \(x, lambda, alpha, beta) { # PDF
  exp(x * log(lambda - lambda^2) + log(x) * (alpha - 1) - beta * x)
  (lambda - lambda^2)^x * x^(alpha-1) * exp(-beta * x)
}

g <- \(x, lambda, alpha, beta) { # CDF
  -x^alpha * (x * (beta - log(lambda - lambda^2)))^-alpha *
    expint::gammainc(alpha, x * (beta - log(lambda - lambda^2)))
}

# Clean pdf
dtau <- \(x, lambda, alpha, beta, log = FALSE) {
  out <- x * log(lambda - lambda^2) + log(x) * (alpha - 1) - beta * x +
    alpha * log(beta) - lgamma(alpha) # Helps with sampling
  if(log) out else exp(out)
  # (lambda - lambda^2)^x * x^(alpha - 1) * exp(-beta * x)
}

# Rejection sampler, using a Gamma proposal
rbg_rs <- \(n, lambda = .5, alpha = 1, beta = .5) {
  out <- numeric(n)
  i <- j <- 1
  # The Gamma density we draw from should be encompassing already
  m <- 1
  # exp(dtau(1e-32, lambda, alpha, beta, log = TRUE) - dgamma(1e-32, alpha, beta, log = TRUE))
  while (i <= n) {
    x <- rgamma(1L, alpha, beta)
    p_acc <- exp(dtau(x, lambda, alpha, beta, log = TRUE) - dgamma(x, alpha, beta, log = TRUE) - log(m))
    if(runif(1L) < p_acc) {
      out[i] <- x
      i <- i + 1
    }
    j <- j + 1
  }
  cat("Tries: ", j, " (", round(j / n, 2), " draws per sample)", sep = "")
  return(out)
}

# Test the number of samples needed
v <- c(.001, .01, .1, .5, 1)
for(i in v) {
  for(j in v) {
    print(paste0("a = ", i, " b =", j))
    rbg_rs(100, alpha = i, beta = j)
  }
}

# Compare densities
lambda <- .5
alpha <- .5
beta <- .1
x <- seq(1e-12, 100, len = 1000)

op <- par(mfrow = c(2, 2))
dtau(x, lambda, 2, 1, log = TRUE) |> plot(x = x, type = "l")
dgamma(x, 2, 1, log = TRUE) |> lines(x = x, lty = 2)

dtau(x, lambda, 1, .1, log = TRUE) |> plot(x = x, type = "l")
dgamma(x, 1, .1, log = TRUE) |> lines(x = x, lty = 2)

dtau(x, lambda, .1, .01, log = TRUE) |> plot(x = x, type = "l")
dgamma(x, .1, .01, log = TRUE) |> lines(x = x, lty = 2)

dtau(x, lambda, .01, .001, log = TRUE) |> plot(x = x, type = "l")
dgamma(x, .01, .001, log = TRUE) |> lines(x = x, lty = 2)
par(op)

# Compare draws and densities
op <- par(mfrow = c(2, 2))
draws1 <- rbg_rs(100000, lambda, .5, .1) |> density(from = 0)
plot(draws1$x, draws1$y, type = "l")
dtau(draws1$x, lambda, .5, .1) |> lines(x = draws1$x, lty = 2)
draws2 <- rbg_rs(100000, lambda, .5, .1) |> density(from = 0)
plot(draws2$x, log(draws2$y), type = "l")
dtau(draws2$x, lambda, .5, .1, log = TRUE) |> lines(x = draws2$x, lty = 2)
draws3 <- rbg_rs(100000, lambda, .1, .01) |> density(from = 0)
plot(draws3$x, draws3$y, type = "l")
dtau(draws3$x, lambda, .1, .01) |> lines(x = draws3$x, lty = 2)
draws4 <- rbg_rs(100000, lambda, .1, .01) |> density(from = 0)
plot(draws4$x, log(draws4$y), type = "l")
dtau(draws4$x, lambda, .1, .01, log = TRUE) |> lines(x = draws4$x, lty = 2)
par(op)


draws1 <- rbg_rs(10000, .9, .5, .1) |> density(from = 0)
plot(draws1$x, log(draws1$y), type = "l")
(dgamma(draws1$x, .5, .1, log = TRUE)) |> lines(x = draws1$x, lty = 2) # Plain
(2 * dgamma(draws1$x, .5, .1, log = TRUE)) |> lines(x = draws1$x, lty = 3)
(2 * (.9 - .9^2)^draws1$x * dgamma(draws1$x, .5, .1, log = TRUE)) |> # Scaled
  lines(x = draws1$x, lty = 2, col = 2)
(dtau(draws1$x, .9, .5, .1, log = TRUE)) |> lines(x = draws1$x, lty = 4)

f(x, lambda, alpha, beta) |> log() |> plot(x = x, type = "l")
f(x, lambda + .2, alpha, beta) |> log() |> lines(x = x, lty = 2, col = 2)
f(x, lambda - .2, alpha, beta) |> log() |> lines(x = x, lty = 3, col = 2)
f(x, lambda, alpha, beta / 10) |> log() |> lines(x = x, lty = 2, col = 3)
f(x, lambda, alpha, beta * 10) |> log() |> lines(x = x, lty = 3, col = 3)
f(x, lambda, alpha / 10, beta) |> log() |> lines(x = x, lty = 2, col = 4)
f(x, lambda, alpha * 2, beta) |> log() |> lines(x = x, lty = 3, col = 4)

g(x, lambda, alpha, beta) |> plot(x = x, type = "l")
g(x, lambda + .2, alpha, beta) |> lines(x = x, lty = 2, col = 2)
g(x, lambda - .2, alpha, beta) |> lines(x = x, lty = 3, col = 2)
g(x, lambda, alpha, beta / 10) |> lines(x = x, lty = 3, col = 3)
g(x, lambda, alpha, beta * 10) |> lines(x = x, lty = 3, col = 3)

# Limit to (0, 1)
x <- seq(1e-16, 2, len = 50)
kappa <- g(1e-16, lambda, alpha, beta)
(-g(x, lambda, alpha, beta) / kappa + 1) |> plot(x = x, type = "l")
