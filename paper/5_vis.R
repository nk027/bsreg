
library("dplyr") # or dplyr

# Beta-Gamma shrinkage ---

n_draw <- 1000
n_beta <- 1000
x_vals <- seq(0, 1, length.out = n_beta + 2)[-c(1, n_beta + 2)]

gamma_coef <- \(mu, var) c("shape" = mu^2 / var, "rate" = 1 / (var / mu))
gamma_shape <- \(mu, var) gamma_coef(mu, var)["shape"]
gamma_rate <- \(mu, var) gamma_coef(mu, var)["rate"]

pars <- expand.grid(
  "mu" = c(0.01, 1, 10, 100),
  "var" = c(0.01, 0.1, 1, 10, 100)
) |> rowwise() |>
  mutate(shape = gamma_shape(mu, var), rate = gamma_rate(mu, var))

out <- vector("list", nrow(pars))
for(i in seq_along(out)) {
  out[[i]] <- matrix(NA_real_, n_draw, n_beta)
  for(j in seq_len(n_draw)) {
    tau <- rgamma(1L, pars[i, ][["shape"]], pars[i, ][["rate"]])
    out[[i]][j, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  }
}

op <- par(mfrow = c(5, 4), mar = c(2, 2, 2, 2))
for(i in seq_along(out)) {
  plot(x = (x_vals - .5) * 2, y = colMeans(out[[i]]), type = "l",
    main = paste0("mu = ", pars[[i, 1]], ", var = ", pars[[i, 2]]),
    # main = paste0("shape = ", pars[[i, 3]], ", rate = ", pars[[i, 4]]),
    ylim = c(0, 5)
    # ylim = range(apply(out[[i]], 2, quantile, c(0.01, 0.99)))
  )
  lines(x = (x_vals - .5) * 2, y = apply(out[[i]], 2, min), lty = 2, col = "darkred")
  lines(x = (x_vals - .5) * 2, y = apply(out[[i]], 2, max), lty = 2, col = "darkgreen")
  lines(x = (x_vals - .5) * 2, y = apply(out[[i]], 2, quantile, 0.9), lty = 3)
  lines(x = (x_vals - .5) * 2, y = apply(out[[i]], 2, quantile, 0.1), lty = 3)
}

op <- par(mfrow = c(5, 4), mar = c(2, 2, 2, 2))
for(i in seq_along(out)) {
  plot(x = (x_vals - .5) * 2, y = colMeans(out[[i]]), type = "l",
    main = paste0("mu = ", pars[[i, 1]], ", var = ", pars[[i, 2]]),
    # main = paste0("shape = ", pars[[i, 3]], ", rate = ", pars[[i, 4]]),
    ylim = c(0, 5)
    # ylim = range(apply(out[[i]], 2, quantile, c(0.01, 0.99)))
  )
  lines(x = (x_vals - .5) * 2, y = apply(out[[i]], 2, min), lty = 2, col = "darkred")
  lines(x = (x_vals - .5) * 2, y = apply(out[[i]], 2, max), lty = 2, col = "darkgreen")
  lines(x = (x_vals - .5) * 2, y = apply(out[[i]], 2, quantile, 0.9), lty = 3)
  lines(x = (x_vals - .5) * 2, y = apply(out[[i]], 2, quantile, 0.1), lty = 3)
}




  fo
mean <- shapes / rates
variance <- shapes / rates^2

x <- matrix(NA_real_, ndraw, nbeta)
out1 <- out2 <- out3 <- list(x, x, x, x, x)

for(i in seq(ndraw)) {
  tau <- rgamma(1, 1, .1)
  out1[[1]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, .5)
  out1[[2]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, 1)
  out1[[3]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, 1.5)
  out1[[4]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, 2)
  out1[[5]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
}
for(i in seq(ndraw)) {
  tau <- rgamma(1, .1, 1)
  out2[[1]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, .5, 1)
  out2[[2]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, 1)
  out2[[3]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1.5, 1)
  out2[[4]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 2, 1)
  out2[[5]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
}
for(i in seq(ndraw)) {
  tau <- rgamma(1, .5, 1)
  out3[[1]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, .5)
  out3[[2]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, 1)
  out3[[3]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 2, 1)
  out3[[4]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, 2)
  out3[[5]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
}

cols <- viridisLite::viridis(5, option = "D", end = 0.8)
ltys <- c(1, 2, 1, 2, 1)

cols <- c("#4B0055", "#006290", "#00AC8E", "#A6DA42")
ltys <- c(4, 3, 2, 1)

pdf("output/beta-gamma.pdf", width = 7, height = 3.4)
png("output/beta-exponential.png", width = 1000, height = 400, pointsize = 24)
op <- par(mfrow = c(1, 1), mar = c(2, 2, 0.1, 0.1), bg = "transparent")

out <- out1
out[[4]] <- NULL
y <- lapply(out, function(x) colSums(x) / ndraw)
plot(y[[1]], x = x_vals * 2 - 1, type = "l", col = cols[1], lty = ltys[1], lwd = 7,
  ylim = c(0, max(sapply(y, max))), ylab = "Density", xlab = "Value")
for(i in seq_along(out)[-1]) {
  lines(y[[i]], x = x_vals * 2 - 1, col = cols[i], lty = ltys[i], lwd = 7)
}
legend("topright", title = "rate", lwd = 7, bg = "transparent",
  legend = c("0.1", "0.5", "1.0", "2.0"), cex = 1.1,
  y.intersp = 0.75, x.intersp = 0.75,
  col = cols, lty = ltys)
dev.off()

out <- out1
y <- lapply(out, function(x) colSums(x) / ndraw)
plot(y[[1]], x = x_vals, type = "l", col = cols[1], lty = ltys[1], lwd = 2,
  ylim = c(0, max(sapply(y, max))), ylab = "Density", xlab = "Value")
for(i in seq_along(out)[-1]) {
  lines(y[[i]], x = x_vals, col = cols[i], lty = ltys[i], lwd = 2)
}
legend("topleft", title = "Rate", lwd = 2, bg = "white",
  legend = c(".1", ".5", "1", "1.5", "2"),
  col = cols, lty = ltys)

out <- out3
y <- lapply(out, function(x) colSums(x) / ndraw)
plot(y[[1]], x = x_vals, type = "l", col = cols[1], lty = ltys[1], lwd = 2,
  ylim = c(0, max(sapply(y, max))), ylab = "Density", xlab = "Value")
for(i in seq_along(out)[-1]) {
  lines(y[[i]], x = x_vals, col = cols[i], lty = ltys[i], lwd = 2)
}
legend("topleft", title = "Shape-Rate", lwd = 2, bg = "white",
  legend = c(".5-1", "1-.5", "1-1", "2-1", "1-2"),
  col = cols, lty = ltys)

