
ndraw <- 10000
nbeta <- 1000
x_vals <- seq(0, 1, length.out = nbeta + 2)[-c(1, nbeta + 2)]

x <- matrix(NA_real_, ndraw, nbeta)
out1 <- out2 <- out3 <- out4 <- list(x, x, x, x, x)

for(i in seq(ndraw)) { # Shape = 1 --- Exponential
  tau <- rgamma(1, 1, .01)
  out1[[1]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, .1)
  out1[[2]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, 1)
  out1[[3]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, 10)
  out1[[4]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, 100)
  out1[[5]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
}
for(i in seq(ndraw)) { # Rate = 1
  tau <- rgamma(1, 100, 1)
  out2[[1]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 10, 1)
  out2[[2]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 1, 1)
  out2[[3]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, .1, 1)
  out2[[4]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, .01, 1)
  out2[[5]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
}
for(i in seq(ndraw)) { # Shape = 10
  tau <- rgamma(1, 10, .01)
  out3[[1]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 10, .1)
  out3[[2]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 10, 1)
  out3[[3]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 10, 10)
  out3[[4]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, 10, 100)
  out3[[5]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
}
for(i in seq(ndraw)) { # Shape = .1
  tau <- rgamma(1, .1, .01)
  out4[[1]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, .1, .1)
  out4[[2]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, .1, 1)
  out4[[3]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, .1, 10)
  out4[[4]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
  tau <- rgamma(1, .1, 100)
  out4[[5]][i, ] <- dbeta(x_vals, 1 + tau, 1 + tau)
}

# Plot -----

cols <- viridisLite::viridis(4, option = "D", end = 0.8)
ltys <- c(3, 1, 5, 1)

# Default Beta(1 + tau, 1 + tau)
taus <- c(0.1, 1, 10, 100)

cairo_pdf("betas.pdf", height = 3.5, width = 7,
  family = "Noto Sans", bg = "transparent", pointsize = 8)
op <- par(mar = c(1, 3.5, 3, 0))
plot.new()
plot.window(xlim = c(-1, 1), ylim = c(0, dbeta(.5, 101, 101)))
axis(1, at = c(-1, 0, 1), labels = FALSE)
axis(3, at = round(qbeta(p = c(.025, .975), 101, 101) * 2 - 1, 2), pos = 12,
  col = cols[4])
axis(3, at = round(qbeta(p = c(.025, .975), 11, 11) * 2 - 1, 2), pos = 11.9,
  col = cols[3], lty = 2)
axis(3, at = round(qbeta(p = c(.025, .975), 2, 2) * 2 - 1, 2), pos = 11.8,
  col = cols[2])
axis(2, las = 1,
  at = round(c(0, dbeta(.5, 2, 2), dbeta(.5, 11, 11), dbeta(.5, 101, 101)), 2))
for(i in seq_along(taus))
  lines(dbeta(x_vals, 1 + taus[i], 1 + taus[i]), x = x_vals * 2 - 1,
    col = cols[i], lty = ltys[i], lwd = 2)
text(0.15, 11, labels = "τ = 100", cex = 1.5)
text(0.35, 2.5, labels = "τ = 10", cex = 1.5)
text(-.45, 1.8, labels = "τ = 1", cex = 1.5)
text(-.9, 1.4, labels = "τ = 0.1", cex = 1.5)
dev.off()

# Exponential on tau

bexp <- lapply(out1, function(x) colSums(x) / ndraw)[-5]
qbgamma <- \(x, i = 1) {
  c(x_vals[which(cumsum(x[[i]]) / sum(x[[i]]) > .025)[1]] * 2 - 1,
    x_vals[which(cumsum(x[[i]]) / sum(x[[i]]) > .975)[1]] * 2 - 1)
}

cairo_pdf("beta_expos.pdf", height = 3.5, width = 7,
  family = "Noto Sans", bg = "transparent", pointsize = 8)
op <- par(mar = c(1, 3.5, 3, 0))
plot.new()
plot.window(xlim = c(-1, 1), ylim = c(0, max(do.call(c, bexp))))
axis(1, at = c(-1, 0, 1), labels = FALSE)
axis(3, at = round(qbgamma(bexp, 1), 2), pos = 10.5,
  col = cols[4])
axis(3, at = round(qbgamma(bexp, 2), 2), pos = 10.4,
  col = cols[3], lty = 2)
axis(3, at = round(qbgamma(bexp, 3), 2), pos = 10.3,
  col = cols[2])
axis(2, las = 1, at = round(c(0, sapply(bexp, max)[-4]), 2))
for(i in seq_along(bexp))
  lines(bexp[[5 - i]], x = x_vals * 2 - 1, col = cols[i], lty = ltys[i], lwd = 2)
text(0.15, 9.5, labels = "β = 0.01", cex = 1.5)
text(0.35, 2.5, labels = "β = 0.1", cex = 1.5)
text(-.45, 1.8, labels = "β = 1", cex = 1.5)
text(-.9, 1.4, labels = "β = 10", cex = 1.5)
dev.off()

# Four variants

cairo_pdf("beta_gammas.pdf", height = 3.5, width = 7,
  family = "Noto Sans", bg = "transparent", pointsize = 8)
op <- par(mar = c(2, 3.5, 2, 0), mfrow = c(2, 2))
out <- list("τ ~ G(ɑ = 1, β)" = out1, "τ ~ G(a, β = 1)" = out2,
  "τ ~ G(ɑ = 10, β)" = out3, "τ ~ G(ɑ = 0.1, β)" = out4)
pos <- list(c(-.5, -.4, -.3), c(-.6, -.5, -.4), c(-1.8, -1.5, -1.2), c(0, .025, .05))
for(i in seq_along(out)) {
  x <- lapply(out[[i]], function(x) colSums(x) / ndraw)[-5]
  ymax <- max(do.call(c, x))
  plot.new()
  plot.window(xlim = c(-1, 1), ylim = c(0, ymax))
  title(names(out)[i])
  # axis(1, at = c(-1, 0, 1), labels = FALSE)
  axis(1, at = round(qbgamma(x, 1), 2), labels = c(round(qbgamma(x, 1), 2)[1], ""),
    pos = pos[[i]][1], col = cols[4])
  axis(1, at = round(qbgamma(x, 2), 2), labels = c("", round(qbgamma(x, 2), 2)[2]),
    pos = pos[[i]][2], col = cols[3], lty = 2)
  if(i != 4)
    axis(1, at = round(qbgamma(x, 3), 2), pos = pos[[i]][3], col = cols[2])
  axis(2, las = 1, at = round(c(0, sapply(x, max)[-4]), 2))
  for(i in seq_along(x))
    lines(x[[5 - i]], x = x_vals * 2 - 1, col = cols[i], lty = ltys[i], lwd = 2)
  # text(0.15, 9.5, labels = "β = 100", cex = 1.5)
  # text(0.35, 2.5, labels = "β = 10", cex = 1.5)
  # text(-.45, 1.8, labels = "β = 1", cex = 1.5)
  # text(-.9, 1.4, labels = "β = 0.1", cex = 1.5)
}
dev.off()
