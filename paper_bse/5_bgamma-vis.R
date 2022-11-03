
# Exploratory plots -----

draw_bg <- \(n_tau = 10000, n_lambda = 1000, alpha = 1, beta = 1) {
  tau <- rgamma(n_tau, alpha, beta)
  lambda <- matrix(NA_real_, n_tau, n_lambda)
  for(i in seq(n_tau)) lambda[i, ] <- rbeta(n_lambda, 1 + tau[i], 1 + tau[i])
  list("tau" = tau, "lambda" = lambda, "alpha" = alpha, "beta" = beta)
}
plot_bg <- \(x, log = FALSE) {
  mu <- x$alpha / x$beta
  p <- round(sum(x$lambda > 0.4 & x$lambda < 0.6) / length(x$lambda), 3)
  dens <- x$lambda |> density(from = 1e-16, to = 1 - 1e-16, n = 500)
  if(isTRUE(log)) {
    dens$y <- log(dens$y)
  }
  plot(dens, main = paste0("lambda | tau ~ G(", x$alpha, ",", x$beta, ")"))
  lines(dens$x, dbeta(dens$x, 1 + mu, 1 + mu, log = log), col = 3, lty = 2)
  text(.3, dens$y[100] / 2, labels = paste0("Info = ", round(mu, 3)))
  text(.7, dens$y[100] / 2, labels = paste0("p(x≈0) = ", p))
}

x <- draw_bg()
plot_bg(x)

op <- par(mfrow = c(3, 3))
# Normal scale ---
draw_bg(alpha = .01, beta = 1e-16) |> plot_bg()
draw_bg(alpha = .01, beta = 1e-12) |> plot_bg()
draw_bg(alpha = .01, beta = 1e-8) |> plot_bg()

draw_bg(alpha = .1, beta = .001) |> plot_bg()
draw_bg(alpha = .1, beta = .01) |> plot_bg()
draw_bg(alpha = .1, beta = .1) |> plot_bg()

draw_bg(alpha = .5, beta = .005) |> plot_bg()
draw_bg(alpha = .5, beta = .05) |> plot_bg()
draw_bg(alpha = .5, beta = .5) |> plot_bg()

# Log scale ---
draw_bg(alpha = .01, beta = 1e-16) |> plot_bg(log = TRUE)
draw_bg(alpha = .01, beta = 1e-12) |> plot_bg(log = TRUE)
draw_bg(alpha = .01, beta = 1e-8) |> plot_bg(log = TRUE)

draw_bg(alpha = .1, beta = .001) |> plot_bg(log = TRUE)
draw_bg(alpha = .1, beta = .01) |> plot_bg(log = TRUE)
draw_bg(alpha = .1, beta = .1) |> plot_bg(log = TRUE)

draw_bg(alpha = .5, beta = .005) |> plot_bg(log = TRUE)
draw_bg(alpha = .5, beta = .05) |> plot_bg(log = TRUE)
draw_bg(alpha = .5, beta = .5) |> plot_bg(log = TRUE)

# Very spiky ---
draw_bg(alpha = .01, beta = 1e-16) |> plot_bg(log = TRUE)
draw_bg(alpha = .01, beta = 1e-12) |> plot_bg(log = TRUE)
draw_bg(alpha = .01, beta = 1e-8) |> plot_bg(log = TRUE)

draw_bg(alpha = .1, beta = 1e-16) |> plot_bg(log = TRUE)
draw_bg(alpha = .1, beta = 1e-12) |> plot_bg(log = TRUE)
draw_bg(alpha = .1, beta = 1e-8) |> plot_bg(log = TRUE)

draw_bg(alpha = 1, beta = 1e-16) |> plot_bg(log = TRUE)
draw_bg(alpha = 1, beta = 1e-12) |> plot_bg(log = TRUE)
draw_bg(alpha = 1, beta = 1e-8) |> plot_bg(log = TRUE)


# Visualisations -----

cols <- viridisLite::viridis(4, option = "D", end = 0.8)
ltys <- c(3, 6, 5, 1)
x_vals <- seq(0, 1, length.out = 2e4)

# Default Beta(1 + tau, 1 + tau) ---
taus <- c(0.1, 1, 10, 100)

cairo_pdf("paper/betas.pdf", height = 3.5, width = 7,
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

# Split Beta(1 + tau, 1 + tau) ---
x_vals_a <- seq(0, .5, length.out = 1e4)
x_vals_b <- seq(.5, 1, length.out = 1e4)

q5 <- c(dbeta(.5, 1.1, 1.1), dbeta(.5, 2, 2), dbeta(.5, 11, 11), dbeta(.5, 101, 101))
q5l <- c(sprintf("%.2f", q5[1:3]), sprintf("%.1f", q5[4]))
q9 <- c(dbeta(.95, 1.1, 1.1, log = TRUE), dbeta(.95, 2, 2, log = TRUE),
  dbeta(.95, 11, 11, log = TRUE), dbeta(.95, 101, 101, log = TRUE))
q9l <- c(sprintf("%.2f", q9[1:2]), sprintf("%.1f", q9[3]), sprintf("%.0f", q9[4]))
q9l[4] <- paste0("↓ ", q9l[4])

# Left with levels
cairo_pdf("paper/betas_left.pdf", height = 3.5, width = 3.6,
  family = "Noto Sans", bg = "transparent", pointsize = 8)
op <- par(mar = c(0, 3.5, 2, 3.5), mfrow = c(1, 1))
plot.new()
plot.window(xlim = c(-1, 0), ylim = c(0, dbeta(.5, 101, 101)))
axis(3, at = c(-1, 0, .9, 1), labels = TRUE)
axis(2, las = 1, at = q5, labels = q5l)
segments(-2, q5[4], 0, q5[4], col = "gray", lwd = 2, lty = 2)
for(i in seq_along(taus))
  lines(dbeta(x_vals_a, 1 + taus[i], 1 + taus[i]), x = x_vals_a * 2 - 1,
    col = cols[i], lty = ltys[i], lwd = 2)
text(-.2, 9, labels = "τ = 100", cex = 1.5, col = cols[4])
text(-.3, 2.8, labels = "τ = 10", cex = 1.5, col = cols[3])
text(-.45, 1.8, labels = "τ = 1", cex = 1.5, col = cols[2])
text(-.85, 1.4, labels = "τ = 0.1", cex = 1.5, col = cols[1])
abline(v = 0, col = "gray", lwd = 2)
dev.off()

# Right with logs
cairo_pdf("paper/betas_right.pdf", height = 3.5, width = 3.6,
  family = "Noto Sans", bg = "transparent", pointsize = 8)
op <- par(mar = c(0, 0, 2, 3.5), mfrow = c(1, 1))
plot.new()
plot.window(xlim = c(0, 1), yl= c(dbeta(.975, 11, 11, log = TRUE),
  dbeta(.5, 101, 101, log = TRUE)))
axis(3, at = c(-1, 0, .9, 1), labels = TRUE)
axis(4, las = 1, at = q9, labels = q9l)
axis(4, las = 1, at = dbeta(.975, 11, 11, log = TRUE), tick = FALSE,
  labels = formatC(dbeta(.95, 101, 101, log = TRUE), digits = 3))
segments(0.9, q9[2], 2, q9[2], col = "gray", lwd = 2, lty = 2)
for(i in seq_along(taus))
  lines(dbeta(x_vals_b, 1 + taus[i], 1 + taus[i], log = TRUE),
    x = x_vals_b * 2 - 1, col = cols[i], lty = ltys[i], lwd = 2)
abline(v = 0.9, col = "gray", lwd = 2)
dev.off()

# Bottom with intervals
qs <- list(c(.25, .75), c(.1, .9), c(.005, .995))

cairo_pdf("paper/betas_range.pdf", height = 1, width = 7,
  family = "Noto Sans", bg = "transparent", pointsize = 8)
op <- par(mar = c(2, .5, .5, .5))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(-.025, .425))
# abline(v = c(0.5, .95), col = "gray", lwd = 2)
axis(1, at = c(0, 1), labels = c(-1, 1))
axis(1, at = c(0.25, .75), labels = c(-0.5, 0.5))
rect(0, 0, 1, .075, col = "white", border = "white")
rect(0, .1, 1, .175, col = "white", border = "white")
rect(0, .2, 1, .275, col = "white", border = "white")
rect(0, .3, 1, .375, col = "white", border = "white")
col_ramp <- colorRampPalette(c(cols[4], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 101, 101)
  rect(pos[1], 0, pos[2], .075, col = col_ramp[j], border = col_ramp[j])
}
col_ramp <- colorRampPalette(c(cols[3], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 11, 11)
  rect(pos[1], .1, pos[2], .175, col = col_ramp[j], border = col_ramp[j])
}
col_ramp <- colorRampPalette(c(cols[2], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 2, 2)
  rect(pos[1], .2, pos[2], .275, col = col_ramp[j], border = col_ramp[j])
}
col_ramp <- colorRampPalette(c(cols[1], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 1.1, 1.1)
  rect(pos[1], .3, pos[2], .375, col = col_ramp[j], border = col_ramp[j])
}

segments(1, y0 = 0, y1 = .4)
segments(0, y0 = 0, y1 = .4)
text(qbeta(.25, 1.1, 1.1) + .025, .4035, cex = 1.2, "50%")
text(qbeta(.005, 1.1, 1.1) + .025, .4035, cex = 1.2, "99%", col = "black")
# Range 1
# text(qbeta(.75, 1.1, 1.1) - .025, .3375, cex = 1.2,
#   sprintf("%.2f", qbeta(.75, 1.1, 1.1) * 2 - 1))
# text(qbeta(.75, 2, 2) - .025, .2375, cex = 1.2,
#   sprintf("%.2f", qbeta(.75, 2, 2) * 2 - 1))
# text(qbeta(.75, 11, 11) - .025, .1375, cex = 1.2,
#   sprintf("%.2f", qbeta(.75, 11, 11) * 2 - 1))
# text(qbeta(.75, 101, 101) - .025, .0375, cex = 1.2,
#   sprintf("%.2f", qbeta(.75, 101, 101) * 2 - 1))
segments(qbeta(c(.25, .75), 1.1, 1.1), 0.29, y1 = .385, lwd = 2)
segments(qbeta(c(.25, .75), 2, 2), 0.19, y1 = .285, lwd = 2)
segments(qbeta(c(.25, .75), 11, 11), 0.09, y1 = .185, lwd = 2)
segments(qbeta(c(.25, .75), 101, 101), -.01, y1 = .085, lwd = 2)
# Range 2
# text(qbeta(.995, 1.1, 1.1) - .025, .3375, cex = 1.2,
#   sprintf("%.2f", qbeta(.995, 1.1, 1.1) * 2 - 1), col = "white")
# text(qbeta(.995, 2, 2) - .025, .2375, cex = 1.2,
#   sprintf("%.2f", qbeta(.995, 2, 2) * 2 - 1), col = "white")
# text(qbeta(.995, 11, 11) - .025, .1375, cex = 1.2,
#   sprintf("%.2f", qbeta(.995, 11, 11) * 2 - 1), col = "white")
# text(qbeta(.995, 101, 101) - .025, .0375, cex = 1.2,
#   sprintf("%.2f", qbeta(.995, 101, 101) * 2 - 1), col = "white")
segments(qbeta(c(.005, .995), 1.1, 1.1), 0.29, y1 = .385, lwd = 2)
segments(qbeta(c(.005, .995), 2, 2), 0.19, y1 = .285, lwd = 2)
segments(qbeta(c(.005, .995), 11, 11), 0.09, y1 = .185, lwd = 2)
segments(qbeta(c(.005, .995), 101, 101), -.01, y1 = .085, lwd = 2)
abline(v = c(0.5, .95), col = "gray", lwd = 2)
dev.off()


# Merged ---

cairo_pdf("paper/betas_lr.pdf", height = 3.8, width = 7.2,
  family = "Noto Sans", bg = "transparent", pointsize = 8)
op <- par(mar = c(0, 3.5, 2, 3.5), mfrow = c(1, 2), fig = c(0, .5815, 0.175, 1))
plot.new()
plot.window(xlim = c(-1, 0), ylim = c(0, dbeta(.5, 101, 101)))
axis(3, at = c(-1, 0, .9, 1), labels = TRUE)
axis(2, las = 1, at = q5, labels = q5l)
segments(-2, q5[4], 0, q5[4], col = "gray", lwd = 2, lty = 1)
for(i in seq_along(taus))
  lines(dbeta(x_vals_a, 1 + taus[i], 1 + taus[i]), x = x_vals_a * 2 - 1,
    col = cols[i], lty = ltys[i], lwd = 2)
text(-.2, 9, labels = "τ = 100", cex = 1.5, col = cols[4])
text(-.3, 2.8, labels = "τ = 10", cex = 1.5, col = cols[3])
text(-.45, 1.8, labels = "τ = 1", cex = 1.5, col = cols[2])
text(-.85, 1.4, labels = "τ = 0.1", cex = 1.5, col = cols[1])
abline(v = 0, col = "gray", lwd = 2)
par(fig = c(0.4185, 1, 0.175, 1), new = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(dbeta(.975, 11, 11, log = TRUE),
  dbeta(.5, 101, 101, log = TRUE)))
axis(3, at = c(-1, .9, 1), labels = TRUE)
axis(4, las = 1, at = q9, labels = q9l)
segments(0.9, q9[2], 2, q9[2], col = "gray", lwd = 2, lty = 1)
for(i in seq_along(taus))
  lines(pmax(dbeta(x_vals_b, 1 + taus[i], 1 + taus[i], log = TRUE),
    dbeta(.98, 11, 11, log = TRUE)),
    x = x_vals_b * 2 - 1, col = cols[i], lty = ltys[i], lwd = 2)
abline(v = c(0, 0.9), col = "gray", lwd = 2)
par(fig = c(0, 1, 0, 0.175), mar = c(.5, 2.5895, .5, 2.5895), new = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(-.025, .435))
axis(1, at = c(0, 1), labels = FALSE)
# rect(0, 0, 1, .075, col = "white", border = "white")
# rect(0, .1, 1, .175, col = "white", border = "white")
# rect(0, .2, 1, .275, col = "white", border = "white")
# rect(0, .3, 1, .375, col = "white", border = "white")
col_ramp <- colorRampPalette(c(cols[4], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 101, 101)
  rect(pos[1], 0, pos[2], .075, col = col_ramp[j], border = col_ramp[j])
}
col_ramp <- colorRampPalette(c(cols[3], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 11, 11)
  rect(pos[1], .1, pos[2], .175, col = col_ramp[j], border = col_ramp[j])
}
col_ramp <- colorRampPalette(c(cols[2], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 2, 2)
  rect(pos[1], .2, pos[2], .275, col = col_ramp[j], border = col_ramp[j])
}
col_ramp <- colorRampPalette(c(cols[1], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 1.1, 1.1)
  rect(pos[1], .3, pos[2], .375, col = col_ramp[j], border = col_ramp[j])
}
# segments(1, y0 = 0, y1 = .4)
# segments(0, y0 = 0, y1 = .4)
text(qbeta(.25, 1.1, 1.1) + .025, .39, cex = 1.2, "50%", col = "#333333")
text(qbeta(.1, 1.1, 1.1) + .025, .39, cex = 1.2, "80%", col = "#333333")
text(qbeta(.005, 1.1, 1.1) + .025, .39, cex = 1.2, "99%", col = "#333333")
# Range 1
# segments(qbeta(c(.25, .75), 1.1, 1.1), 0.29, y1 = .385, lwd = 2, col = "#333333")
# segments(qbeta(c(.25, .75), 2, 2), 0.19, y1 = .285, lwd = 2, col = "#333333")
# segments(qbeta(c(.25, .75), 11, 11), 0.09, y1 = .185, lwd = 2, col = "#333333")
# segments(qbeta(c(.25, .75), 101, 101), -.01, y1 = .085, lwd = 2, col = "#333333")
# Range 2
# segments(qbeta(c(.1, .9), 1.1, 1.1), 0.29, y1 = .385, lwd = 2, col = "#333333")
# segments(qbeta(c(.1, .9), 2, 2), 0.19, y1 = .285, lwd = 2, col = "#333333")
# segments(qbeta(c(.1, .9), 11, 11), 0.09, y1 = .185, lwd = 2, col = "#333333")
# segments(qbeta(c(.1, .9), 101, 101), -.01, y1 = .085, lwd = 2, col = "#333333")
# Range 3
segments(qbeta(c(.005, .995), 1.1, 1.1), 0.29, y1 = .385, lwd = 2, col = "#333333")
segments(qbeta(c(.005, .995), 2, 2), 0.19, y1 = .285, lwd = 2, col = "#333333")
segments(qbeta(c(.005, .995), 11, 11), 0.09, y1 = .185, lwd = 2, col = "#333333")
segments(qbeta(c(.005, .995), 101, 101), -.01, y1 = .085, lwd = 2, col = "#333333")
abline(v = c(0.5, .95), col = "gray", lwd = 2)
dev.off()

taus <- c(1, 10, 100)
cols <- viridisLite::viridis(3, option = "D", end = 0.8)
ltys <- c(3, 5, 1)

x_vals_a <- seq(0, .5, length.out = 1e4)
x_vals_b <- seq(.5, 1, length.out = 1e4)

q5 <- c(dbeta(.5, 2, 2), dbeta(.5, 11, 11), dbeta(.5, 101, 101))
q5l <- c(sprintf("%.2f", q5[1:2]), sprintf("%.1f", q5[3]))
q9 <- c(dbeta(.95, 2, 2, log = TRUE),
  dbeta(.95, 11, 11, log = TRUE), dbeta(.95, 101, 101, log = TRUE))
q9l <- c(sprintf("%.2f", q9[1]), sprintf("%.1f", q9[2]), sprintf("%.0f", q9[3]))
q9l[3] <- paste0("↓ ", q9l[3])

qs <- list(c(.25, .75), c(.1, .9), c(.005, .995))

cairo_pdf("paper/betas_lr.pdf", height = 3.8, width = 7.2,
  family = "Noto Sans", bg = "transparent", pointsize = 8)
op <- par(mar = c(0, 3.5, 2, 3.5), mfrow = c(1, 2), fig = c(0, .5815, 0.175, 1))
plot.new()
plot.window(xlim = c(-1, 0), ylim = c(0, dbeta(.5, 101, 101)))
axis(3, at = c(-1, 0, .9, 1), labels = TRUE, cex = 1.5)
axis(2, las = 1, at = q5, labels = q5l, cex = 1.5)
segments(-2, q5[3], 0, q5[3], col = "gray", lwd = 2, lty = 1)
for(i in seq_along(rev(taus)))
  lines(dbeta(x_vals_a, 1 + taus[i], 1 + taus[i]), x = x_vals_a * 2 - 1,
    col = cols[i], lty = ltys[i], lwd = 2)
text(-.2, 9, labels = "τ = 100", cex = 1.8, col = cols[3])
text(-.3, 2.8, labels = "τ = 10", cex = 1.8, col = cols[2])
text(-.6, 1.75, labels = "τ = 1", cex = 1.8, col = cols[1])
# text(-.85, 1.4, labels = "τ = 0.1", cex = 1.5, col = cols[1])
abline(v = 0, col = "gray", lwd = 2)
par(fig = c(0.4185, 1, 0.175, 1), new = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(dbeta(.975, 11, 11, log = TRUE),
  dbeta(.5, 101, 101, log = TRUE)))
axis(3, at = c(-1, .9, 1), labels = TRUE, cex = 1.5)
axis(4, las = 1, at = c(q9[1:2], dbeta(.975, 11, 11, log = TRUE)),
  labels = q9l, cex = 1.5)
segments(0.9, q9[1], 2, q9[1], col = "gray", lwd = 2, lty = 1)
for(i in seq_along(taus))
  lines(pmax(dbeta(x_vals_b, 1 + taus[i], 1 + taus[i], log = TRUE),
    dbeta(.98, 11, 11, log = TRUE)),
    x = x_vals_b * 2 - 1, col = cols[i], lty = ltys[i], lwd = 2)
abline(v = c(0, 0.9), col = "gray", lwd = 2)
par(fig = c(0, 1, 0, 0.175), mar = c(.5, 2.5895, 0, 2.5895), new = TRUE)
plot.new()
axis(1, at = c(0, .25, .5, .75, 1), labels = FALSE)
plot.window(xlim = c(0, 1), ylim = c(-.025, .335))
# axis(1, at = c(0, 1), labels = FALSE)
# rect(0, 0, 1, .075, col = "white", border = "white")
# rect(0, .1, 1, .175, col = "white", border = "white")
# rect(0, .2, 1, .275, col = "white", border = "white")
# rect(0, .3, 1, .375, col = "white", border = "white")
col_ramp <- colorRampPalette(c(cols[3], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 101, 101)
  rect(pos[1], 0, pos[2], .075, col = col_ramp[j], border = col_ramp[j])
}
col_ramp <- colorRampPalette(c(cols[2], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 11, 11)
  rect(pos[1], .1, pos[2], .175, col = col_ramp[j], border = col_ramp[j])
}
col_ramp <- colorRampPalette(c(cols[1], "white"))(6)[5:2]
for(j in seq(length(qs), 1)) {
  pos <- qbeta(qs[[j]], 2, 2)
  rect(pos[1], .2, pos[2], .275, col = col_ramp[j], border = col_ramp[j])
}
segments(1, y0 = 0, y1 = .275)
segments(0, y0 = 0, y1 = .275)
text(qbeta(.25, 2, 2) + .03, .3, cex = 1.5, "50%", col = "#333333")
text(qbeta(.1, 2, 2) + .03, .3, cex = 1.5, "80%", col = "#333333")
text(qbeta(.005, 2, 2) + .03, .3, cex = 1.5, "99%", col = "#333333")

segments(qbeta(c(.005, .995), 2, 2), 0.19, y1 = .285, lwd = 2, col = "#333333")
segments(qbeta(c(.005, .995), 11, 11), 0.09, y1 = .185, lwd = 2, col = "#333333")
segments(qbeta(c(.005, .995), 101, 101), -.01, y1 = .085, lwd = 2, col = "#333333")
abline(v = c(0.5, .95), col = "gray", lwd = 2)
dev.off()
