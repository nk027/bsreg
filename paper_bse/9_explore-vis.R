
# Lambda priors -----

library("ggplot2")
library("ggdist")

# Visualise ranges ---
p <- data.frame(group = c("100", "10", "1", "0.1"), lambda = c(101, 11, 2, 1.1)) |>
  ggplot(aes(x = group)) +
  scale_color_ramp_discrete(range = c(1, .1), na.translate = FALSE) +
  scale_color_viridis_d(option = "D", end = .7) +
  coord_flip() +
  scale_y_continuous(breaks = c(0, .5, 1), labels = c("-1", "0", "1")) +
  labs(title = "", subtitle = "", fill = "τ", fill_ramp = "Interval",
    color = "τ", color_ramp = "Interval") +
  theme_minimal()

p + stat_histinterval(aes(ydist = distributional::dist_beta(lambda, lambda), fill = stat(pdf))) +
  scale_fill_viridis_c(option = "D", end = .7)
p + stat_slab(aes(ydist = distributional::dist_beta(lambda, lambda), fill = group,
    fill_ramp = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .9, .95, 99), labels = scales::percent_format(accuracy = 1)))),
    position = "dodgejust", height = 2, show_interval = FALSE) +
  scale_fill_viridis_d(option = "D", end = .7)
p + stat_interval(aes(ydist = distributional::dist_beta(lambda, lambda),
    color_ramp = stat(level), color = group),
    .width = c(0.5, 0.9, 0.95, .99), size = 10, position = "dodge")


# Test shrinkage of \lambda_a versus \lambda_b
n <- 100000
x_2 <- x_4 <- x_a <- numeric(n + 1)
y_2 <- y_4 <- y_a <- numeric(n)
x_2[1] <- x_4[1] <- x_a[1] <- .5
for(i in seq(n)) {
  # Absolute
  y_a[i] <- rbeta(1L, 1 + 10 * abs(2 * x_a[i] - 1), 1 + 10 * abs(2 * x_a[i] - 1))
  x_a[i + 1] <- rbeta(1L, 1 + 10 * abs(2 * y_a[i] - 1), 1 + 10 * abs(2 * y_a[i] - 1))
  # Squared
  y_2[i] <- rbeta(1L, 1 + 10 * (2 * x_2[i] - 1)^2, 1 + 10 * (2 * x_2[i] - 1)^2)
  x_2[i + 1] <- rbeta(1L, 1 + 10 * (2 * y_2[i] - 1)^2, 1 + 10 * (2 * y_2[i] - 1)^2)
  # To the fourth
  y_4[i] <- rbeta(1L, 1 + 10 * (2 * x_4[i] - 1)^4, 1 + 10 * (2 * x_4[i] - 1)^4)
  x_4[i + 1] <- rbeta(1L, 1 + 10 * (2 * y_4[i] - 1)^4, 1 + 10 * (2 * y_4[i] - 1)^4)
}

op <- par(mfrow = c(2, 2))
plot(x_a[-1], y_a, main = "|x|")
plot(x_2[-1], y_2, main = "x^2")
plot(x_4[-1], y_4, main = "x^4")
plot(density(x_a))
lines(density(x_2), lty = 2)
lines(density(x_4), lty = 3)
par(op)

# Add contours
op <- par(mfrow = c(2, 2))
plot(x_a[-1], y_a, main = "|x|")
contour(MASS::kde2d(x_a[-1], y_a), lwd = 2, add = TRUE, col = viridis::viridis(12))
plot(x_2[-1], y_2, main = "x^2")
contour(MASS::kde2d(x_2[-1], y_2), lwd = 2, add = TRUE, col = viridis::viridis(12))
plot(x_4[-1], y_4, main = "x^4")
contour(MASS::kde2d(x_4[-1], y_4), lwd = 2, add = TRUE, col = viridis::viridis(12))
plot(density(x_a))
lines(density(x_2), lty = 2)
lines(density(x_4), lty = 3)
par(op)

s <- seq(.01, .99, len = 101)
op <- par(mfrow = c(2, 1))
plot(s, 1 + 10 * abs(2 * s - 1), ylim = c(1, 11), main = "Shrinkage level")
lines(s, 1 + 10 * (2 * s - 1)^2)
points(s, 1 + 10 * (2 * s - 1)^4, pch = 2)
lines(s, 1 + 10 * (2 * s - 1)^6, lty = 2)
points(s, 1 + 10 * (2 * s - 1)^8, pch = 3)
lines(s, 1 + 10 * (2 * s - 1)^16, lty = 3)
points(s, 1 + 10 * (2 * s - 1)^32, pch = 4)

plot(s, log(1 + 10 * abs(2 * s - 1)), ylim = log(c(1, 11)), main = "Log")
lines(s, log(1 + 10 * (2 * s - 1)^2))
points(s, log(1 + 10 * (2 * s - 1)^4), pch = 2)
lines(s, log(1 + 10 * (2 * s - 1)^6), lty = 2)
points(s, log(1 + 10 * (2 * s - 1)^8), pch = 3)
lines(s, log(1 + 10 * (2 * s - 1)^16), lty = 3)
points(s, log(1 + 10 * (2 * s - 1)^32), pch = 4)
par(op)


# knn priors -----

x <- seq(0, 100)

# Beta Binomial and Beta Negative-Binomial
op <- par(mfrow = c(2, 1), mar = c(2, 2, 2, .5))
plot(x, extraDistr::dbbinom(x, size = 100L, alpha = 2, beta = 5),
  pch = 15, col = 1, ylim = c(0, .05), main = "BB")
points(x, extraDistr::dbbinom(x, size = 100L, alpha = 2, beta = 10),
  pch = 15, col = 3)
points(x, extraDistr::dbbinom(x, size = 100L, alpha = 5, beta = 10),
  pch = 18, col = 5)
points(x, extraDistr::dbbinom(x, size = 100L, alpha = 5, beta = 20),
  pch = 18, col = 6)

plot(x, extraDistr::dbnbinom(x, size = 100L, alpha = 5, beta = 2),
  pch = 15, col = 1, ylim = c(0, .05), main = "BNB")
points(x, extraDistr::dbnbinom(x, size = 100L, alpha = 10, beta = 2),
  pch = 15, col = 3)
points(x, extraDistr::dbnbinom(x, size = 100L, alpha = 10, beta = 5),
  pch = 18, col = 5)
points(x, extraDistr::dbnbinom(x, size = 100L, alpha = 20, beta = 5),
  pch = 18, col = 6)
par(op)

