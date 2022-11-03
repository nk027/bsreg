
N <- 100

# Simulate connectivity
xy <- cbind(runif(N), runif(N))
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
Psi <- function(x) {
  W <- 1 / dist ^ x
  W / max(eigen(W)$values)
}

i_lambda <- seq(-1 + 1e-12, 1 - 1e-12, len = 50)
i_delta <- seq(1e-12, 10, len = 20)
ldets <- cbind("ldet" = NA_real_, as.matrix(expand.grid(i_lambda, i_delta)))

for(i in seq_along(i_delta)) {
  ev <- eigen(Psi(i_delta[i]), symmetric = TRUE, only.values = TRUE)$values
  ldets[seq((i - 1) * length(i_lambda) + 1, i * length(i_lambda)), 3] <- i_delta[i]
  ldets[seq((i - 1) * length(i_lambda) + 1, i * length(i_lambda)), 2] <- i_lambda
  ldets[seq((i - 1) * length(i_lambda) + 1, i * length(i_lambda)), 1] <- vapply(i_lambda,  \(lambda) {
    Re(sum(log(1 - lambda * ev))) * 1L
  }, numeric(1L))
}

i_delta_oos <- seq(10.5, 15, len = 10)
ldets_oos <- cbind("ldet" = NA_real_, as.matrix(expand.grid(i_lambda, i_delta_oos)))

for(i in seq_along(i_delta_oos)) {
  ev <- eigen(Psi(i_delta_oos[i]), symmetric = TRUE, only.values = TRUE)$values
  ldets_oos[seq((i - 1) * length(i_lambda) + 1, i * length(i_lambda)), 3] <- i_delta_oos[i]
  ldets_oos[seq((i - 1) * length(i_lambda) + 1, i * length(i_lambda)), 2] <- i_lambda
  ldets_oos[seq((i - 1) * length(i_lambda) + 1, i * length(i_lambda)), 1] <- vapply(i_lambda,  \(lambda) {
    Re(sum(log(1 - lambda * ev))) * 1L
  }, numeric(1L))
}

gp <- GauPro::GauPro(ldets[, 2:3], ldets[, 1], parallel = FALSE)
# The log-determinent is predicted using the Gaussian process
get_ldet <- function(lambda, delta, ...) {
  gp$predict(cbind(lambda, delta))[[1]]
}

plot.new()
plot.window(xlim = c(-1, 1), ylim = c(0, 15))
axis(1); axis(2)
points(ldets[, 2], ldets[, 3])

library("ggplot2")

# General
approx <- cbind("ldet" = NA_real_,
  as.matrix(expand.grid(seq(-1, 1, len = 250), seq(1e-12, 15, len = 100))))
# Compare out-of-sample
approx <- cbind("ldet" = NA_real_,
  as.matrix(expand.grid(i_lambda, c(i_delta, i_delta_oos))))

approx[, 1] <- apply(approx, 1, \(x) get_ldet(x[2], x[3]))

colnames(ldets) <- colnames(approx) <- colnames(ldets_oos) <- c("ldet", "lambda", "delta")
diff <- merge(ldets_oos, approx, by = c("delta", "lambda"))

ldets |> as.data.frame() |>
  ggplot(aes(x = lambda, y = delta)) +
  geom_tile(data = as.data.frame(approx), aes(fill = ldet)) +
  geom_point(aes(fill = ldet), colour = "black", pch = 22, size = 4) +
  geom_point(data = as.data.frame(ldets_oos),
    aes(fill = ldet), colour = "black", pch = 23, size = 4) +
  geom_point(data = as.data.frame(diff),
    aes(col = abs((ldet.x - ldet.y))), pch = 16, size = 5) +
  coord_cartesian(ylim = c(0, 15)) +
  scale_fill_viridis_c(option = "magma", name = "Log-determinant") +
  scale_color_viridis_c(name = "|Error|") +
  theme_minimal()


# Consider approximation using different grids ---

x <- pracma::gaussLegendre(100, -1 + 1e-12, 1 - 1e-12)
plot(x$x, x$w)
x <- pracma::gaussLegendre(100, 0, 10)
plot(x$x, x$w)

x <- laguerre.quadrature.rules(20)
plot(x[[20]][, 1], x[[20]][, 2])
lines(x[[20]][, 1], x[[20]][, 2], col = "darkgray")
points(x[[10]][, 1], x[[10]][, 2], pch = 2, col = 2)
points(x[[5]][, 1], x[[5]][, 2], pch = 3, col = 3)
