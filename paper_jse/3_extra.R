
library("ggplot2")

# Lazy BIC for comparing fixed SLX(delta) variants ---
bics <- c(
  "cont" = BIC(lm(y ~ cbind(X, X_cont))),
  "dist2" = BIC(lm(y ~ cbind(X, Psi(2) %*% X_lag))),
  "dist3" = BIC(lm(y ~ cbind(X, Psi(3) %*% X_lag))),
  "dist4" = BIC(lm(y ~ cbind(X, Psi(4) %*% X_lag)))
)
bics_s <- bics - min(bics) # Standardise for numerics
exp(bics_s / -2) / sum(exp(bics_s / -2)) # Posterior probabilities

# Visualise inverse-distance connectivity
distances <- seq(1, 5, length.out = 200)
deltas <- seq(1, 5, length.out = 200)
d_conn <- expand.grid(distance = distances, delta = deltas)
idd <- function(delta, distance) distance^(-delta)
d_conn$connectivity <- NA_real_
for(i in seq_len(nrow(d_conn)))
  d_conn[i, "connectivity"] <- idd(d_conn[i, "delta"], d_conn[i, "distance"])

# Plot the connectivity as function of delta and the distance
p_idd <- ggplot(d_conn, aes(x = distance, y = delta, fill = connectivity)) +
  geom_tile(col = "transparent", alpha = 1) +
  ggtitle("Distance-decay connectivity strength") +
  ylab("δ") + xlab("distance") +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = "transparent"),
    plot.title = element_text(color = "#333333", size = 18, face = "bold"),
    axis.title.x = element_text(color = "#333333", size = 14, face = "bold"),
    axis.title.y = element_text(color = "#333333", size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "#333333"),
    legend.text = element_text(color = "#333333", size = 14, face = "bold"),
    legend.position = c(0.90, 0.85),
    text = element_text(family = "Helvetica")
  ) +
  scale_fill_viridis_c()

# Save the plot
p_idd
# png("dist-decay.png", # Base device to avoid bugged outputs
#   height = 1800, width = 1800, res = 300, bg = "transparent")
# p_idd
# dev.off()
