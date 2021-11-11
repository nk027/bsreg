
devtools::load_all()

df <- readxl::read_excel("~/repos/bse/data/cigarette+2var.xls")
contig_data <- readxl::read_excel("~/repos/bse/data/Spat-Sym-US.xls", col_names = FALSE)
xy_data <- readxl::read_excel("~/repos/bse/data/cigar_states.xls", col_names = TRUE)

# Utah has a wrong longitude
xy_data_fixed <- xy_data
xy_data_fixed[40, 3] <- 112

x <- as.matrix(df[, c("logc", "logp", "logy")])
y <- x[, 1]
X <- cbind(1, x[, -1])
# Add TFE
n_time <- length(unique(df$year))
X <- cbind(X, kronecker(diag(n_time), matrix(1, nrow(x) / n_time))[, -1])
# Add IFE
X <- cbind(X, kronecker(matrix(1, n_time), diag(nrow(x) / n_time))[, -1])

# Use Frisch-Waugh-Lovell theorem to get rid of fixed effects
# variables <- 1:3
# Q_fwl <- qr.Q(qr(X[, -variables, drop = FALSE]))
# y <- y - Q_fwl %*% crossprod(Q_fwl, y)
# X <- X[, variables, drop = FALSE] - Q_fwl %*% # We don't care about the intercept
#   crossprod(Q_fwl, X[, variables, drop = FALSE])

# Connectivity
W <- as.matrix(contig_data)
W <- kronecker(diag(n_time), W / rowSums(W))
# W <- kronecker(diag(n_time), W / max(Re(eigen(W)$values)))
dist <- as.matrix(dist(xy_data))
dist <- as.matrix(dist(xy_data_fixed))
diag(dist) <- Inf
Psi <- function(x) {
  W <- 1 / dist ^ x
  W <- W / max(Re(eigen(W[1:46, 1:46])$values))
  kronecker(diag(n_time), W)
}
# W <- Psi(3)

# Reproduce -----

n_save <- 1000L
n_burn <- 500L

# Row-stochastic binary contiguity ---

X_LX <- cbind(X, W %*% X[, 2:3]) # Easier for lm and spatialreg

(out_blm <- blm(y ~ X - 1, n_save = n_save, n_burn = n_burn))
(out_lm <- lm(y ~ X - 1))

(out_bslx <- bslx(y ~ X - 1, W = W, X_SLX = X[, 2:3],
  n_save = n_save, n_burn = n_burn))
(out_slx <- lm(y ~ X_LX - 1))

(out_bsar <- bsar(y ~ X - 1, W = W, n_save = n_save, n_burn = n_burn,
  ldet_SAR = list(reps = n_time)))
(out_sar <- spatialreg::lagsarlm(y ~ X - 1, listw = spdep::mat2listw(W)))

(out_bsem <- bsem(y ~ X, W = W, n_save = n_save, n_burn = n_burn,
  ldet_SEM = list(reps = n_time)))
(out_sem <- spatialreg::errorsarlm(y ~ X - 1, listw = spdep::mat2listw(W)))

(out_bsdm <- bsdm(y ~ X, W = W, X_SLX = X[, 2:3],
  n_save = n_save, n_burn = n_burn, ldet_SAR = list(reps = n_time)))
(out_sdm <- spatialreg::lagsarlm(y ~ X_LX - 1, listw = spdep::mat2listw(W)))

(out_bsdem <- bsdem(y ~ X, W = W, X_SLX = X[, 2:3],
  n_save = n_save, n_burn = n_burn, ldet_SEM = list(reps = n_time)))
(out_sdem <- spatialreg::errorsarlm(y ~ X_LX - 1, listw = spdep::mat2listw(W)))

# Plot total effects
library("dplyr")
library("ggplot2")
library("ggdist")
d1 <- rbind(
  as_tibble(out_blm$draws) %>% transmute(model = "LM",
    price = beta2, income = beta3),
  as_tibble(out_bslx$draws) %>% transmute(model = "SLX",
    price = beta2 + beta78, income = beta3 + beta79),
  as_tibble(out_bsar$draws) %>% transmute(model = "SAR",
    price = beta2 / (1 - lambda_SAR), income = beta3 / (1 - lambda_SAR)),
  as_tibble(out_bsem$draws) %>% transmute(model = "SEM",
    price = beta2, income = beta3),
  as_tibble(out_bsdm$draws) %>% transmute(model = "SDM",
    price = (beta2 + beta78) / (1 - lambda_SAR),
    income = (beta3 + beta79) / (1 - lambda_SAR)),
  as_tibble(out_bsdem$draws) %>% transmute(model = "SDEM",
    price = beta2 + beta78, income = beta3 + beta79)
) %>% tidyr::pivot_longer(cols = 2:3)
d2 <- rbind(
  tibble(model = "LM", price = coef(out_lm)[2], income = coef(out_lm)[3]),
  tibble(model = "SLX",
    price = sum(coef(out_slx)[c(2, 78)]), income = sum(coef(out_slx)[c(3, 79)])),
  tibble(model = "SAR", price = coef(out_sar)[3] / (1 - coef(out_sar)[1]),
    income = coef(out_sar)[4] / (1 - coef(out_sar)[1])),
  tibble(model = "SEM", price = coef(out_sem)[3], income = coef(out_sem)[4]),
  tibble(model = "SDM",
    price = sum(coef(out_sdm)[c(3, 79)]) / (1 - coef(out_sdm)[1]),
    income = sum(coef(out_sdm)[c(4, 80)]) / (1 - coef(out_sdm)[1])),
  tibble(model = "SDEM",
    price = sum(coef(out_sdem)[c(3, 79)]), income = sum(coef(out_sdem)[c(4, 80)]))
) %>% tidyr::pivot_longer(cols = 2:3)

d1 <- d1 %>% mutate(model = factor(model,
  levels = c("LM", "SEM", "SDEM", "SLX", "SDM", "SAR")))

d1 %>% filter(model != "SDM") %>% group_by(name) %>%
  summarise(min = min(value), max = max(value))

p1 <- d1 %>% filter(name == "price") %>% ggplot() +
  stat_dots(aes(x = model, y = value, fill = model, col = model), quantiles = 250,
    width = .75, justification = -0.2) +
  ggplot2::geom_boxplot(aes(x = model, y = value, fill = model), col = "#444444",
    alpha = 0.75, width = .2, size = .8, outlier.color = NA) +
  geom_point(data = d2,
    aes(x = model, y = value), shape = 4, stroke = 1.5, size = 3) +
  ggthemes::theme_stata(base_size = 14) +
  theme(
    plot.title = element_text(color = "#333333", size = 18, face = "bold"),
    axis.title.x = element_text(color = "#333333", size = 14, face = "bold"),
    axis.title.y = element_text(color = "#333333", size = 14, face = "bold"),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(-2.2, 0)) +
  ggtitle("Total effect of price and income by model") +
  scale_colour_manual(values = ggthemes::colorblind_pal()(7)[-1]) +
  scale_fill_manual(values = ggthemes::colorblind_pal()(7)[-1])

p2 <- d1 %>% filter(name == "income") %>% ggplot() +
  stat_dots(aes(x = model, y = value, fill = model, col = model), quantiles = 250,
    width = .75, justification = -0.2) +
  ggplot2::geom_boxplot(aes(x = model, y = value, fill = model), col = "#444444",
    alpha = 0.75, width = .2, size = .8, outlier.color = NA) +
  geom_point(data = d2,
    aes(x = model, y = value), shape = 4, stroke = 1.5, size = 3) +
  ggthemes::theme_stata(base_size = 14) +
  theme(
    plot.title = element_text(color = "#333333", size = 18, face = "bold"),
    axis.title.x = element_text(color = "#333333", size = 14, face = "bold"),
    axis.title.y = element_text(color = "#333333", size = 14, face = "bold"),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(-.2, 1.8)) +
  ggtitle("Total effect of price and income by model") +
  scale_colour_manual(values = ggthemes::colorblind_pal()(7)[-1]) +
  scale_fill_manual(values = ggthemes::colorblind_pal()(7)[-1])

gridExtra::grid.arrange(p1, p2, nrow = 2, ncol = 1)

# Table of coefficients


# Inverse-distance decay ---

out_slxd2 <- bslx(y ~ X, W = Psi(2), X_SLX = X[, 2:3],
  n_save = n_save, n_burn = n_burn)
print(out_slxd2)

out_slxd3 <- bslx(y ~ X, W = Psi(3), X_SLX = X[, 2:3],
  n_save = n_save, n_burn = n_burn)
print(out_slxd3)

out_slxd4 <- bslx(y ~ X, W = Psi(4), X_SLX = X[, 2:3],
  n_save = n_save, n_burn = n_burn)
print(out_slxd4)

out_slxdx <- bslx(y ~ X, W = Psi, X_SLX = X[, 2:3],
  n_save = n_save, n_burn = n_burn,
  options = set_options(SLX = set_SLX(delta = 3, delta_scale = 0.05)))
print(out_slxdx)

save.image()
