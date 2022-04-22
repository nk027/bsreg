
# Generate outputs ---
library("dplyr")
library("spatialreg")
library("ggplot2")
library("ggdist") # Dotplots
library("ggthemes")

load("paper/cigarettes_dist.Rda")

# Prepare object ---
d_te <- rbind(
  as_tibble(out_bslx$draws) %>% transmute(model = "SLX",
    price = beta2, price_i = beta78,
    income = beta3, income_i = beta79, delta = NA_real_),
  as_tibble(out_slxd2$draws) %>% transmute(model = "SLX(2)",
    price = beta2, price_i = beta78,
    income = beta3, income_i = beta79, delta = 2),
  as_tibble(out_slxd3$draws) %>% transmute(model = "SLX(3)",
    price = beta2, price_i = beta78,
    income = beta3, income_i = beta79, delta = 3),
  as_tibble(out_slxd4$draws) %>% transmute(model = "SLX(4)",
    price = beta2, price_i = beta78,
    income = beta3, income_i = beta79, delta = 4),
  as_tibble(out_slxdx$draws) %>% transmute(model = "SLX(delta)",
    price = beta2, price_i = beta78,
    income = beta3, income_i = beta79, delta = delta_SLX)
)

# Build a table of coefficients ---
tbl <- d_te %>% filter(model != "SLX") %>%
  tidyr::pivot_longer(cols = 2:6) %>%
  group_by(model, name) %>%
  summarise_all(list(mean = mean, sd = sd, # Report mean, sd, and quantiles
    qu01 = \(x) quantile(x, .01), qu99 = \(x) quantile(x, .99))) %>%
  tidyr::pivot_longer(cols = 3:6, names_to = "measure") %>%
  mutate(value = round(value, 3)) %>%
  mutate(value = ifelse(measure == "sd", # Brackets around sd
    gsub("(.*)", "(\\1)", value), value)) %>%
  mutate(value = ifelse(measure == "qu01", # Square brackets around quantiles
    gsub("(.*)", "[\\1, ", value), value)) %>%
  mutate(value = ifelse(measure == "qu99",
    gsub("(.*)", " \\1]", value), value)) %>%
  tidyr::pivot_wider(names_from = model)

# Fix order
tbl <- tbl[c(13:16, 5:8, 17:20, 9:12, 1:4), ]
# Move the quantiles to one row with the credible interval
for(r in c(3, 7, 11, 15, 19)) {
  tbl[r, 2:6] <- as.list(c("ci", paste0(tbl[r, 3:6], tbl[r + 1, 3:6])))
}
tbl <- tbl[-c(4, 8, 12, 16, 20), ]

# Export to LaTeX
tbl %>% as.data.frame() %>% memisc:::toLatex.data.frame()

# Total effect plots ---

# The size of indirect effects depends on delta via Psi
N <- nrow(X)
delta <- seq(2, 4, 0.01)
avg_weight <- sapply(delta, \(x) sum(Psi(x)) / N)
s <- splinefun(delta, avg_weight)

# Scaling factor for the table
round(sapply(c(2, 3, 4, 3.01, 2.54, 3.74), s), 2)

d_te <- d_te %>% mutate(
  price_t = price + price_i * s(delta),
  income_t = income + income_i * s(delta)) %>%
  mutate( # Row-stochastic is simpler
    price_t = ifelse(model == "SLX", price + price_i, price_t),
    income_t = ifelse(model == "SLX", income + income_i, income_t)
  )

# Check summaries to manually limit the y axis
d_te %>% tidyr::pivot_longer(cols = 7:8) %>%
  group_by(model, name) %>%
  summarise(min = min(value), max = max(value),
    q9 = quantile(value, .99), q1 = quantile(value, .01))

# Fix colour of SLX with binary contiguity
cols <- ggthemes::colorblind_pal()(7)[c(5, 6, 7, 4, 3, 2)]

# Plots
p1 <- d_te %>% transmute(model, price_t) %>%
  tidyr::pivot_longer(cols = 2) %>%
  # Use unicode delta
  mutate(model = ifelse(model == "SLX(delta)", "SLX(δ)", model)) %>%
  ggplot(aes(x = model, y = value, fill = model, col = model)) +
  geom_hline(yintercept = 0, lwd = 1.5, col = "#a4a4a4") +
  stat_dots(quantiles = 250, width = .75, justification = -0.2) +
  geom_boxplot(col = "#444444", alpha = 1, # No points for outliers
    width = .2, size = .8, outlier.color = NA) +
  coord_cartesian(ylim = c(-1.5, -0.6)) + # Manually set limits
  ggtitle("Price effect by connectivity") + ylab("total effect") + xlab("") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(color = "#333333", size = 18, face = "bold"),
    axis.title.x = element_text(color = "#333333", size = 14, face = "bold"),
    axis.title.y = element_text(color = "#333333", size = 14, face = "bold"),
    text = element_text(family = "Helvetica"),
    legend.position = "none"
  ) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols)

p2 <- d_te %>% transmute(model, income_t) %>%
  tidyr::pivot_longer(cols = 2) %>%
  # Use unicode delta
  mutate(model = ifelse(model == "SLX(delta)", "SLX(δ)", model)) %>%
  ggplot(aes(x = model, y = value, fill = model, col = model)) +
  geom_hline(yintercept = 0, lwd = 1.5, col = "#a4a4a4") +
  stat_dots(quantiles = 250, width = .75, justification = -0.2) +
  geom_boxplot(col = "#444444", alpha = 1, # No points for outliers
    width = .2, size = .8, outlier.color = NA) +
  coord_cartesian(ylim = c(0.1, 0.8)) +
  ggtitle("Income effect by connectivity") + ylab("total effect") + xlab("") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(color = "#333333", size = 18, face = "bold"),
    axis.title.x = element_text(color = "#333333", size = 14, face = "bold"),
    axis.title.y = element_text(color = "#333333", size = 14, face = "bold"),
    text = element_text(family = "Helvetica"),
    legend.position = "none",
  ) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols)

# Merge and save the plots
gridExtra::grid.arrange(p1, p2, ncol = 2)
ggsave(file = "cigar-dist_te.eps", # Cairo for unicode delta
  plot = gridExtra::arrangeGrob(p1, p2, nrow = 1, ncol = 2),
  device = cairo_pdf, width = 9, height = 4)

# Diagnostic plot for delta ---

# Visualise connectivity ---
d_delta <- as_tibble(list(
  "iteration" = seq_len(nrow(out_slxdx$draws)),
  "delta" = out_slxdx$draws[, "delta_SLX"]
))

# Trace plot
p1 <- ggplot(d_delta) +
  geom_line(aes(y = delta, x = iteration)) +
  ggtitle("Connectivity parameter") + ylab("δ") + xlab("iteration") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(color = "#333333", size = 18, face = "bold"),
    axis.title.x = element_text(color = "#333333", size = 14, face = "bold"),
    axis.title.y = element_text(color = "#333333", size = 14, face = "bold"),
    text = element_text(family = "Helvetica"),
    legend.position = "none"
  )
# Density plot
p2 <- ggplot(d_delta) +
  stat_halfeye(aes(y = delta)) +
  ggtitle("") + ylab("δ") + xlab("density") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(color = "#333333", size = 18, face = "bold"),
    axis.title.x = element_text(color = "#333333", size = 14, face = "bold"),
    axis.title.y = element_text(color = "#333333", size = 14, face = "bold"),
    text = element_text(family = "Helvetica"),
    legend.position = "none"
  )

# Merge and save the plots
gridExtra::grid.arrange(p1, p2, ncol = 2)
ggsave(file = "cigar-dist_delta.eps", # Cairo for unicode delta
  plot = gridExtra::arrangeGrob(p1, p2, nrow = 1, ncol = 2),
    device = cairo_pdf, width = 9, height = 4)
