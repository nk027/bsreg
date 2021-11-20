
# Generate outputs ---
library("ggplot2")
library("ggdist") # Dotplots

# Compute total effects ---
d_te1 <- rbind(
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
) %>% tidyr::pivot_longer(cols = 2:3) %>%
  mutate(model = factor(model, # Set ordering
    levels = c("LM", "SEM", "SDEM", "SLX", "SDM", "SAR")))
d_te2 <- rbind(
  tibble(model = "LM",
    price = coef(out_lm)[2], income = coef(out_lm)[3]),
  tibble(model = "SLX",
    price = sum(coef(out_slx)[c(2, 78)]),
    income = sum(coef(out_slx)[c(3, 79)])),
  tibble(model = "SAR",
    price = coef(out_sar)[3] / (1 - coef(out_sar)[1]),
    income = coef(out_sar)[4] / (1 - coef(out_sar)[1])),
  tibble(model = "SEM",
    price = coef(out_sem)[3], income = coef(out_sem)[4]),
  tibble(model = "SDM",
    price = sum(coef(out_sdm)[c(3, 79)]) / (1 - coef(out_sdm)[1]),
    income = sum(coef(out_sdm)[c(4, 80)]) / (1 - coef(out_sdm)[1])),
  tibble(model = "SDEM",
    price = sum(coef(out_sdem)[c(3, 79)]),
    income = sum(coef(out_sdem)[c(4, 80)]))
) %>% tidyr::pivot_longer(cols = 2:3)

# Manually limit the y axis, ignoring the SDM's fat tails
d_te1 %>% filter(model != "SDM") %>%
  group_by(name) %>% summarise(min = min(value), max = max(value),
  q9 = quantile(value, .999), q1 = quantile(value, .001))

# Total effect plots ---
p1 <- d_te1 %>% filter(name == "price") %>%
  ggplot(aes(x = model, y = value, fill = model)) +
  geom_hline(yintercept = 0, lwd = 1.5, col = "#a4a4a4") +
  stat_dots(aes(col = model), quantiles = 250,
    width = .75, justification = -0.2) + # 250 dots, narrower & shifted right
  geom_boxplot(col = "#444444", alpha = 0.75, # No points for outliers
    width = .2, size = .8, outlier.color = NA) +
  geom_point(data = d_te2, aes(x = model, y = value), # Add freq. estimates
    col = "#444444", shape = 4, stroke = 1.5, size = 3) +
  coord_cartesian(ylim = c(-1.8, -.7)) + # Manually set limits
  ggtitle("Price effect by model") + ylab("total effect") + xlab("") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(color = "#333333", size = 18, face = "bold"),
    axis.title.x = element_text(color = "#333333", size = 14, face = "bold"),
    axis.title.y = element_text(color = "#333333", size = 14, face = "bold"),
    text = element_text(family = "Helvetica"),
    legend.position = "none"
  ) +
  scale_colour_manual(values = ggthemes::colorblind_pal()(7)[-1]) +
  scale_fill_manual(values = ggthemes::colorblind_pal()(7)[-1])

p2 <- d_te1 %>% filter(name == "income") %>%
  ggplot(aes(x = model, y = value, fill = model)) +
  geom_hline(yintercept = 0, lwd = 1.5, col = "#a4a4a4") +
  stat_dots(aes(col = model), quantiles = 250,
    width = .75, justification = -0.2) + # 250 dots, narrower & shifted right
  geom_boxplot(col = "#444444", alpha = 0.75, # No points for outliers
    width = .2, size = .8, outlier.color = NA) +
  geom_point(data = d_te2, aes(x = model, y = value), # Add frequentist estimates
    col = "#444444", shape = 4, stroke = 1.5, size = 3) +
  coord_cartesian(ylim = c(-0.1, 1.1)) + # Manually set limits
  ggtitle("Income effect by model") + ylab("total effect") + xlab("") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(color = "#333333", size = 18, face = "bold"),
    axis.title.x = element_text(color = "#333333", size = 14, face = "bold"),
    axis.title.y = element_text(color = "#333333", size = 14, face = "bold"),
    text = element_text(family = "Helvetica"),
    legend.position = "none"
  ) +
  scale_colour_manual(values = ggthemes::colorblind_pal()(7)[-1]) +
  scale_fill_manual(values = ggthemes::colorblind_pal()(7)[-1])

# Merge and save the plots
gridExtra::grid.arrange(p1, p2, nrow = 2)
# ggsave(file = "cigar-contig_te.pdf",
  # plot = gridExtra::arrangeGrob(p1, p2, nrow = 2), width = 9, height = 8)

# Investigate distribution of total effects ---
d_qq1 <- d_te1 %>% # Focus on price in the SLX, SDM, and SAR models
  filter(model %in% c("SLX", "SDM", "SAR"), name == "price")
set.seed(27) # We plot a random subset of all posterior values
d_qq2 <- d_qq1 %>% slice(c(sample(1:25000, 5000), # Assuming there's 25,000 draws
  sample(25001:50000, 5000), sample(50001:75000, 5000)))

# QQ plots for the distribution of total effects ---
pq1 <- d_qq1 %>%
  ggplot(mapping = aes(sample = value, col = model, fill = model)) +
  facet_grid(. ~ model, scale = "free") + # A model per column
  qqplotr::stat_qq_point(data = d_qq2, pch = 4, size = 0.5, col = "#333333") +
  qqplotr::stat_qq_band(alpha = 0.5, conf = 0.999, band = "pointwise") +
  qqplotr:stat_qq_line() +
  ggtitle("QQ plot of total price effect") +
  ylab("sample quantiles") + xlab("theoretical quantiles") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(color = "#333333", size = 18, face = "bold"),
    axis.title.x = element_text(color = "#333333", size = 14, face = "bold"),
    axis.title.y = element_text(color = "#333333", size = 14, face = "bold"),
    text = element_text(family = "Helvetica"),
    legend.position = "none"
  ) +
  scale_colour_manual(values = ggthemes::colorblind_pal()(7)[-1:-4]) +
  scale_fill_manual(values = ggthemes::colorblind_pal()(7)[-1:-4])

# Save the plot
pq1
# ggsave(file = "cigar-contig_qq.pdf", plot = pq1, width = 9, height = 4)

# Build a table of coefficients ---
d_tab <- rbind(
  as_tibble(out_blm$draws) %>% transmute(model = "LM",
    price = beta2, price_ind = NA,
    income = beta3, income_ind = NA, lambda = NA),
  as_tibble(out_bslx$draws) %>% transmute(model = "SLX",
    price = beta2, price_ind = beta78,
    income = beta3, income_ind = beta79, lambda = NA),
  as_tibble(out_bsar$draws) %>% transmute(model = "SAR",
    price = beta2, price_ind = NA,
    income = beta3, income_ind = NA, lambda = lambda_SAR),
  as_tibble(out_bsem$draws) %>% transmute(model = "SEM",
    price = beta2, price_ind = NA,
    income = beta3, income_ind = NA, lambda = lambda_SEM),
  as_tibble(out_bsdm$draws) %>% transmute(model = "SDM",
    price = beta2, price_ind = beta78,
    income = beta3, income_ind = beta79, lambda = lambda_SAR),
  as_tibble(out_bsdem$draws) %>% transmute(model = "SDEM",
    price = beta2, price_ind = beta78,
    income = beta3, income_ind = beta79, lambda = lambda_SEM)
) %>% tidyr::pivot_longer(cols = 2:6)

# Prepare the table
tbl <- d_tab %>% group_by(model, name) %>%
  summarise_all(list(mean = mean, sd = sd)) %>% # Report mean and sd
  tidyr::pivot_longer(cols = 3:4, names_to = "measure") %>%
  mutate(value = round(value, 3)) %>% # Three digits
  mutate(value = ifelse(measure == "sd",
    gsub("(.*)", "(\\1)", value), value)) %>% # Brackets around sd
  tidyr::pivot_wider(names_from = model)

# Fix order
tbl <- tbl[c(7, 8, 1, 2, 9, 10, 3, 4, 5, 6),
  c("name", "LM", "SEM", "SDEM", "SLX", "SDM", "SAR")]
# Remove NAs
tbl[is.na(tbl)] <- ""

# Export to LaTeX
tbl %>% as.data.frame() %>% memisc:::toLatex.data.frame()
