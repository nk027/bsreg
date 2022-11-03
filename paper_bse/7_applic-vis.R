
# Generate outputs ---
library("dplyr")
library("ggplot2")
library("ggdist") # Dotplots
library("ggthemes")

load("paper/te.Rda")

# Compute total effects ---
d <- tibble(
  "SAR(2)" = te_sar1a,
  "SAR(1.4)" = te_sar1b,
  "SAR(δ)" = te_sar2
)


tidyr::pivot_longer(d, 1:4) |>
  mutate(name = factor(name, levels = c("SAR(δ)", "SAR(1)", "SAR(1.5)", "SAR(2)"))) |>
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_hline(yintercept = 0, lwd = 1.5, col = "gray") +
  stat_dots(aes(col = name), quantiles = 250,
    width = .75, justification = -0.2) + # 250 dots, narrower & shifted right
  geom_boxplot(col = "#444444", alpha = 1, # No points for outliers
    width = .2, size = .8, outlier.color = NA) +
  scale_y_continuous(limits = c(-22.5, .1)) + # Manually set limits
  # scale_x_discrete(expand = c(0, 0)) + # Manually set limits
  ylab("average partial effect") + xlab("") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 14),
    text = element_text(family = "Noto Sans"),
    # plot.margin = unit(c(0, 0, 0, .5), "cm"),
    panel.spacing = unit(20, "lines"),
    legend.position = "none"
  ) +
  scale_colour_manual(values = c("#008080", "#aaaaaa", "#aaaaaa", "#aaaaaa")) +
  scale_fill_manual(values = c("#008080", "#aaaaaa", "#aaaaaa", "#aaaaaa"))

ggsave(file = "paper/applic_te.pdf", width = 7, height = 3.1, device = cairo_pdf)



plot.ts((x2$draws[, c("beta4", "lambda_SAR", "delta_SAR")]))

op <- par(mfrow = c(2, 2))
plot(density(x2$draws[, c("beta4")]))
lines(density(x1b$draws[, c("beta4")]), lty = 2)

plot(density(x2$draws[, c("lambda_SAR")]))
lines(density(x1b$draws[, c("lambda_SAR")]), lty = 2)

plot(density(x2$draws[, c("delta_SAR")]))
par(op)

quantile(x2$draws[, c("delta_SAR")], c(.1, .9))

sum(qr.solve(diag(n_reg) - mean(x2$draws[, "lambda_SAR"]) *
  Psi(1.333)[seq(n_reg), seq(n_reg)])) / n_reg
sum(qr.solve(diag(n_reg) - mean(x2$draws[, "lambda_SAR"]) *
  Psi(1.5)[seq(n_reg), seq(n_reg)])) / n_reg
