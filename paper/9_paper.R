
set.seed(42)

library("bsreg")

x <- blm(log(sales) ~ log(price), data = cigarettes)

x <- blm(log(sales) ~ log(price), data = cigarettes,
  options = set_options("Conjugate",
    NG = set_NG(prec = 1e-4, shape = 1, rate = 1)))

xy <- cigarettes[cigarettes[["year"]] == 1980, c("longitude", "latitude")]
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
delta <- 2 # Decay parameter
W_decay <- dist ^ -delta
W_scaled <- W_decay / max(eigen(W_decay, symmetric = TRUE)[["values"]])
n_time <- length(unique(cigarettes[["year"]]))
W <- kronecker(diag(n_time), W_scaled) # Repeated for every year

x <- bsdm(log(sales) ~ log(price), data = cigarettes,
  W = W, options = set_options(
    SAR = set_SAR(lambda_a = 1, lambda_b = 1)),
  n_save = 5000L, n_burn = 1000L)
