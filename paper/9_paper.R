
set.seed(42)
library("bsreg")

# Code in the paper ---

# Plain linear model (Independent Normal-Gamma)
x <- blm(log(sales) ~ log(price), data = cigarettes)

# Conjugate linear model (Dependent Normal-Gamma)
x <- blm(log(sales) ~ log(price), data = cigarettes,
  options = set_options("Conjugate",
    NG = set_NG(prec = 1e-4, shape = 1, rate = 1)))

# Construct an inverse distance-decay matrix
xy <- cigarettes[cigarettes[["year"]] == 1980, c("longitude", "latitude")]
dist <- as.matrix(dist(xy))
diag(dist) <- Inf
delta <- 3 # Decay parameter
W_decay <- dist ^ -delta
W_scaled <- W_decay / max(eigen(W_decay, symmetric = TRUE)[["values"]])
n_time <- length(unique(cigarettes[["year"]]))
W <- kronecker(diag(n_time), W_scaled) # Repeated for every year

# Spatial Durbin model (Uniform lambda)
x <- bsdm(log(sales) ~ log(price), data = cigarettes,
  W = W, options = set_options(
    SAR = set_SAR(lambda_a = 1, lambda_b = 1)),
  n_save = 5000L, n_burn = 1000L)

# Extra code ---

# The results in Section 5 use the data by Halleck-Vega and Elhorst (2015),
# including their choices with respect to spatial connectivity. These results
# can be reproduced with external code attached.

# However, we can reproduce the results just with 'bsreg' to some extent

# The linear model reproduces
blm(log(sales) ~ log(price / cpi) + log(ndi / cpi) +
  factor(name) + factor(year), data = cigarettes)

# There are some discrepancies for inverse distance-decay
y <- log(cigarettes$sales)
X <- model.matrix(log(sales) ~ log(price / cpi) + log(ndi / cpi) +
  factor(name) + factor(year), data = cigarettes)
X_lag <- X[, 2:3] # To lag with the W from above
bslx(y ~ X, X_SLX = X_lag, W = W, data = cigarettes)
# Other packages are needed to construct a contiguity matrix
