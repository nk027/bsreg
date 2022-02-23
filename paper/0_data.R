
# Prepare data ---
library("dplyr")
library("readxl")

# Data kindly made available by Paul Elhorst at <https://spatial-panels.com/>
download.file(
  "https://spatial-panels.com/wp-content/uploads/2017/06/Files-SLX-paper.zip",
  destfile = "paper/data/halleckvega2015.zip")
unzip("paper/data/halleckvega2015.zip", exdir = "paper/data",
  files = c("cigarette+2var.xls", "Spat-Sym-US.xls", "cigar_states.xls"))

# Prepare the data
df <- read_excel("paper/data/cigarette+2var.xls") %>%
  mutate(year = factor(year + 1963), state = factor(state))
contig <- read_excel("paper/data/Spat-Sym-US.xls",
  col_names = paste0(1:46)) %>% as.matrix()
xy <- read_excel("paper/data/cigar_states.xls", col_names = TRUE)

# Utah has a wrong longitude
xy %>% filter(Name == "UTAH") # 11.9 is off the Irish coast
xy <- xy %>% mutate(longitude = ifelse(Name == "UTAH", 111.7, longitude))

# Connectivities ---
n_time <- length(unique(df$year))
# Contiguity matrix
W_cont <- kronecker(diag(n_time), contig / rowSums(contig))
# Inverse-distance decay function
dist <- as.matrix(dist(xy))
diag(dist) <- Inf # Diagonal elements will be 0
Psi <- function(delta) {
  W_dist <- dist ^ (-delta) # Build
  W_dist <- W_dist / max(eigen(W_dist, symmetric = TRUE)$values) # Scale
  kronecker(diag(n_time), W_dist)
}
W_dist <- Psi(3)

# Prepare variables ---
y <- cbind(logc = df$logc)
X <- model.matrix(logc ~ logp + logy + year + state, data = df)
X_lag <- X[, c("logp", "logy")]
colnames(X_lag) <- c("wlogp", "wlogy")
X_cont <- W_cont %*% X_lag
