# -----------------------------
# Parameters
# -----------------------------
K <- 3
R <- 4
p <- 100
n <- 100

# Group probabilities
p1 <- c(0.3, 0.3, 0.2, 0.2)

# Signal strength grid
ol <- seq(0, 2, 0.2)

# Storage
results_array <- array(NA, dim = c(R, n, length(ol)))

# -----------------------------
# Covariance structure
# -----------------------------
pho1 <- 0
sigma1 <- matrix(pho1, p, p)
diag(sigma1) <- 1

mu <- rep(0, p)

# -----------------------------
# Main simulation loop (single run demo)
# -----------------------------
lli <- 1  # choose one signal level for example

# Sample group labels
lp <- sample(1:R, n, replace = TRUE, prob = p1)

# Group sizes
nk <- sapply(1:R, function(v) sum(lp == v))

# Mean shifts
u1 <- c(1, 1, rep(0, p - 2)) * ol[lli]
u2 <- c(2, 2, rep(0, p - 2)) * ol[lli]

# Factor model components
BX <- matrix(runif(p * K, -2, 2), nrow = p)
fX <- matrix(rnorm(K * n, 0, 1), nrow = n)

# -----------------------------
# Generate data matrix
# -----------------------------
X <- matrix(NA, n, p)

# Base noise
base_noise <- matrix(rnorm(n * p), n, p)

# Assign groups
for (i in 1:n) {
  if (lp[i] == 1) {
    X[i, ] <- base_noise[i, ]
  }
  if (lp[i] == 2) {
    X[i, ] <- base_noise[i, ] + u1
  }
  if (lp[i] == 3) {
    X[i, ] <- base_noise[i, ] + u2
  }
  if (lp[i] == 4) {
    X[i, ] <- base_noise[i, ] + 0.5 * u2
  }
}

# Add factor structure
X <- X + fX %*% t(BX)

# -----------------------------
# Construct group index list
# -----------------------------
groups <- list(
  which(lp == 1),
  which(lp == 2),
  which(lp == 3),
  which(lp == 4)
)

# -----------------------------
# Run variance test
# -----------------------------
result <- run_variance_test(X, groups)

# -----------------------------
# Output
# -----------------------------
cat("Final Test Result:\n")
print(result$final_result)

cat("\nFull Result Object:\n")
print(result)