# ============================================================
# Factor Removal using Robust Huber Estimation
# ============================================================

library(psych)

remove_factor <- function(X, K_hat) {
  
  if (K_hat <= 0) return(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Step 1: Estimate factor loadings
  fa_result <- fa(X, nfactors = K_hat, fm = "minres", rotate = "varimax")
  B_hat <- as.matrix(fa_result$loadings)
  
  # Step 2: Define Huber loss
  huber_loss <- function(f, x, b, delta = 0.5) {
    r <- x - sum(b * f)
    if (abs(r) <= delta) return(0.5 * r^2)
    return(delta * (abs(r) - 0.5 * delta))
  }
  
  # Step 3: Estimate factor scores
  F_hat <- matrix(0, n, K_hat)
  
  for (i in 1:n) {
    
    objective <- function(f) {
      sum(sapply(1:p, function(j)
        huber_loss(f, X[i, j], B_hat[j, ])
      ))
    }
    
    res <- optim(rep(0, K_hat), objective, method = "L-BFGS-B")
    F_hat[i, ] <- res$par
  }
  
  # Step 4: Remove factor component
  X_hat <- F_hat %*% t(B_hat)
  
  return(X - X_hat)
}