# ============================================================
# Main Pipeline
# Variance Test with Factor Removal + Power Enhancement
# ============================================================

source("R/factor_number_test.R")
source("R/factor_removal.R")
source("R/variance_statistic.R")
source("R/power_enhance.R")

run_variance_test <- function(X, groups, kmax = 8, alpha = 0.05, B = 200) {
  
  cat("Step 1: Estimating number of factors...\n")
  K_hat <- estimate_factor_number(X, kmax, alpha, B)
  
  cat("Estimated number of factors:", K_hat, "\n")
  
  cat("Step 2: Removing common factors...\n")
  X_star <- remove_factor(X, K_hat)
  
  cat("Step 3: Computing variance statistic...\n")
  stat_info <- compute_variance_statistic(X_star, groups)
  
  cat("Step 4: Applying power enhancement...\n")
  final <- apply_power_enhancement(
    stat = stat_info$raw_stat,
    standardized_stat = stat_info$standardized_stat,
    X_star = X_star,
    groups = groups
  )
  
  return(list(
    K_hat = K_hat,
    X_star = X_star,
    statistic = stat_info,
    final_result = final
  ))
}