# ============================================================
# Variance Test Statistic Construction
# ============================================================

compute_variance_statistic <- function(X_star, groups) {
  
  R <- length(groups)
  p <- ncol(X_star)
  nk <- sapply(groups, length)
  
  X_list <- lapply(groups, function(idx) X_star[idx, , drop = FALSE])
  
  sum1 <- matrix(0, R, p)
  sum2 <- array(0, c(R, R, p))
  
  # Within-group term
  for (l in 1:R) {
    for (j in 1:p) {
      vec <- X_list[[l]][, j]
      mat <- vec %*% t(vec)
      sum1[l, j] <- (sum(mat) - sum(diag(mat))) / (nk[l]*(nk[l]-1))
    }
  }
  
  # Between-group term
  for (l in 1:(R-1)) {
    for (ll in (l+1):R) {
      for (j in 1:p) {
        mat <- X_list[[l]][, j] %*% t(X_list[[ll]][, j])
        sum2[l,ll,j] <- sum(mat)/(nk[l]*nk[ll])
      }
    }
  }
  
  # Final statistic
  M_j <- rep(0, p)
  
  for (j in 1:p) {
    for (l in 1:(R-1)) {
      for (ll in (l+1):R) {
        M_j[j] <- M_j[j] + nk[l]*nk[ll] *
          (sum1[l,j] + sum1[ll,j] - sum2[l,ll,j])
      }
    }
  }
  
  # Variance (simple version)
  var_est <- rep(var(M_j), p)
  
  standardized <- M_j / sqrt(var_est)
  
  return(list(
    raw_stat = M_j,
    standardized_stat = standardized
  ))
}


# ============================================================
# Power Enhancement Step
# ============================================================

apply_power_enhancement <- function(stat, standardized_stat, X_star, groups) {
  
  p <- length(stat)
  n <- nrow(X_star)
  
  threshold <- 1.9 * log(p) * log(log(n))
  
  # Detect sparse signals
  signal <- (sqrt(2)*abs(standardized_stat) + 1) > threshold
  
  if (sum(signal) == 0) {
    
    # No enhancement needed
    final_stat <- sum(standardized_stat) / sqrt(p)
    
    return(list(
      enhanced = FALSE,
      statistic = final_stat
    ))
    
  } else {
    
    # Apply enhancement
    W <- diag(abs(standardized_stat), p, p)
    
    final_stat <- sum(standardized_stat * signal)
    
    return(list(
      enhanced = TRUE,
      statistic = final_stat,
      W = W
    ))
  }
}