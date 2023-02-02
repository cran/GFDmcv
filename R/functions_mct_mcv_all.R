
# rr

# test statistics MCT
# x - a list of length k with elements being n_i\times d matrices of data
mct_test_stat_rr <- function(x, h) {
  kk <- length(x)
  nn <- sum(sapply(x, nrow))
  p <- ncol(x[[1]])
  r <- nrow(h)
  
  mcv_c_est_temp <- sapply(x, mcv_est_rr)
  mcv_b_est_temp <- sapply(x, function(y) lapply(mcv_est_rr(y), function(z) 1 / z))
  
  sigma2_cb_est_mat <- sapply(x, sigma2_est_rr, nn = nn)
  sigma_rr_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_rr_b_est <- diag(sigma2_cb_est_mat[2, ])
  
  if (r > 1) {
    rr_c_stat <- sqrt(nn) * max(abs(diag((diag(h %*% sigma_rr_c_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_c_est_temp))))
    rr_b_stat <- sqrt(nn) * max(abs(diag((diag(h %*% sigma_rr_b_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_b_est_temp))))
  } else if (r == 1) {
    rr_c_stat <- sqrt(nn) * max(abs((diag(h %*% sigma_rr_c_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_c_est_temp))))
    rr_b_stat <- sqrt(nn) * max(abs((diag(h %*% sigma_rr_b_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_b_est_temp))))
  }
  
  return(list(rr_c_stat, rr_b_stat,
              sigma_rr_c_est, sigma_rr_b_est))
}

# test statistic for each contrast separately
mct_test_stat_sep_rr <- function(x, h) {
  kk <- length(x)
  nn <- sum(sapply(x, nrow))
  p <- ncol(x[[1]])
  r <- nrow(h)
  
  mcv_c_est_temp <- sapply(x, mcv_est_rr)
  mcv_b_est_temp <- sapply(x, function(y) lapply(mcv_est_rr(y), function(z) 1 / z))
  
  sigma2_cb_est_mat <- sapply(x, sigma2_est_rr, nn = nn)
  sigma_rr_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_rr_b_est <- diag(sigma2_cb_est_mat[2, ])
  
  if (r > 1) {
    rr_c_stat <- sqrt(nn) * abs(diag((diag(h %*% sigma_rr_c_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_c_est_temp)))
    rr_b_stat <- sqrt(nn) * abs(diag((diag(h %*% sigma_rr_b_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_b_est_temp)))
  } else if (r == 1) {
    rr_c_stat <- sqrt(nn) * abs((diag(h %*% sigma_rr_c_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_c_est_temp)))
    rr_b_stat <- sqrt(nn) * abs((diag(h %*% sigma_rr_b_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_b_est_temp)))
  }
  
  return(list(rr_c_stat[, 1], rr_b_stat[, 1]))
}

mct_asympt_test_rr <- function(x, h, alpha = 0.05) {
  r <- nrow(h)
  temp <- mct_test_stat_rr(x, h)
  sigma_rr_c_est <- temp[[3]]
  sigma_rr_b_est <- temp[[4]]
  d_sigma_d_rr_c <- h %*% sigma_rr_c_est %*% t(h)
  d_sigma_d_rr_b <- h %*% sigma_rr_b_est %*% t(h)
  
  if (r > 1) {
    asympt_cov_rr_c <- diag((diag(d_sigma_d_rr_c))^(-0.5)) %*% d_sigma_d_rr_c %*% diag((diag(d_sigma_d_rr_c))^(-0.5))
    asympt_cov_rr_b <- diag((diag(d_sigma_d_rr_b))^(-0.5)) %*% d_sigma_d_rr_b %*% diag((diag(d_sigma_d_rr_b))^(-0.5))
  } else if (r == 1) {
    asympt_cov_rr_c <- (diag(d_sigma_d_rr_c))^(-0.5) %*% d_sigma_d_rr_c %*% (diag(d_sigma_d_rr_c))^(-0.5)
    asympt_cov_rr_b <- (diag(d_sigma_d_rr_b))^(-0.5) %*% d_sigma_d_rr_b %*% (diag(d_sigma_d_rr_b))^(-0.5)
  }
  
  p_val_c <- mvtnorm::pmvnorm(lower = rep(-Inf, nrow(asympt_cov_rr_c)), upper = temp[[1]], sigma = asympt_cov_rr_c)
  p_val_b <- mvtnorm::pmvnorm(lower = rep(-Inf, nrow(asympt_cov_rr_c)), upper = temp[[2]], sigma = asympt_cov_rr_b)
  
  return(c(2 * min(p_val_c, 1 - p_val_c),
           2 * min(p_val_b, 1 - p_val_b),
           mvtnorm::qmvnorm(1 - alpha, sigma = asympt_cov_rr_c, tail = "both.tails")$quantile,
           mvtnorm::qmvnorm(1 - alpha, sigma = asympt_cov_rr_b, tail = "both.tails")$quantile))
}

# vv

# test statistics MCT
# x - a list of length k with elements being n_i\times d matrices of data
mct_test_stat_vv <- function(x, h) {
  kk <- length(x)
  nn <- sum(sapply(x, nrow))
  p <- ncol(x[[1]])
  r <- nrow(h)
  
  mcv_c_est_temp <- sapply(x, mcv_est_vv)
  mcv_b_est_temp <- sapply(x, function(y) lapply(mcv_est_vv(y), function(z) 1 / z))
  
  sigma2_cb_est_mat <- sapply(x, sigma2_est_vv, nn = nn)
  sigma_vv_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_vv_b_est <- diag(sigma2_cb_est_mat[2, ])
  
  if (r > 1) {
    vv_c_stat <- sqrt(nn) * max(abs(diag((diag(h %*% sigma_vv_c_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_c_est_temp))))
    vv_b_stat <- sqrt(nn) * max(abs(diag((diag(h %*% sigma_vv_b_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_b_est_temp))))
  } else if (r == 1) {
    vv_c_stat <- sqrt(nn) * max(abs((diag(h %*% sigma_vv_c_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_c_est_temp))))
    vv_b_stat <- sqrt(nn) * max(abs((diag(h %*% sigma_vv_b_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_b_est_temp))))
  }
  
  return(list(vv_c_stat, vv_b_stat,
              sigma_vv_c_est, sigma_vv_b_est))
}

# test statistic for each contrast separately
mct_test_stat_sep_vv <- function(x, h) {
  kk <- length(x)
  nn <- sum(sapply(x, nrow))
  p <- ncol(x[[1]])
  r <- nrow(h)
  
  mcv_c_est_temp <- sapply(x, mcv_est_vv)
  mcv_b_est_temp <- sapply(x, function(y) lapply(mcv_est_vv(y), function(z) 1 / z))
  
  sigma2_cb_est_mat <- sapply(x, sigma2_est_vv, nn = nn)
  sigma_vv_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_vv_b_est <- diag(sigma2_cb_est_mat[2, ])
  
  if (r > 1) {
    vv_c_stat <- sqrt(nn) * abs(diag((diag(h %*% sigma_vv_c_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_c_est_temp)))
    vv_b_stat <- sqrt(nn) * abs(diag((diag(h %*% sigma_vv_b_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_b_est_temp)))
  } else if (r == 1) {
    vv_c_stat <- sqrt(nn) * abs((diag(h %*% sigma_vv_c_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_c_est_temp)))
    vv_b_stat <- sqrt(nn) * abs((diag(h %*% sigma_vv_b_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_b_est_temp)))
  }
  
  return(list(vv_c_stat[, 1], vv_b_stat[, 1]))
}

mct_asympt_test_vv <- function(x, h, alpha = 0.05) {
  r <- nrow(h)
  temp <- mct_test_stat_vv(x, h)
  sigma_vv_c_est <- temp[[3]]
  sigma_vv_b_est <- temp[[4]]
  d_sigma_d_vv_c <- h %*% sigma_vv_c_est %*% t(h)
  d_sigma_d_vv_b <- h %*% sigma_vv_b_est %*% t(h)
  
  if (r > 1) {
    asympt_cov_vv_c <- diag((diag(d_sigma_d_vv_c))^(-0.5)) %*% d_sigma_d_vv_c %*% diag((diag(d_sigma_d_vv_c))^(-0.5))
    asympt_cov_vv_b <- diag((diag(d_sigma_d_vv_b))^(-0.5)) %*% d_sigma_d_vv_b %*% diag((diag(d_sigma_d_vv_b))^(-0.5))
  } else if (r == 1) {
    asympt_cov_vv_c <- (diag(d_sigma_d_vv_c))^(-0.5) %*% d_sigma_d_vv_c %*% (diag(d_sigma_d_vv_c))^(-0.5)
    asympt_cov_vv_b <- (diag(d_sigma_d_vv_b))^(-0.5) %*% d_sigma_d_vv_b %*% (diag(d_sigma_d_vv_b))^(-0.5)
  }
  
  p_val_c <- mvtnorm::pmvnorm(lower = rep(-Inf, nrow(asympt_cov_vv_c)), upper = temp[[1]], sigma = asympt_cov_vv_c)
  p_val_b <- mvtnorm::pmvnorm(lower = rep(-Inf, nrow(asympt_cov_vv_c)), upper = temp[[2]], sigma = asympt_cov_vv_b)
  
  return(c(2 * min(p_val_c, 1 - p_val_c),
           2 * min(p_val_b, 1 - p_val_b),
           mvtnorm::qmvnorm(1 - alpha, sigma = asympt_cov_vv_c, tail = "both.tails")$quantile,
           mvtnorm::qmvnorm(1 - alpha, sigma = asympt_cov_vv_b, tail = "both.tails")$quantile))
}

# vn

# test statistics MCT
# x - a list of length k with elements being n_i\times d matrices of data
mct_test_stat_vn <- function(x, h) {
  kk <- length(x)
  nn <- sum(sapply(x, nrow))
  p <- ncol(x[[1]])
  r <- nrow(h)
  
  mcv_c_est_temp <- sapply(x, mcv_est_vn)
  mcv_b_est_temp <- sapply(x, function(y) lapply(mcv_est_vn(y), function(z) 1 / z))
  
  sigma2_cb_est_mat <- sapply(x, sigma2_est_vn, nn = nn)
  sigma_vn_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_vn_b_est <- diag(sigma2_cb_est_mat[2, ])
  
  if (r > 1) {
    vn_c_stat <- sqrt(nn) * max(abs(diag((diag(h %*% sigma_vn_c_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_c_est_temp))))
    vn_b_stat <- sqrt(nn) * max(abs(diag((diag(h %*% sigma_vn_b_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_b_est_temp))))
  } else if (r == 1) {
    vn_c_stat <- sqrt(nn) * max(abs((diag(h %*% sigma_vn_c_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_c_est_temp))))
    vn_b_stat <- sqrt(nn) * max(abs((diag(h %*% sigma_vn_b_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_b_est_temp))))
  }
  
  return(list(vn_c_stat, vn_b_stat,
              sigma_vn_c_est, sigma_vn_b_est))
}

# test statistic for each contrast separately
mct_test_stat_sep_vn <- function(x, h) {
  kk <- length(x)
  nn <- sum(sapply(x, nrow))
  p <- ncol(x[[1]])
  r <- nrow(h)
  
  mcv_c_est_temp <- sapply(x, mcv_est_vn)
  mcv_b_est_temp <- sapply(x, function(y) lapply(mcv_est_vn(y), function(z) 1 / z))
  
  sigma2_cb_est_mat <- sapply(x, sigma2_est_vn, nn = nn)
  sigma_vn_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_vn_b_est <- diag(sigma2_cb_est_mat[2, ])
  
  if (r > 1) {
    vn_c_stat <- sqrt(nn) * abs(diag((diag(h %*% sigma_vn_c_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_c_est_temp)))
    vn_b_stat <- sqrt(nn) * abs(diag((diag(h %*% sigma_vn_b_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_b_est_temp)))
  } else if (r == 1) {
    vn_c_stat <- sqrt(nn) * abs((diag(h %*% sigma_vn_c_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_c_est_temp)))
    vn_b_stat <- sqrt(nn) * abs((diag(h %*% sigma_vn_b_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_b_est_temp)))
  }
  
  return(list(vn_c_stat[, 1], vn_b_stat[, 1]))
}

mct_asympt_test_vn <- function(x, h, alpha = 0.05) {
  r <- nrow(h)
  temp <- mct_test_stat_vn(x, h)
  sigma_vn_c_est <- temp[[3]]
  sigma_vn_b_est <- temp[[4]]
  d_sigma_d_vn_c <- h %*% sigma_vn_c_est %*% t(h)
  d_sigma_d_vn_b <- h %*% sigma_vn_b_est %*% t(h)
  
  if (r > 1) {
    asympt_cov_vn_c <- diag((diag(d_sigma_d_vn_c))^(-0.5)) %*% d_sigma_d_vn_c %*% diag((diag(d_sigma_d_vn_c))^(-0.5))
    asympt_cov_vn_b <- diag((diag(d_sigma_d_vn_b))^(-0.5)) %*% d_sigma_d_vn_b %*% diag((diag(d_sigma_d_vn_b))^(-0.5))
  } else if (r == 1) {
    asympt_cov_vn_c <- (diag(d_sigma_d_vn_c))^(-0.5) %*% d_sigma_d_vn_c %*% (diag(d_sigma_d_vn_c))^(-0.5)
    asympt_cov_vn_b <- (diag(d_sigma_d_vn_b))^(-0.5) %*% d_sigma_d_vn_b %*% (diag(d_sigma_d_vn_b))^(-0.5)
  }
  
  p_val_c <- mvtnorm::pmvnorm(lower = rep(-Inf, nrow(asympt_cov_vn_c)), upper = temp[[1]], sigma = asympt_cov_vn_c)
  p_val_b <- mvtnorm::pmvnorm(lower = rep(-Inf, nrow(asympt_cov_vn_c)), upper = temp[[2]], sigma = asympt_cov_vn_b)
  
  return(c(2 * min(p_val_c, 1 - p_val_c),
           2 * min(p_val_b, 1 - p_val_b),
           mvtnorm::qmvnorm(1 - alpha, sigma = asympt_cov_vn_c, tail = "both.tails")$quantile,
           mvtnorm::qmvnorm(1 - alpha, sigma = asympt_cov_vn_b, tail = "both.tails")$quantile))
}

# az

# test statistics MCT
# x - a list of length k with elements being n_i\times d matrices of data
mct_test_stat_az <- function(x, h) {
  kk <- length(x)
  nn <- sum(sapply(x, nrow))
  p <- ncol(x[[1]])
  r <- nrow(h)
  
  mcv_c_est_temp <- sapply(x, mcv_est_az)
  mcv_b_est_temp <- sapply(x, function(y) lapply(mcv_est_az(y), function(z) 1 / z))
  
  sigma2_cb_est_mat <- sapply(x, sigma2_est_az, nn = nn)
  sigma_az_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_az_b_est <- diag(sigma2_cb_est_mat[2, ])
  
  if (r > 1) {
    az_c_stat <- sqrt(nn) * max(abs(diag((diag(h %*% sigma_az_c_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_c_est_temp))))
    az_b_stat <- sqrt(nn) * max(abs(diag((diag(h %*% sigma_az_b_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_b_est_temp))))
  } else if (r == 1) {
    az_c_stat <- sqrt(nn) * max(abs((diag(h %*% sigma_az_c_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_c_est_temp))))
    az_b_stat <- sqrt(nn) * max(abs((diag(h %*% sigma_az_b_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_b_est_temp))))
  }
  
  return(list(az_c_stat, az_b_stat,
              sigma_az_c_est, sigma_az_b_est))
}

# test statistic for each contrast separately
mct_test_stat_sep_az <- function(x, h) {
  kk <- length(x)
  nn <- sum(sapply(x, nrow))
  p <- ncol(x[[1]])
  r <- nrow(h)
  
  mcv_c_est_temp <- sapply(x, mcv_est_az)
  mcv_b_est_temp <- sapply(x, function(y) lapply(mcv_est_az(y), function(z) 1 / z))
  
  sigma2_cb_est_mat <- sapply(x, sigma2_est_az, nn = nn)
  sigma_az_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_az_b_est <- diag(sigma2_cb_est_mat[2, ])
  
  if (r > 1) {
    az_c_stat <- sqrt(nn) * abs(diag((diag(h %*% sigma_az_c_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_c_est_temp)))
    az_b_stat <- sqrt(nn) * abs(diag((diag(h %*% sigma_az_b_est %*% t(h)))^(-0.5)) %*% h %*% diag(diag(mcv_b_est_temp)))
  } else if (r == 1) {
    az_c_stat <- sqrt(nn) * abs((diag(h %*% sigma_az_c_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_c_est_temp)))
    az_b_stat <- sqrt(nn) * abs((diag(h %*% sigma_az_b_est %*% t(h)))^(-0.5) %*% h %*% diag(diag(mcv_b_est_temp)))
  }
  
  return(list(az_c_stat[, 1], az_b_stat[, 1]))
}

mct_asympt_test_az <- function(x, h, alpha = 0.05) {
  r <- nrow(h)
  temp <- mct_test_stat_az(x, h)
  sigma_az_c_est <- temp[[3]]
  sigma_az_b_est <- temp[[4]]
  d_sigma_d_az_c <- h %*% sigma_az_c_est %*% t(h)
  d_sigma_d_az_b <- h %*% sigma_az_b_est %*% t(h)
  
  if (r > 1) {
    asympt_cov_az_c <- diag((diag(d_sigma_d_az_c))^(-0.5)) %*% d_sigma_d_az_c %*% diag((diag(d_sigma_d_az_c))^(-0.5))
    asympt_cov_az_b <- diag((diag(d_sigma_d_az_b))^(-0.5)) %*% d_sigma_d_az_b %*% diag((diag(d_sigma_d_az_b))^(-0.5))
  } else if (r == 1) {
    asympt_cov_az_c <- (diag(d_sigma_d_az_c))^(-0.5) %*% d_sigma_d_az_c %*% (diag(d_sigma_d_az_c))^(-0.5)
    asympt_cov_az_b <- (diag(d_sigma_d_az_b))^(-0.5) %*% d_sigma_d_az_b %*% (diag(d_sigma_d_az_b))^(-0.5)
  }
  
  p_val_c <- mvtnorm::pmvnorm(lower = rep(-Inf, nrow(asympt_cov_az_c)), upper = temp[[1]], sigma = asympt_cov_az_c)
  p_val_b <- mvtnorm::pmvnorm(lower = rep(-Inf, nrow(asympt_cov_az_c)), upper = temp[[2]], sigma = asympt_cov_az_b)
  
  return(c(2 * min(p_val_c, 1 - p_val_c),
           2 * min(p_val_b, 1 - p_val_b),
           mvtnorm::qmvnorm(1 - alpha, sigma = asympt_cov_az_c, tail = "both.tails")$quantile,
           mvtnorm::qmvnorm(1 - alpha, sigma = asympt_cov_az_b, tail = "both.tails")$quantile))
}

# pooled bootstrap for all mcvs
mct_pool_boot <- function(x, h, n_boot = 1000, alpha = 0.05, 
                          parallel = FALSE, n_cores = NULL) {
  temp_rr <- mct_test_stat_rr(x, h)
  sigma_rr_c_est <- temp_rr[[3]]
  sigma_rr_b_est <- temp_rr[[4]]
  temp_vv <- mct_test_stat_vv(x, h)
  sigma_vv_c_est <- temp_vv[[3]]
  sigma_vv_b_est <- temp_vv[[4]]
  temp_vn <- mct_test_stat_vn(x, h)
  sigma_vn_c_est <- temp_vn[[3]]
  sigma_vn_b_est <- temp_vn[[4]]
  temp_az <- mct_test_stat_az(x, h)
  sigma_az_c_est <- temp_az[[3]]
  sigma_az_b_est <- temp_az[[4]]
  
  kk <- length(x)
  n_i <- sapply(x, nrow)
  nn <- sum(n_i)
  p <- ncol(x[[1]])
  r <- nrow(h)
  
  if (p == 1) {
    # stop("the dimension is equal to 1")
  }
  
  x_mat <- x[[1]]
  for (i_kk in 2:kk) {
    x_mat <- rbind(x_mat, x[[i_kk]])
  }
  c_pool_rr <- mcv_est_rr(x_mat)$rr
  b_pool_rr <- 1 / c_pool_rr
  c_pool_vv <- mcv_est_vv(x_mat)$vv
  b_pool_vv <- 1 / c_pool_vv
  c_pool_vn <- mcv_est_vn(x_mat)$vn
  b_pool_vn <- 1 / c_pool_vn
  c_pool_az <- mcv_est_az(x_mat)$az
  b_pool_az <- 1 / c_pool_az
  
  groups <- rep(1:kk, n_i)
  if (!parallel) {
    rr_c_stat <- numeric(n_boot)
    rr_b_stat <- numeric(n_boot)
    vv_c_stat <- numeric(n_boot)
    vv_b_stat <- numeric(n_boot)
    vn_c_stat <- numeric(n_boot)
    vn_b_stat <- numeric(n_boot)
    az_c_stat <- numeric(n_boot)
    az_b_stat <- numeric(n_boot)
    for (i_b in seq_len(n_boot)) {
      if (p > 1) {
        x_mat_boot <- x_mat[sample(nn, replace = TRUE), ]
      } else if (p == 1) {
        x_mat_boot <- matrix(x_mat[sample(nn, replace = TRUE), ], ncol = 1)
      }
      x_boot_list <- vector("list", kk)
      for (j_kk in 1:kk) {
        if (p > 1) {
          x_boot_list[[j_kk]] <- x_mat_boot[groups == j_kk, ]
        } else if (p == 1) {
          x_boot_list[[j_kk]] <- matrix(x_mat_boot[groups == j_kk, ], ncol = 1)
        }
      }
      mcv_c_est_temp_rr <- unlist(sapply(x_boot_list, mcv_est_rr))
      mcv_c_est_temp_vv <- unlist(sapply(x_boot_list, mcv_est_vv))
      mcv_c_est_temp_vn <- unlist(sapply(x_boot_list, mcv_est_vn))
      mcv_c_est_temp_az <- unlist(sapply(x_boot_list, mcv_est_az))
      mcv_b_est_temp_rr <- unlist(sapply(x_boot_list, function(y) lapply(mcv_est_rr(y), function(z) 1 / z)))
      mcv_b_est_temp_vv <- unlist(sapply(x_boot_list, function(y) lapply(mcv_est_vv(y), function(z) 1 / z)))
      mcv_b_est_temp_vn <- unlist(sapply(x_boot_list, function(y) lapply(mcv_est_vn(y), function(z) 1 / z)))
      mcv_b_est_temp_az <- unlist(sapply(x_boot_list, function(y) lapply(mcv_est_az(y), function(z) 1 / z)))
      temp_boot_rr <- mct_test_stat_rr(x_boot_list, h)
      temp_boot_vv <- mct_test_stat_vv(x_boot_list, h)
      temp_boot_vn <- mct_test_stat_vn(x_boot_list, h)
      temp_boot_az <- mct_test_stat_az(x_boot_list, h)
      sigma_rr_c_est_b <- temp_boot_rr[[3]]
      sigma_vv_c_est_b <- temp_boot_vv[[3]]
      sigma_vn_c_est_b <- temp_boot_vn[[3]]
      sigma_az_c_est_b <- temp_boot_az[[3]]
      sigma_rr_b_est_b <- temp_boot_rr[[4]]
      sigma_vv_b_est_b <- temp_boot_vv[[4]]
      sigma_vn_b_est_b <- temp_boot_vn[[4]]
      sigma_az_b_est_b <- temp_boot_az[[4]]
      sigma_rr_c_est_b_m05 <- diag(1 / sqrt(diag(sigma_rr_c_est_b)))
      sigma_vv_c_est_b_m05 <- diag(1 / sqrt(diag(sigma_vv_c_est_b)))
      sigma_vn_c_est_b_m05 <- diag(1 / sqrt(diag(sigma_vn_c_est_b)))
      sigma_az_c_est_b_m05 <- diag(1 / sqrt(diag(sigma_az_c_est_b)))
      sigma_rr_b_est_b_m05 <- diag(1 / sqrt(diag(sigma_rr_b_est_b)))
      sigma_vv_b_est_b_m05 <- diag(1 / sqrt(diag(sigma_vv_b_est_b)))
      sigma_vn_b_est_b_m05 <- diag(1 / sqrt(diag(sigma_vn_b_est_b)))
      sigma_az_b_est_b_m05 <- diag(1 / sqrt(diag(sigma_az_b_est_b)))
      rr_c_stat_temp <- numeric(r)
      vv_c_stat_temp <- numeric(r)
      vn_c_stat_temp <- numeric(r)
      az_c_stat_temp <- numeric(r)
      rr_b_stat_temp <- numeric(r)
      vv_b_stat_temp <- numeric(r)
      vn_b_stat_temp <- numeric(r)
      az_b_stat_temp <- numeric(r)
      for (i in seq_len(r)) {
        rr_c_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_rr_c_est) %*% sigma_rr_c_est_b_m05 %*% (mcv_c_est_temp_rr - c_pool_rr)) /
                                sqrt(t(h[i, ]) %*% sigma_rr_c_est %*% h[i, ]))[1, 1]
        rr_b_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_rr_b_est) %*% sigma_rr_b_est_b_m05 %*% (mcv_b_est_temp_rr - b_pool_rr)) /
                                sqrt(t(h[i, ]) %*% sigma_rr_b_est %*% h[i, ]))[1, 1]
        vv_c_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_vv_c_est) %*% sigma_vv_c_est_b_m05 %*% (mcv_c_est_temp_vv - c_pool_vv)) /
                                sqrt(t(h[i, ]) %*% sigma_vv_c_est %*% h[i, ]))[1, 1]
        vv_b_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_vv_b_est) %*% sigma_vv_b_est_b_m05 %*% (mcv_b_est_temp_vv - b_pool_vv)) /
                                sqrt(t(h[i, ]) %*% sigma_vv_b_est %*% h[i, ]))[1, 1]
        vn_c_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_vn_c_est) %*% sigma_vn_c_est_b_m05 %*% (mcv_c_est_temp_vn - c_pool_vn)) /
                                sqrt(t(h[i, ]) %*% sigma_vn_c_est %*% h[i, ]))[1, 1]
        vn_b_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_vn_b_est) %*% sigma_vn_b_est_b_m05 %*% (mcv_b_est_temp_vn - b_pool_vn)) /
                                sqrt(t(h[i, ]) %*% sigma_vn_b_est %*% h[i, ]))[1, 1]
        az_c_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_az_c_est) %*% sigma_az_c_est_b_m05 %*% (mcv_c_est_temp_az - c_pool_az)) /
                                sqrt(t(h[i, ]) %*% sigma_az_c_est %*% h[i, ]))[1, 1]
        az_b_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_az_b_est) %*% sigma_az_b_est_b_m05 %*% (mcv_b_est_temp_az - b_pool_az)) /
                                sqrt(t(h[i, ]) %*% sigma_az_b_est %*% h[i, ]))[1, 1]
      }
      rr_c_stat[i_b] <- sqrt(nn) * max(abs(rr_c_stat_temp))
      rr_b_stat[i_b] <- sqrt(nn) * max(abs(rr_b_stat_temp))
      vv_c_stat[i_b] <- sqrt(nn) * max(abs(vv_c_stat_temp))
      vv_b_stat[i_b] <- sqrt(nn) * max(abs(vv_b_stat_temp))
      vn_c_stat[i_b] <- sqrt(nn) * max(abs(vn_c_stat_temp))
      vn_b_stat[i_b] <- sqrt(nn) * max(abs(vn_b_stat_temp))
      az_c_stat[i_b] <- sqrt(nn) * max(abs(az_c_stat_temp))
      az_b_stat[i_b] <- sqrt(nn) * max(abs(az_b_stat_temp))
    }
  } else {
    RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion")
    wyniki_stat <- foreach(i_n_boot = 1:n_boot, .combine = rbind, 
                           .packages = c("Rcpp", "MASS", "mvtnorm", "GFDmcv")) %dopar%
      {
        if (p > 1) {
          x_mat_boot <- x_mat[sample(nn, replace = TRUE), ]
        } else if (p == 1) {
          x_mat_boot <- matrix(x_mat[sample(nn, replace = TRUE), ], ncol = 1)
        }
        x_boot_list <- vector("list", kk)
        for (j_kk in 1:kk) {
          if (p > 1) {
            x_boot_list[[j_kk]] <- x_mat_boot[groups == j_kk, ]
          } else if (p == 1) {
            x_boot_list[[j_kk]] <- matrix(x_mat_boot[groups == j_kk, ], ncol = 1)
          }
        }
        mcv_c_est_temp_rr <- unlist(sapply(x_boot_list, mcv_est_rr))
        mcv_c_est_temp_vv <- unlist(sapply(x_boot_list, mcv_est_vv))
        mcv_c_est_temp_vn <- unlist(sapply(x_boot_list, mcv_est_vn))
        mcv_c_est_temp_az <- unlist(sapply(x_boot_list, mcv_est_az))
        mcv_b_est_temp_rr <- unlist(sapply(x_boot_list, function(y) lapply(mcv_est_rr(y), function(z) 1 / z)))
        mcv_b_est_temp_vv <- unlist(sapply(x_boot_list, function(y) lapply(mcv_est_vv(y), function(z) 1 / z)))
        mcv_b_est_temp_vn <- unlist(sapply(x_boot_list, function(y) lapply(mcv_est_vn(y), function(z) 1 / z)))
        mcv_b_est_temp_az <- unlist(sapply(x_boot_list, function(y) lapply(mcv_est_az(y), function(z) 1 / z)))
        temp_boot_rr <- mct_test_stat_rr(x_boot_list, h)
        temp_boot_vv <- mct_test_stat_vv(x_boot_list, h)
        temp_boot_vn <- mct_test_stat_vn(x_boot_list, h)
        temp_boot_az <- mct_test_stat_az(x_boot_list, h)
        sigma_rr_c_est_b <- temp_boot_rr[[3]]
        sigma_vv_c_est_b <- temp_boot_vv[[3]]
        sigma_vn_c_est_b <- temp_boot_vn[[3]]
        sigma_az_c_est_b <- temp_boot_az[[3]]
        sigma_rr_b_est_b <- temp_boot_rr[[4]]
        sigma_vv_b_est_b <- temp_boot_vv[[4]]
        sigma_vn_b_est_b <- temp_boot_vn[[4]]
        sigma_az_b_est_b <- temp_boot_az[[4]]
        sigma_rr_c_est_b_m05 <- diag(1 / sqrt(diag(sigma_rr_c_est_b)))
        sigma_vv_c_est_b_m05 <- diag(1 / sqrt(diag(sigma_vv_c_est_b)))
        sigma_vn_c_est_b_m05 <- diag(1 / sqrt(diag(sigma_vn_c_est_b)))
        sigma_az_c_est_b_m05 <- diag(1 / sqrt(diag(sigma_az_c_est_b)))
        sigma_rr_b_est_b_m05 <- diag(1 / sqrt(diag(sigma_rr_b_est_b)))
        sigma_vv_b_est_b_m05 <- diag(1 / sqrt(diag(sigma_vv_b_est_b)))
        sigma_vn_b_est_b_m05 <- diag(1 / sqrt(diag(sigma_vn_b_est_b)))
        sigma_az_b_est_b_m05 <- diag(1 / sqrt(diag(sigma_az_b_est_b)))
        rr_c_stat_temp <- numeric(r)
        vv_c_stat_temp <- numeric(r)
        vn_c_stat_temp <- numeric(r)
        az_c_stat_temp <- numeric(r)
        rr_b_stat_temp <- numeric(r)
        vv_b_stat_temp <- numeric(r)
        vn_b_stat_temp <- numeric(r)
        az_b_stat_temp <- numeric(r)
        for (i in seq_len(r)) {
          rr_c_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_rr_c_est) %*% sigma_rr_c_est_b_m05 %*% (mcv_c_est_temp_rr - c_pool_rr)) /
                                  sqrt(t(h[i, ]) %*% sigma_rr_c_est %*% h[i, ]))[1, 1]
          rr_b_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_rr_b_est) %*% sigma_rr_b_est_b_m05 %*% (mcv_b_est_temp_rr - b_pool_rr)) /
                                  sqrt(t(h[i, ]) %*% sigma_rr_b_est %*% h[i, ]))[1, 1]
          vv_c_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_vv_c_est) %*% sigma_vv_c_est_b_m05 %*% (mcv_c_est_temp_vv - c_pool_vv)) /
                                  sqrt(t(h[i, ]) %*% sigma_vv_c_est %*% h[i, ]))[1, 1]
          vv_b_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_vv_b_est) %*% sigma_vv_b_est_b_m05 %*% (mcv_b_est_temp_vv - b_pool_vv)) /
                                  sqrt(t(h[i, ]) %*% sigma_vv_b_est %*% h[i, ]))[1, 1]
          vn_c_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_vn_c_est) %*% sigma_vn_c_est_b_m05 %*% (mcv_c_est_temp_vn - c_pool_vn)) /
                                  sqrt(t(h[i, ]) %*% sigma_vn_c_est %*% h[i, ]))[1, 1]
          vn_b_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_vn_b_est) %*% sigma_vn_b_est_b_m05 %*% (mcv_b_est_temp_vn - b_pool_vn)) /
                                  sqrt(t(h[i, ]) %*% sigma_vn_b_est %*% h[i, ]))[1, 1]
          az_c_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_az_c_est) %*% sigma_az_c_est_b_m05 %*% (mcv_c_est_temp_az - c_pool_az)) /
                                  sqrt(t(h[i, ]) %*% sigma_az_c_est %*% h[i, ]))[1, 1]
          az_b_stat_temp[i] <- ((t(h[i, ]) %*% sqrt(sigma_az_b_est) %*% sigma_az_b_est_b_m05 %*% (mcv_b_est_temp_az - b_pool_az)) /
                                  sqrt(t(h[i, ]) %*% sigma_az_b_est %*% h[i, ]))[1, 1]
        }
        c(sqrt(nn) * max(abs(rr_c_stat_temp)),
          sqrt(nn) * max(abs(rr_b_stat_temp)),
          sqrt(nn) * max(abs(vv_c_stat_temp)),
          sqrt(nn) * max(abs(vv_b_stat_temp)),
          sqrt(nn) * max(abs(vn_c_stat_temp)),
          sqrt(nn) * max(abs(vn_b_stat_temp)),
          sqrt(nn) * max(abs(az_c_stat_temp)),
          sqrt(nn) * max(abs(az_b_stat_temp)))
      }
    
    rr_c_stat <- wyniki_stat[, 1]
    rr_b_stat <- wyniki_stat[, 2]
    vv_c_stat <- wyniki_stat[, 3]
    vv_b_stat <- wyniki_stat[, 4]
    vn_c_stat <- wyniki_stat[, 5]
    vn_b_stat <- wyniki_stat[, 6]
    az_c_stat <- wyniki_stat[, 7]
    az_b_stat <- wyniki_stat[, 8]
  }
  
  p_val_rr_c <- mean(rr_c_stat > temp_rr[[1]])
  p_val_rr_b <- mean(rr_b_stat > temp_rr[[2]])
  p_val_vv_c <- mean(vv_c_stat > temp_vv[[1]])
  p_val_vv_b <- mean(vv_b_stat > temp_vv[[2]])
  p_val_vn_c <- mean(vn_c_stat > temp_vn[[1]])
  p_val_vn_b <- mean(vn_b_stat > temp_vn[[2]])
  p_val_az_c <- mean(az_c_stat > temp_az[[1]])
  p_val_az_b <- mean(az_b_stat > temp_az[[2]])

  return(c(p_val_rr_c, p_val_rr_b,
           quantile(rr_c_stat, 1 - alpha),
           quantile(rr_b_stat, 1 - alpha),
           p_val_vv_c, p_val_vv_b,
           quantile(vv_c_stat, 1 - alpha),
           quantile(vv_b_stat, 1 - alpha),
           p_val_vn_c, p_val_vn_b,
           quantile(vn_c_stat, 1 - alpha),
           quantile(vn_b_stat, 1 - alpha),
           p_val_az_c, p_val_az_b,
           quantile(az_c_stat, 1 - alpha),
           quantile(az_b_stat, 1 - alpha)))
}
