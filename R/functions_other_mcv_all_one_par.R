
# function vec from ks package
vec_from_ks_package <- function (x, byrow = FALSE) {
  if (is.vector(x))
    return(x)
  if (byrow) 
    x <- t(x)
  d <- ncol(x)
  vecx <- vector()
  for (j in 1:d) vecx <- c(vecx, x[, j])
  return(vecx)
}

# all mcv in one functions
A_mcv_mu_sigma <- function(mu, sigma) {
  D_tilde_rcpp_mu <- D_tilde_rcpp(mu)
  tmu_mu <- ((t(mu) %*% mu)[1, 1])
  # rr
  dd <- length(mu)
  temp_1 <- -2 * det(sigma) * t(mu) * dd / tmu_mu^(dd + 1)
  temp_2 <- matrix(det(sigma) * vec_from_ks_package(solve(sigma)) / tmu_mu^dd,
                   nrow = 1)
  a_rr <- cbind(temp_1 + temp_2 %*% D_tilde_rcpp_mu, temp_2)
  # vv
  temp_1 <- -2 * sum(diag(sigma)) * t(mu) / tmu_mu^2
  temp_2 <- matrix(vec_from_ks_package(diag(nrow(sigma))) / tmu_mu, nrow = 1)
  a_vv <- cbind(temp_1 + temp_2 %*% D_tilde_rcpp_mu, temp_2)
  # vn
  temp_1 <- 2 * t(mu) %*% solve(sigma)
  temp_2 <- -kronecker(temp_1 / 2, temp_1 / 2)
  a_vn <- cbind(temp_1 + temp_2 %*% D_tilde_rcpp_mu, temp_2)
  # az
  temp_1 <- (-2 * t(mu) %*% sigma %*% (2 * mu %*% t(mu) / tmu_mu - 
                                         diag(nrow(sigma)))) / tmu_mu^2
  temp_2 <- kronecker(t(mu), t(mu)) / tmu_mu^2
  a_az <- cbind(temp_1 + temp_2 %*% D_tilde_rcpp_mu, temp_2)
  
  return(list(a_rr = a_rr, a_vv = a_vv, a_vn = a_vn, a_az = a_az))
}

# sigma for confidence intervals of mcvs and their inverses
sigma2_est_ci <- function(x) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  A_est <- A_mcv_mu_sigma(mu_est, sigma_est)
  Psi_est_mat <- Psi_est_rcpp(x)
  weight_mat_est <- rbind(cbind(sigma_est, t(Psi_est_mat[[1]])), 
                          cbind(Psi_est_mat[[1]], Psi_est_mat[[2]]))
  mcv_est_temp <- mcv_est_sigma2(mu_est, sigma_est, dd)
  
  temp_var_rr <- (1 / (4 * (dd^2))) * A_est$a_rr %*% weight_mat_est %*% t(A_est$a_rr)
  sigma_2_rr_C_est <- (mcv_est_temp$rr^(2 * dd))^(1 / dd - 2) * temp_var_rr
  sigma_2_rr_B_est <- (mcv_est_temp$rr^(2 * dd))^(-1 / dd - 2) * temp_var_rr
  
  temp_var_vv <- (1 / 4) * A_est$a_vv %*% weight_mat_est %*% t(A_est$a_vv)
  sigma_2_vv_C_est <- (1 / mcv_est_temp$vv^2) * temp_var_vv
  sigma_2_vv_B_est <- ((mcv_est_temp$vv^2)^(-3)) * temp_var_vv
  
  temp_var_vn <- (1 / 4) * A_est$a_vn %*% weight_mat_est %*% t(A_est$a_vn)
  sigma_2_vn_C_est <- ((mcv_est_temp$vn^2)^3) * temp_var_vn
  sigma_2_vn_B_est <- ((mcv_est_temp$vn)^2) * temp_var_vn
  
  temp_var_az <- (1 / 4) * A_est$a_az %*% weight_mat_est %*% t(A_est$a_az)
  sigma_2_az_C_est <- (1 / mcv_est_temp$az^2) * temp_var_az
  sigma_2_az_B_est <- ((mcv_est_temp$az^2)^(-3)) * temp_var_az
  
  return(list(rr_c = sigma_2_rr_C_est[1, 1],
              rr_b = sigma_2_rr_B_est[1, 1],
              vv_c = sigma_2_vv_C_est[1, 1],
              vv_b = sigma_2_vv_B_est[1, 1],
              vn_c = sigma_2_vn_C_est[1, 1],
              vn_b = sigma_2_vn_B_est[1, 1],
              az_c = sigma_2_az_C_est[1, 1],
              az_b = sigma_2_az_B_est[1, 1]))
}

#' Estimators and confidence intervals of four multivariate coefficients of variation and their reciprocals
#' 
#' Calculates the estimators with respective \eqn{(1-\alpha)}-confidence intervals for the four different variants of the multivariate coefficients (MCV) and their reciprocals
#' by Reyment (1960), Van Valen (1974), Voinov and Nikulin (1996) and Albert and Zhang (2010). 
#'    
#' @param x a matrix of data of size \eqn{n\times d}.
#' @param conf_level a confidence level. By default, it is equal to 0.95.
#' 
#' @details The function \code{e_mcv()} calculates four different variants of multivariate coefficient of variation for \eqn{d}-dimensional data. These variant were introduced by
#' by Reyment (1960, RR), Van Valen (1974, VV), Voinov and Nikulin (1996, VN) and Albert and Zhang (2010, AZ):
#'     \deqn{
#'     {\widehat C}^{RR}=\sqrt{\frac{(\det\mathbf{\widehat\Sigma})^{1/d}}{\boldsymbol{\widehat\mu}^{\top}\boldsymbol{\widehat\mu}}},\ 
#'     {\widehat C}^{VV}=\sqrt{\frac{\mathrm{tr}\mathbf{\widehat\Sigma}}{\boldsymbol{\widehat\mu}^{\top}\boldsymbol{\widehat\mu}}},\ 
#'     {\widehat C}^{VN}=\sqrt{\frac{1}{\boldsymbol{\widehat\mu}^{\top}\mathbf{\widehat\Sigma}^{-1}\boldsymbol{\widehat\mu}}},\ 
#'     {\widehat C}^{AZ}=\sqrt{\frac{\boldsymbol{\widehat\mu}^{\top}\mathbf{\widehat\Sigma}\boldsymbol{\widehat\mu}}{(\boldsymbol{\widehat\mu}^{\top}\boldsymbol{\widehat\mu})^2}},
#'      }
#'  where \eqn{n} is the sample size, \eqn{\boldsymbol{\widehat\mu}} is the empirical mean vector and \eqn{\mathbf{\widehat \Sigma}} is the empirical covariance matrix:
#'    \deqn{
#'      \boldsymbol{\widehat\mu}_i = \frac{1}{n}\sum_{j=1}^{n} \mathbf{X}_{j},\; \mathbf{\widehat \Sigma} =\frac{1}{n}\sum_{j=1}^{n} (\mathbf{X}_{j} - \boldsymbol{\widehat \mu})(\mathbf{X}_{j} - \boldsymbol{\widehat \mu})^{\top}.
#'    }
#' In the univariate case (\eqn{d=1}), all four variants reduce to coefficient of variation. Furthermore, their reciprocals, the so-called standardized means, are determined:     
#'    \deqn{
#'      {\widehat B}^{RR}=\sqrt{\frac{\boldsymbol{\widehat\mu}^{\top}\boldsymbol{\widehat\mu}}{(\det\mathbf{\widehat\Sigma})^{1/d}}},\ 
#'      {\widehat B}^{VV}=\sqrt{\frac{\boldsymbol{\widehat\mu}^{\top}\boldsymbol{\widehat\mu}}{\mathrm{tr}\mathbf{\widehat\Sigma}}},\ 
#'      {\widehat B}^{VN}=\sqrt{\boldsymbol{\widehat\mu}^{\top}\mathbf{\widehat\Sigma}^{-1}\boldsymbol{\widehat\mu}},\ 
#'      {\widehat B}^{AZ}=\sqrt{\frac{(\boldsymbol{\widehat\mu}^{\top}\boldsymbol{\widehat\mu})^2}{\boldsymbol{\widehat\mu}^{\top}\mathbf{\widehat\Sigma}\boldsymbol{\widehat\mu}}}.
#'    }
#'  In addition to the estimators, the respective confidence intervals [\code{C_lwr}, \code{C_upr}] for a given confidence level \eqn{1-\alpha} are calculated by the \code{e_mcv()} function. 
#'  These confidence intervals are based on an asymptotic approximation by a normal distribution, see Ditzhaus and Smaga (2023) for the technical details. These approximations
#' do not rely on any specific (semi-)parametric assumption on the distribution and are valid nonparametrically, even for tied data. 

#' @return When \eqn{d>1} (respectively \eqn{d=1}) a data frame with four rows (one row) corresponding to the four  MCVs (the univariate CV)
#'  and six columns containing the estimators \code{C_est} for the MCV (CV) and the estimators \code{B_est} for their reciprocals  as well as the upper and lower bounds of the corresponding
#'  confidence intervals  [\code{C_lwr}, \code{C_upr}] and [\code{B_lwr}, \code{B_upr}].
#'  
#' @references Albert A., Zhang L. (2010) A novel definition of the multivariate coefficient of variation. Biometrical Journal 52:667-675.
#' 
#' Ditzhaus M., Smaga L. (2023) Inference for all variants of the multivariate coefficient of variation in factorial designs. Preprint https://arxiv.org/abs/2301.12009.
#'  
#' Reyment R.A. (1960) Studies on Nigerian Upper Cretaceous and Lower Tertiary Ostracoda: part 1. Senonian and Maastrichtian Ostracoda, Stockholm Contributions in Geology, vol 7.
#'  
#' Van Valen L. (1974) Multivariate structural statistics in natural history. Journal of Theoretical Biology 45:235-247.
#'  
#' Voinov V., Nikulin M. (1996) Unbiased Estimators and Their Applications, Vol. 2, Multivariate Case. Kluwer, Dordrecht.
#'  
#' @examples 
#' # d > 1 (MCVs)
#' data_set <- lapply(list(iris[iris$Species == "setosa", 1:3],
#'                         iris[iris$Species == "versicolor", 1:3],
#'                         iris[iris$Species == "virginica", 1:3]),
#'                    as.matrix)
#' lapply(data_set, e_mcv)
#' # d = 1 (CV)
#' data_set <- lapply(list(iris[iris$Species == "setosa", 1],
#'                         iris[iris$Species == "versicolor", 1],
#'                         iris[iris$Species == "virginica", 1]),
#'                    as.matrix)
#' lapply(data_set, e_mcv)
#'
#' @importFrom stats cov qnorm
#'  
#' @export
e_mcv <- function(x, conf_level = 0.95) {
  alpha <- 1 - conf_level
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  tmu_mu <- ((t(mu_est) %*% mu_est)[1, 1])
  
  RR <- sqrt((det(sigma_est)^(1 / dd)) / tmu_mu)
  VV <- sqrt(sum(diag(sigma_est)) / tmu_mu)
  VN <- (sqrt(1 / (t(mu_est) %*% solve(sigma_est) %*% mu_est)))[1, 1]
  AZ <- (sqrt(t(mu_est) %*% sigma_est %*% mu_est / tmu_mu^2))[1, 1]
  
  zz <- stats::qnorm(1 - alpha / 2)
  RR_C_l <- RR[1] - (1 / sqrt(n_i)) * sigma2_est_ci(x)[[1]] * zz
  RR_C_u <- RR[1] + (1 / sqrt(n_i)) * sigma2_est_ci(x)[[1]] * zz
  RR_B_l <- 1 / RR[1] - (1 / sqrt(n_i)) * sigma2_est_ci(x)[[2]] * zz
  RR_B_u <- 1 / RR[1] + (1 / sqrt(n_i)) * sigma2_est_ci(x)[[2]] * zz
  VV_C_l <- VV[1] - (1 / sqrt(n_i)) * sigma2_est_ci(x)[[3]] * zz
  VV_C_u <- VV[1] + (1 / sqrt(n_i)) * sigma2_est_ci(x)[[3]] * zz
  VV_B_l <- 1 / VV[1] - (1 / sqrt(n_i)) * sigma2_est_ci(x)[[4]] * zz
  VV_B_u <- 1 / VV[1] + (1 / sqrt(n_i)) * sigma2_est_ci(x)[[4]] * zz
  VN_C_l <- VN[1] - (1 / sqrt(n_i)) * sigma2_est_ci(x)[[5]] * zz
  VN_C_u <- VN[1] + (1 / sqrt(n_i)) * sigma2_est_ci(x)[[5]] * zz
  VN_B_l <- 1 / VN[1] - (1 / sqrt(n_i)) * sigma2_est_ci(x)[[6]] * zz
  VN_B_u <- 1 / VN[1] + (1 / sqrt(n_i)) * sigma2_est_ci(x)[[6]] * zz
  AZ_C_l <- AZ[1] - (1 / sqrt(n_i)) * sigma2_est_ci(x)[[7]] * zz
  AZ_C_u <- AZ[1] + (1 / sqrt(n_i)) * sigma2_est_ci(x)[[7]] * zz
  AZ_B_l <- 1 / AZ[1] - (1 / sqrt(n_i)) * sigma2_est_ci(x)[[8]] * zz
  AZ_B_u <- 1 / AZ[1] + (1 / sqrt(n_i)) * sigma2_est_ci(x)[[8]] * zz
  
  # res <- list(RR_C = RR, RR_C_l = RR_C_l, RR_C_u = RR_C_u,
  #             RR_B = 1 / RR[1], RR_B_l = RR_B_l, RR_B_u = RR_B_u,
  #             VV_C = VV, VV_C_l = VV_C_l, VV_C_u = VV_C_u,
  #             VV_B = 1 / VV[1], VV_B_l = VV_B_l, VV_B_u = VV_B_u,
  #             VN_C = VN, VN_C_l = VN_C_l, VN_C_u = VN_C_u,
  #             VN_B = 1 / VN[1], VN_B_l = VN_B_l, VN_B_u = VN_B_u,
  #             AZ_C = AZ, AZ_C_l = AZ_C_l, AZ_C_u = AZ_C_u,
  #             AZ_B = 1 / AZ[1], AZ_B_l = AZ_B_l, AZ_B_u = AZ_B_u)
  
  if (dd > 1) {
    res <- data.frame(row.names = c("RR", "VV", "VN", "AZ"),
                      C_est =  c(RR, VV, VN, AZ),
                      C_lwr =  c(RR_C_l, VV_C_l, VN_C_l, AZ_C_l),
                      C_upr =  c(RR_C_u, VV_C_u, VN_C_u, AZ_C_u),
                      B_est =  c(1 / RR[1], 1 / VV[1], 1 / VN[1], 1 / AZ[1]),
                      B_lwr =  c(RR_B_l, VV_B_l, VN_B_l, AZ_B_l),
                      B_upr =  c(RR_B_u, VV_B_u, VN_B_u, AZ_B_u))
  } else if (dd == 1) {
    res <- data.frame(row.names = c("CV"),
                      C_est =  c(RR),
                      C_lwr =  c(RR_C_l),
                      C_upr =  c(RR_C_u),
                      B_est =  c(1 / RR[1]),
                      B_lwr =  c(RR_B_l),
                      B_upr =  c(RR_B_u))
  }
  
  return(res)
}
# @return A list containing the following components:
#  \item{RR_C}{a value of the estimator of MCV of Reyment (1960).}
#  \item{RR_C_l}{a value of the lower limit of confidence interval of MCV of Reyment (1960).}
#  \item{RR_C_u}{a value of the upper limit of confidence interval of MCV of Reyment (1960).}
#  \item{RR_B}{a value of the estimator of reciprocal of MCV of Reyment (1960).}
#  \item{RR_B_l}{a value of the lower limit of confidence interval of reciprocal of MCV of Reyment (1960).}
#  \item{RR_B_u}{a value of the upper limit of confidence interval of reciprocal of MCV of Reyment (1960).}
#  \item{VV_C}{a value of the estimator of MCV of Van Valen (1974).}
#  \item{VV_C_l}{a value of the lower limit of confidence interval of MCV of Van Valen (1974).}
#  \item{VV_C_u}{a value of the upper limit of confidence interval of MCV of Van Valen (1974).}
#  \item{VV_B}{a value of the estimator of reciprocal of MCV of Van Valen (1974).}
#  \item{VV_B_l}{a value of the lower limit of confidence interval of reciprocal of MCV of Van Valen (1974).}
#  \item{VV_B_u}{a value of the upper limit of confidence interval of reciprocal of MCV of Van Valen (1974).}
#  \item{VN_C}{a value of the estimator of MCV of Voinov and Nikulin (1996).}
#  \item{VN_C_l}{a value of the lower limit of confidence interval of MCV of Voinov and Nikulin (1996).}
#  \item{VN_C_u}{a value of the upper limit of confidence interval of MCV of Voinov and Nikulin (1996).}
#  \item{VN_B}{a value of the estimator of reciprocal of MCV of Voinov and Nikulin (1996).}
#  \item{VN_B_l}{a value of the lower limit of confidence interval of reciprocal of MCV of Voinov and Nikulin (1996).}
#  \item{VN_B_u}{a value of the upper limit of confidence interval of reciprocal of MCV of Voinov and Nikulin (1996).}
#  \item{AZ_C}{a value of the estimator of MCV of Albert and Zhang (2010).}
#  \item{AZ_C_l}{a value of the lower limit of confidence interval of MCV of Albert and Zhang (2010).}
#  \item{AZ_C_u}{a value of the upper limit of confidence interval of MCV of Albert and Zhang (2010).}
#  \item{AZ_B}{a value of the estimator of reciprocal of MCV of Albert and Zhang (2010).}
#  \item{AZ_B_l}{a value of the lower limit of confidence interval of reciprocal of MCV of Albert and Zhang (2010).}
#  \item{AZ_B_u}{a value of the upper limit of confidence interval of reciprocal of MCV of Albert and Zhang (2010).}

# mcv estimators for test procedures
mcv_est <- function(x) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  tmu_mu <- ((t(mu_est) %*% mu_est)[1, 1])
  return(list(RR = sqrt((det(sigma_est)^(1 / dd)) / tmu_mu),
              VV = sqrt(sum(diag(sigma_est)) / tmu_mu),
              VN = (sqrt(1 / (t(mu_est) %*% solve(sigma_est) %*% mu_est)))[1, 1],
              AZ = (sqrt(t(mu_est) %*% sigma_est %*% mu_est / tmu_mu^2))[1, 1]))
}

# for sigma2_est
mcv_est_sigma2 <- function(mu_est, sigma_est, dd) {
  tmu_mu <- ((t(mu_est) %*% mu_est)[1, 1])
  return(list(rr = sqrt((det(sigma_est)^(1 / dd)) / tmu_mu),
              vv = sqrt(sum(diag(sigma_est)) / tmu_mu),
              vn = (sqrt(1 / (t(mu_est) %*% solve(sigma_est) %*% mu_est)))[1, 1],
              az = (sqrt(t(mu_est) %*% sigma_est %*% mu_est / tmu_mu^2))[1, 1]))
}

# sigma estimators for test procedures
sigma2_est <- function(x, nn) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  A_est <- A_mcv_mu_sigma(mu_est, sigma_est)
  Psi_est_mat <- Psi_est_rcpp(x)
  weight_mat_est <- rbind(cbind(sigma_est, t(Psi_est_mat[[1]])), 
                          cbind(Psi_est_mat[[1]], Psi_est_mat[[2]]))
  mcv_est_temp <- mcv_est_sigma2(mu_est, sigma_est, dd)
  
  temp_var_rr <- (nn / (4 * (dd^2) * n_i)) * A_est$a_rr %*% weight_mat_est %*% t(A_est$a_rr)
  sigma_2_rr_C_est <- (mcv_est_temp$rr^(2 * dd))^(1 / dd - 2) * temp_var_rr
  sigma_2_rr_B_est <- (mcv_est_temp$rr^(2 * dd))^(-1 / dd - 2) * temp_var_rr
  
  temp_var_vv <- (nn / (4 * n_i)) * A_est$a_vv %*% weight_mat_est %*% t(A_est$a_vv)
  sigma_2_vv_C_est <- (1 / mcv_est_temp$vv^2) * temp_var_vv
  sigma_2_vv_B_est <- ((mcv_est_temp$vv^2)^(-3)) * temp_var_vv
  
  temp_var_vn <- (nn / (4 * n_i)) * A_est$a_vn %*% weight_mat_est %*% t(A_est$a_vn)
  sigma_2_vn_C_est <- ((mcv_est_temp$vn^2)^3) * temp_var_vn
  sigma_2_vn_B_est <- ((mcv_est_temp$vn)^2) * temp_var_vn
  
  temp_var_az <- (nn / (4 * n_i)) * A_est$a_az %*% weight_mat_est %*% t(A_est$a_az)
  sigma_2_az_C_est <- (1 / mcv_est_temp$az^2) * temp_var_az
  sigma_2_az_B_est <- ((mcv_est_temp$az^2)^(-3)) * temp_var_az
  
  return(list(rr_c = sigma_2_rr_C_est[1, 1],
              rr_b = sigma_2_rr_B_est[1, 1],
              vv_c = sigma_2_vv_C_est[1, 1],
              vv_b = sigma_2_vv_B_est[1, 1],
              vn_c = sigma_2_vn_C_est[1, 1],
              vn_b = sigma_2_vn_B_est[1, 1],
              az_c = sigma_2_az_C_est[1, 1],
              az_b = sigma_2_az_B_est[1, 1]))
}

# test statistics W_nC(T) and W_nB(T)
# x - a list of length k with elements being n_i\times d matrices of data
wald_test_stat <- function(x, TT = NULL) {
  kk <- length(x)
  nn <- sum(sapply(x, nrow))
  if (is.null(TT)) TT <- diag(rep(1, kk)) - (1 / kk) * matrix(1, ncol = kk, nrow = kk)
  
  mcv_c_est_temp <- sapply(x, mcv_est)
  mcv_b_est_temp <- sapply(x, function(y) lapply(mcv_est(y), function(z) 1 / z))
  
  sigma2_cb_est_mat <- sapply(x, sigma2_est, nn = nn)
  
  sigma_rr_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_rr_b_est <- diag(sigma2_cb_est_mat[2, ])
  
  sigma_vv_c_est <- diag(sigma2_cb_est_mat[3, ])
  sigma_vv_b_est <- diag(sigma2_cb_est_mat[4, ])
  
  sigma_vn_c_est <- diag(sigma2_cb_est_mat[5, ])
  sigma_vn_b_est <- diag(sigma2_cb_est_mat[6, ])
  
  sigma_az_c_est <- diag(sigma2_cb_est_mat[7, ])
  sigma_az_b_est <- diag(sigma2_cb_est_mat[8, ])
  
  
  rr_c_stat <- nn * (t(TT %*% diag(diag(mcv_c_est_temp[1, ]))) %*% 
                       MASS::ginv(TT %*% sigma_rr_c_est %*% t(TT)) %*% 
                       TT %*% diag(diag(mcv_c_est_temp[1, ])))[1, 1]
  
  rr_b_stat <- nn * (t(TT %*% diag(diag(mcv_b_est_temp[1, ]))) %*% 
                       MASS::ginv(TT %*% sigma_rr_b_est %*% t(TT)) %*% 
                       TT %*% diag(diag(mcv_b_est_temp[1, ])))[1, 1]
  
  vv_c_stat <- nn * (t(TT %*% diag(diag(mcv_c_est_temp[2, ]))) %*% 
                       MASS::ginv(TT %*% sigma_vv_c_est %*% t(TT)) %*% 
                       TT %*% diag(diag(mcv_c_est_temp[2, ])))[1, 1]
  
  vv_b_stat <- nn * (t(TT %*% diag(diag(mcv_b_est_temp[2, ]))) %*% 
                       MASS::ginv(TT %*% sigma_vv_b_est %*% t(TT)) %*% 
                       TT %*% diag(diag(mcv_b_est_temp[2, ])))[1, 1]
  
  vn_c_stat <- nn * (t(TT %*% diag(diag(mcv_c_est_temp[3, ]))) %*% 
                       MASS::ginv(TT %*% sigma_vn_c_est %*% t(TT)) %*% 
                       TT %*% diag(diag(mcv_c_est_temp[3, ])))[1, 1]
  
  vn_b_stat <- nn * (t(TT %*% diag(diag(mcv_b_est_temp[3, ]))) %*% 
                       MASS::ginv(TT %*% sigma_vn_b_est %*% t(TT)) %*% 
                       TT %*% diag(diag(mcv_b_est_temp[3, ])))[1, 1]
  
  az_c_stat <- nn * (t(TT %*% diag(diag(mcv_c_est_temp[4, ]))) %*% 
                       MASS::ginv(TT %*% sigma_az_c_est %*% t(TT)) %*% 
                       TT %*% diag(diag(mcv_c_est_temp[4, ])))[1, 1]
  
  az_b_stat <- nn * (t(TT %*% diag(diag(mcv_b_est_temp[4, ]))) %*% 
                       MASS::ginv(TT %*% sigma_az_b_est %*% t(TT)) %*% 
                       TT %*% diag(diag(mcv_b_est_temp[4, ])))[1, 1]
  
  return(c(rr_c_stat, rr_b_stat,
           vv_c_stat, vv_b_stat,
           vn_c_stat, vn_b_stat,
           az_c_stat, az_b_stat))
}

perm_test_wald <- function(x, TT = NULL, n_perm = 1000,
                           parallel = FALSE, n_cores = NULL) {
  pp <- ncol(x[[1]])
  if (pp == 1) {
    # stop("the dimension is equal to 1")
  }
  kk <- length(x)
  if (is.null(TT)) TT <- diag(rep(1, kk)) - (1 / kk) * matrix(1, ncol = kk, nrow = kk)
  n_i_vec <- sapply(x, nrow)
  nn <- sum(n_i_vec)
  x_mat <- x[[1]]
  for (i_kk in 2:kk) {
    x_mat <- rbind(x_mat, x[[i_kk]])
  }
  groups <- rep(1:kk, n_i_vec)
  if (!parallel) {
    test_stat_perm <- matrix(0, nrow = n_perm, ncol = 8)
    for (i_perm in seq_len(n_perm)) {
      if (pp > 1) {
        x_mat_perm <- x_mat[sample(nn), ]
      } else if (pp == 1) {
        x_mat_perm <- matrix(x_mat[sample(nn), ], ncol = 1)
      }
      x_perm_list <- vector("list", kk)
      for (j_kk in 1:kk) {
        if (pp > 1) {
          x_perm_list[[j_kk]] <- x_mat_perm[groups == j_kk, ]
        } else if (pp == 1) {
          x_perm_list[[j_kk]] <- matrix(x_mat_perm[groups == j_kk, ], ncol = 1)
        }
      }
      test_stat_perm[i_perm, ] <- wald_test_stat(x_perm_list, TT = TT)
    }
  } else {
    RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion")
    test_stat_perm <- foreach(i_perm = 1:n_perm, .combine = rbind, 
                              .packages = c("Rcpp", "MASS", "mvtnorm", "GFDmcv")) %dopar%
      {
        if (pp > 1) {
          x_mat_perm <- x_mat[sample(nn), ]
        } else if (pp == 1) {
          x_mat_perm <- matrix(x_mat[sample(nn), ], ncol = 1)
        }
        x_perm_list <- vector("list", kk)
        for (j_kk in 1:kk) {
          if (pp > 1) {
            x_perm_list[[j_kk]] <- x_mat_perm[groups == j_kk, ]
          } else if (pp == 1) {
            x_perm_list[[j_kk]] <- matrix(x_mat_perm[groups == j_kk, ], ncol = 1)
          }
        }
        wald_test_stat(x_perm_list, TT = TT)
      }
  }
  return(colMeans(test_stat_perm > matrix(wald_test_stat(x, TT = TT), 
                                          nrow = n_perm, ncol = 8, byrow = TRUE)))
}

boot_test_wald <- function(x, TT = NULL, n_boot = 1000,
                           parallel = FALSE, n_cores = NULL) {
  pp <- ncol(x[[1]])
  if (pp == 1) {
    # stop("the dimension is equal to 1")
  }
  kk <- length(x)
  if (is.null(TT)) TT <- diag(rep(1, kk)) - (1 / kk) * matrix(1, ncol = kk, nrow = kk)
  n_i_vec <- sapply(x, nrow)
  nn <- sum(n_i_vec)
  x_mat <- x[[1]]
  for (i_kk in 2:kk) {
    x_mat <- rbind(x_mat, x[[i_kk]])
  }
  groups <- rep(1:kk, n_i_vec)
  if (!parallel) {
    test_stat_boot <- matrix(0, nrow = n_boot, ncol = 8)
    for (i_boot in seq_len(n_boot)) {
      if (pp > 1) {
        x_mat_boot <- x_mat[sample(nn, replace = TRUE), ]
      } else if (pp == 1) {
        x_mat_boot <- matrix(x_mat[sample(nn, replace = TRUE), ], ncol = 1)
      }
      x_boot_list <- vector("list", kk)
      for (j_kk in 1:kk) {
        if (pp > 1) {
          x_boot_list[[j_kk]] <- x_mat_boot[groups == j_kk, ]
        } else if (pp == 1) {
          x_boot_list[[j_kk]] <- matrix(x_mat_boot[groups == j_kk, ], ncol = 1)
        }
      }
      test_stat_boot[i_boot, ] <- wald_test_stat(x_boot_list, TT = TT)
    }
  } else {
    RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion")
    test_stat_boot <- foreach(i_n_boot = 1:n_boot, .combine = rbind, 
                              .packages = c("Rcpp", "MASS", "mvtnorm", "GFDmcv")) %dopar%
      {
        if (pp > 1) {
          x_mat_boot <- x_mat[sample(nn, replace = TRUE), ]
        } else if (pp == 1) {
          x_mat_boot <- matrix(x_mat[sample(nn, replace = TRUE), ], ncol = 1)
        }
        x_boot_list <- vector("list", kk)
        for (j_kk in 1:kk) {
          if (pp > 1) {
            x_boot_list[[j_kk]] <- x_mat_boot[groups == j_kk, ]
          } else if (pp == 1) {
            x_boot_list[[j_kk]] <- matrix(x_mat_boot[groups == j_kk, ], ncol = 1)
          }
        }
        wald_test_stat(x_boot_list, TT = TT)
      }
  }
  return(colMeans(test_stat_boot > matrix(wald_test_stat(x, TT = TT), 
                                          nrow = n_boot, ncol = 8, byrow = TRUE)))
}

# separate for mct tests

# rr
A_mcv_mu_sigma_rr <- function(mu, sigma) {
  D_tilde_rcpp_mu <- D_tilde_rcpp(mu)
  tmu_mu <- ((t(mu) %*% mu)[1, 1])
  # rr
  dd <- length(mu)
  temp_1 <- -2 * det(sigma) * t(mu) * dd / tmu_mu^(dd + 1)
  temp_2 <- matrix(det(sigma) * vec_from_ks_package(solve(sigma)) / tmu_mu^dd,
                   nrow = 1)
  a_rr <- cbind(temp_1 + temp_2 %*% D_tilde_rcpp_mu, temp_2)
  return(list(a_rr = a_rr))
}

mcv_est_rr <- function(x) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  tmu_mu <- ((t(mu_est) %*% mu_est)[1, 1])
  return(list(rr = sqrt((det(sigma_est)^(1 / dd)) / tmu_mu)))
}

# for sigma2_est_rr
mcv_est_sigma2_rr <- function(mu_est, sigma_est, dd) {
  tmu_mu <- ((t(mu_est) %*% mu_est)[1, 1])
  return(list(rr = sqrt((det(sigma_est)^(1 / dd)) / tmu_mu)))
}

sigma2_est_rr <- function(x, nn) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  A_est <- A_mcv_mu_sigma_rr(mu_est, sigma_est)
  Psi_est_mat <- Psi_est_rcpp(x)
  weight_mat_est <- rbind(cbind(sigma_est, t(Psi_est_mat[[1]])), 
                          cbind(Psi_est_mat[[1]], Psi_est_mat[[2]]))
  mcv_est_temp <- mcv_est_sigma2_rr(mu_est, sigma_est, dd)
  
  temp_var_rr <- (nn / (4 * (dd^2) * n_i)) * A_est$a_rr %*% weight_mat_est %*% t(A_est$a_rr)
  sigma_2_rr_C_est <- (mcv_est_temp$rr^(2 * dd))^(1 / dd - 2) * temp_var_rr
  sigma_2_rr_B_est <- (mcv_est_temp$rr^(2 * dd))^(-1 / dd - 2) * temp_var_rr
  
  return(list(rr_c = sigma_2_rr_C_est[1, 1],
              rr_b = sigma_2_rr_B_est[1, 1]))
}

# vv
A_mcv_mu_sigma_vv <- function(mu, sigma) {
  D_tilde_rcpp_mu <- D_tilde_rcpp(mu)
  tmu_mu <- ((t(mu) %*% mu)[1, 1])
  # vv
  temp_1 <- -2 * sum(diag(sigma)) * t(mu) / tmu_mu^2
  temp_2 <- matrix(vec_from_ks_package(diag(nrow(sigma))) / tmu_mu, nrow = 1)
  a_vv <- cbind(temp_1 + temp_2 %*% D_tilde_rcpp_mu, temp_2)
  
  return(list(a_vv = a_vv))
}

mcv_est_vv <- function(x) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  tmu_mu <- ((t(mu_est) %*% mu_est)[1, 1])
  return(list(vv = sqrt(sum(diag(sigma_est)) / tmu_mu)))
}

# for sigma2_est_vv
mcv_est_sigma2_vv <- function(mu_est, sigma_est, dd) {
  tmu_mu <- ((t(mu_est) %*% mu_est)[1, 1])
  return(list(vv = sqrt(sum(diag(sigma_est)) / tmu_mu)))
}

sigma2_est_vv <- function(x, nn) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  A_est <- A_mcv_mu_sigma_vv(mu_est, sigma_est)
  Psi_est_mat <- Psi_est_rcpp(x)
  weight_mat_est <- rbind(cbind(sigma_est, t(Psi_est_mat[[1]])), 
                          cbind(Psi_est_mat[[1]], Psi_est_mat[[2]]))
  mcv_est_temp <- mcv_est_sigma2_vv(mu_est, sigma_est, dd)
  
  temp_var_vv <- (nn / (4 * n_i)) * A_est$a_vv %*% weight_mat_est %*% t(A_est$a_vv)
  sigma_2_vv_C_est <- (1 / mcv_est_temp$vv^2) * temp_var_vv
  sigma_2_vv_B_est <- ((mcv_est_temp$vv^2)^(-3)) * temp_var_vv
  
  return(list(vv_c = sigma_2_vv_C_est[1, 1],
              vv_b = sigma_2_vv_B_est[1, 1]))
}

# vn
A_mcv_mu_sigma_vn <- function(mu, sigma) {
  D_tilde_rcpp_mu <- D_tilde_rcpp(mu)
  # vn
  temp_1 <- 2 * t(mu) %*% solve(sigma)
  temp_2 <- -kronecker(temp_1 / 2, temp_1 / 2)
  a_vn <- cbind(temp_1 + temp_2 %*% D_tilde_rcpp_mu, temp_2)
  
  return(list(a_vn = a_vn))
}

mcv_est_vn <- function(x) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  return(list(vn = (sqrt(1 / (t(mu_est) %*% solve(sigma_est) %*% mu_est)))[1, 1]))
}

# for sigma2_est_vn
mcv_est_sigma2_vn <- function(mu_est, sigma_est, dd) {
  return(list(vn = sqrt(1 / (t(mu_est) %*% solve(sigma_est) %*% mu_est))))
}

sigma2_est_vn <- function(x, nn) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  A_est <- A_mcv_mu_sigma_vn(mu_est, sigma_est)
  Psi_est_mat <- Psi_est_rcpp(x)
  weight_mat_est <- rbind(cbind(sigma_est, t(Psi_est_mat[[1]])), 
                          cbind(Psi_est_mat[[1]], Psi_est_mat[[2]]))
  mcv_est_temp <- mcv_est_sigma2_vn(mu_est, sigma_est, dd)
  
  temp_var_vn <- (nn / (4 * n_i)) * A_est$a_vn %*% weight_mat_est %*% t(A_est$a_vn)
  sigma_2_vn_C_est <- ((mcv_est_temp$vn^2)^3) * temp_var_vn
  sigma_2_vn_B_est <- ((mcv_est_temp$vn)^2) * temp_var_vn
  
  return(list(vn_c = sigma_2_vn_C_est[1, 1],
              vn_b = sigma_2_vn_B_est[1, 1]))
}

# az
A_mcv_mu_sigma_az <- function(mu, sigma) {
  D_tilde_rcpp_mu <- D_tilde_rcpp(mu)
  tmu_mu <- ((t(mu) %*% mu)[1, 1])
  # az
  temp_1 <- (-2 * t(mu) %*% sigma %*% (2 * mu %*% t(mu) / tmu_mu - 
                                         diag(nrow(sigma)))) / tmu_mu^2
  temp_2 <- kronecker(t(mu), t(mu)) / tmu_mu^2
  a_az <- cbind(temp_1 + temp_2 %*% D_tilde_rcpp_mu, temp_2)
  
  return(list(a_az = a_az))
}

mcv_est_az <- function(x) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  tmu_mu <- ((t(mu_est) %*% mu_est)[1, 1])
  return(list(az = (sqrt(t(mu_est) %*% sigma_est %*% mu_est / tmu_mu^2))[1, 1]))
}

# for sigma2_est_az
mcv_est_sigma2_az <- function(mu_est, sigma_est, dd) {
  tmu_mu <- ((t(mu_est) %*% mu_est)[1, 1])
  return(list(az = sqrt(t(mu_est) %*% sigma_est %*% mu_est / tmu_mu^2)))
}

sigma2_est_az <- function(x, nn) {
  dd <- ncol(x)
  n_i <- nrow(x)
  mu_est <- colMeans(x)
  sigma_est <- ((n_i - 1) / n_i) * cov(x)
  A_est <- A_mcv_mu_sigma_az(mu_est, sigma_est)
  Psi_est_mat <- Psi_est_rcpp(x)
  weight_mat_est <- rbind(cbind(sigma_est, t(Psi_est_mat[[1]])), 
                          cbind(Psi_est_mat[[1]], Psi_est_mat[[2]]))
  mcv_est_temp <- mcv_est_sigma2_az(mu_est, sigma_est, dd)
  
  temp_var_az <- (nn / (4 * n_i)) * A_est$a_az %*% weight_mat_est %*% t(A_est$a_az)
  sigma_2_az_C_est <- (1 / mcv_est_temp$az^2) * temp_var_az
  sigma_2_az_B_est <- ((mcv_est_temp$az^2)^(-3)) * temp_var_az
  
  return(list(az_c = sigma_2_az_C_est[1, 1],
              az_b = sigma_2_az_B_est[1, 1]))
}
