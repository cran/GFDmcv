#' Inference for four multivariate coefficients of variation
#' 
#' The function \code{GFDmcv()} calculates the Wald-type statistic for global null hypotheses 
#' and max-type statistics for multiple local null hypotheses, both in terms of the four variants 
#' of the multivariate coefficient of variation. Respective \eqn{p}-values
#' are obtained by a \eqn{\chi^2}-approximation, a pooled bootstrap strategy and a pooled permutation approach (only for the
#' Wald-type statistic), respectively.
#' 
#' @param x a list of length \eqn{k} with elements being \eqn{n_i\times d} matrices of data, \eqn{i=1,\dots,k}.
#' @param h_mct a \eqn{r\times k} contrast matrix \eqn{\mathbf{H}} of full row rank for multiple contrast tests. 
#' Remember to specify it correctly taking into account the order of elements of the list \code{x}.
#' @param h_wald a \eqn{q\times k} contrast matrix \eqn{\mathbf{H}} of full row rank for the Wald-type tests. 
#' Remember to specify it correctly taking into account the order of elements of the list \code{x}. 
#' @param alpha a significance level (then \code{1-alpha} is the confidence level).
#' @param n_perm a number of permutation replicates.
#' @param n_boot a number of bootstrap replicates.
#' @param parallel a logical indicating whether to use parallelization.
#' @param n_cores if \code{parallel = TRUE}, a number of processes used in parallel computation. 
#' Its default value means that it will be equal to a number of cores of a computer used.
#' 
#' @details The function \code{GFDmcv()} calculates the Wald-type statistic for global null hypotheses of the form
#' \deqn{ \mathcal H_0: \mathbf{H} (C_1,\ldots,C_k)^\top = \mathbf{0}\,\,\text{and}\,\,\mathcal H_0: \mathbf{H} (B_1,\ldots,B_k)^\top = \mathbf{0},}
#' where \eqn{\mathbf{H}} is a contrast matrix reflecting the research question of interest and
#' \eqn{C_i} (\eqn{B_i}) are the subgroup-specific MCVs (and their reciprocal) by Reyment (1960, RR), Van Valen (1974, VV), 
#' Voinov and Nikulin (1996, VN) or Albert and Zhang (2010, AZ), respectively. 
#' We refer to the function \code{e_mcv()} for the detailed definitions of the different
#' variants. The \eqn{p}-value of the Wald-type statistic relies on a \eqn{\chi^2}-approximation,
#' a (pooled) bootstrap or permutation approach.
#' 
#' Furthermore, the function \code{GFDmcv()} calculates a max-type test statistic for the multiple
#' comparison of \eqn{q} local null hypotheses:
#'  \deqn{ \mathcal H_{0,\ell}: \mathbf{h_\ell}^\top \mathbf{C} = \mathbf{0}\,\,
#'  \text{or}\,\,\mathcal H_{0,\ell}: \mathbf{h_\ell}^\top \mathbf{B} = \mathbf{0}, \,\,\ell=1,\ldots,q,}
#' where \eqn{\mathbf{C}=(C_1,\ldots,C_k)^\top} and  \eqn{\mathbf{B}=(B_1,\ldots,B_k)^\top}. The \eqn{p}-values are determined by a Gaussian approximation and a bootstrap approach, respectively.
#' In addition to the local test decisions, multiple adjusted confidence intervals for the contrasts \eqn{\mathbf{h_{\ell}^{\top}\pmb{C}}} and
#' \eqn{\mathbf{h_{\ell}^{\top}\pmb{B}}}, respectively, are calculated.
#' 
#' Please have a look on the plot and summary functions designed for this package. They can be used
#' to simplify the output of \code{GFDmcv()}.
#' 
#' @return A list of class \code{gfdmcv} containing the following components:
#' \item{overall_res}{a list of two elements representing the results for testing
#' the global null hypothesis. The first one is a matrix \code{test_stat} of the
#' test statistics, while the second is a matrix \code{p_values} of the respective \eqn{p}-values.}
#' \item{mct_res}{all results of MCT tests for particular hypothesis in \code{h_mct}, i.e., 
#' the estimators and simultaneous confidence intervals for \eqn{\mathbf{h_{\ell}^{\top}\pmb{C}}} and 
#' for \eqn{\mathbf{h_{\ell}^{\top}\pmb{B}}}, the test statistics and critical values as well as the decisions.}
#' \item{h_mct}{an argument \code{h_mct}.}
#' \item{h_wald}{an argument \code{h_wald}.}
#' \item{alpha}{an argument \code{alpha}.}
#' 
#' @examples
#' # Some of the examples may run some time.
#' 
#' # one-way analysis for MCV and CV
#' # d > 1 (MCV)
#' data_set <- lapply(list(iris[iris$Species == "setosa", 1:3],
#'                         iris[iris$Species == "versicolor", 1:3],
#'                         iris[iris$Species == "virginica", 1:3]),
#'                    as.matrix)
#' # estimators and confidence intervals of MCVs and their reciprocals
#' lapply(data_set, e_mcv)
#' # contrast matrices
#' k <- length(data_set)
#' # Tukey's contrast matrix
#' h_mct <- contr_mat(k, type = "Tukey")
#' # centering matrix P_k
#' h_wald <- contr_mat(k, type = "center")
#' \donttest{
#' # testing without parallel computing
#' res <- GFDmcv(data_set, h_mct, h_wald)
#' summary(res, digits = 3)
#' oldpar <- par(mar = c(4, 5, 2, 0.3))
#' plot(res)
#' par(oldpar)}
#' \donttest{
#' # testing with parallel computing
#' library(doParallel)
#' res <- GFDmcv(data_set, h_mct, h_wald, parallel = TRUE, n_cores = 2)
#' summary(res, digits = 3)
#' oldpar <- par(mar = c(4, 5, 2, 0.3))
#' plot(res)
#' par(oldpar)}
#' \dontshow{
#' data_set_2 <- lapply(data_set, function(x) x[1:5, ])
#' res <- GFDmcv(data_set, h_mct, h_wald, n_perm = 2, n_boot = 2)
#' summary(res, digits = 3)
#' oldpar <- par(mar = c(4, 5, 2, 0.3))
#' plot(res)
#' par(oldpar)}
#' # two-way analysis for CV (based on the example in Ditzhaus and Smaga, 2022)
#' library(HSAUR)
#' data_set <- lapply(list(BtheB$bdi.pre[BtheB$drug == "No" & BtheB$length == "<6m"],
#'                         BtheB$bdi.pre[BtheB$drug == "No" & BtheB$length == ">6m"],
#'                         BtheB$bdi.pre[BtheB$drug == "Yes" & BtheB$length == "<6m"],
#'                         BtheB$bdi.pre[BtheB$drug == "Yes" & BtheB$length == ">6m"]), 
#'                    as.matrix)
#' # estimators and confidence intervals of CV and its reciprocal
#' lapply(data_set, e_mcv)
#' \donttest{
#' # interaction
#' h_mct <- contr_mat(4, type = "Tukey")
#' h_wald <- kronecker(contr_mat(2, type = "center"), 
#'                     contr_mat(2, type = "center"))
#' res <- GFDmcv(data_set, h_mct, h_wald)
#' summary(res, digits = 3)
#' oldpar <- par(mar = c(4, 6, 2, 0.1))
#' plot(res)
#' par(oldpar)}
#' \donttest{
#' # main effect drug
#' h_mct <- matrix(c(1, 1, -1, -1), nrow = 1)
#' h_wald <- kronecker(contr_mat(2, type = "center"), 0.5 * matrix(1, 1, 2))
#' res <- GFDmcv(data_set, h_mct, h_wald)
#' summary(res, digits = 3)
#' oldpar <- par(mar = c(4, 6, 2, 0.1))
#' plot(res)
#' par(oldpar)}
#' \donttest{
#' # main effect length
#' h_mct <- matrix(c(1, -1, 1, -1), nrow = 1)
#' h_wald <- kronecker(0.5 * matrix(1, 1, 2), contr_mat(2, type = "center"))
#' res <- GFDmcv(data_set, h_mct, h_wald)
#' summary(res, digits = 3)
#' oldpar <- par(mar = c(4, 6, 2, 0.1))
#' plot(res)
#' par(oldpar)}
#' 
#' @references Albert A., Zhang L. (2010) A novel definition of the multivariate coefficient of variation. Biometrical Journal 52:667-675.
#' 
#' Ditzhaus M., Smaga L. (2022) Permutation test for the multivariate coefficient of variation in factorial designs. 
#' Journal of Multivariate Analysis 187, 104848.
#' 
#' Ditzhaus M., Smaga L. (2023) Inference for all variants of the multivariate coefficient of variation in factorial designs. Preprint https://arxiv.org/abs/2301.12009.
#' 
#' Reyment R.A. (1960) Studies on Nigerian Upper Cretaceous and Lower Tertiary Ostracoda: part 1. 
#' Senonian and Maastrichtian Ostracoda, Stockholm Contributions in Geology, vol 7.
#' 
#' Van Valen L. (1974) Multivariate structural statistics in natural history. Journal of Theoretical Biology 45:235-247.
#' 
#' Voinov V., Nikulin M. (1996) Unbiased Estimators and Their Applications, Vol. 2, Multivariate Case. Kluwer, Dordrecht.
#' 
#' @importFrom stats cov pchisq p.adjust quantile
#' @importFrom utils installed.packages
#' @importFrom mvtnorm pmvnorm qmvnorm
#' @importFrom Matrix rankMatrix
#' @import doParallel
#' @import foreach
#' @import Rcpp
#' @import stringr
#' @import HSAUR
#' 
#' @export

# implementation of all tests for:
# global hypothesis (h_wald matrix for Wald tests, h_mct for MCT) - asym, perm, boot, mct_asym, mct_boot
# local hypothesis (rows of h_mct matrix) - mct_asym, mct_boot

# Let \eqn{r} be the number of contrasts in \code{h}.
# By default, it is equal to the centering matrix \eqn{\mathbf{P}_k = \mathbf{I}_k  - \mathbf{J}_k/k}, which 
# is the difference of the unity matrix \eqn{\mathbf{I}_k} and the scaled version of the matrix 
# \eqn{\mathbf{J}_k=\mathbf{1}\mathbf{1}^{\top} \in R^{k\times k}} consisting of \eqn{1}'s only. This matrix 
# leads to the null hypotheses of no group effect, i.e.,
# \deqn{\mathcal H_{0,C^v}: \{ \mathbf{P}_k \mathbf{C}^v = \mathbf{0} \} = \{C_1^v= \ldots = C_k^v\},}
# \eqn{v\in\{RR,VV,VN,AZ\}}.

# main argument:
# x - a list of length k with elements being n_i\times d matrices of data
# values:
# global_tests - p-values of global tests
# local_tests - particular contrast: 1 denotes rejection, 0 - nonrejection
GFDmcv <- function(x, h_mct, h_wald, alpha = 0.05, n_perm = 1000, n_boot = 1000,
                   parallel = FALSE, n_cores = NULL) {
  if (parallel) {
    if (!("doParallel" %in% rownames(installed.packages()))) {
      stop("Please install package 'doParallel'")
    }
    # require(foreach, quietly = TRUE)
    requireNamespace("foreach", quietly = TRUE)
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores()
    } else if (!(n_cores > 1)) {
      stop("n_cores should be greater than 1")
    }
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  }
  
  r <- nrow(h_mct)
  dd <- ncol(x[[1]])
  
  # global tests
  temp_mct <- c(mct_asympt_test_rr(x, h_mct, alpha),
                mct_asympt_test_vv(x, h_mct, alpha),
                mct_asympt_test_vn(x, h_mct, alpha),
                mct_asympt_test_az(x, h_mct, alpha),
                mct_pool_boot(x, h_mct, n_boot, alpha, parallel, n_cores))
  
  p_val_mct <- temp_mct[rep(c(TRUE, TRUE, FALSE, FALSE), 8)]
  temp_mct_q <- temp_mct[rep(c(FALSE, FALSE, TRUE, TRUE), 8)]
  
  p_values <- c(1 - pchisq(wald_test_stat(x, TT = h_wald), 
                           Matrix::rankMatrix(h_wald)),
                perm_test_wald(x, TT = h_wald, n_perm = n_perm,
                               parallel, n_cores),
                boot_test_wald(x, TT = h_wald, n_boot = n_boot,
                               parallel, n_cores),
                p_val_mct)
  
  results <- data.frame(asym_Wald_type_test = p_values[1:8],
                        perm_Wald_type_test = p_values[9:16],
                        boot_Wald_type_test = p_values[17:24],
                        asym_MCT_test = p_values[25:32],
                        boot_MCT_test = p_values[33:40])
  rownames(results) <- c("RR_C", "RR_B", "VV_C", "VV_B", "VN_C", "VN_B", "AZ_C", "AZ_B")
  overall_p_val <- t(results)
  
  # values of test statistics for global hypotheses
  overall_test_stat <- rbind(wald_test_stat(x, TT = h_wald),
                            c(mct_test_stat_rr(x, h_mct)[[1]],
                              mct_test_stat_rr(x, h_mct)[[2]],
                              mct_test_stat_vv(x, h_mct)[[1]],
                              mct_test_stat_vv(x, h_mct)[[2]],
                              mct_test_stat_vv(x, h_mct)[[1]],
                              mct_test_stat_vv(x, h_mct)[[2]],
                              mct_test_stat_az(x, h_mct)[[1]],
                              mct_test_stat_az(x, h_mct)[[2]]))
  rownames(overall_test_stat) <- c("Wald_type_test", "MCT_test")
  colnames(overall_test_stat) <- c("RR_C", "RR_B", "VV_C", "VV_B", "VN_C", "VN_B", "AZ_C", "AZ_B")
  
  # simultaneous confidence intervals
  nn <- sum(sapply(x, nrow))
  nn_m05 <- 1 / sqrt(nn)
  c_est_rr <- matrix(unlist(sapply(x, mcv_est_rr)), ncol = 1)
  b_est_rr <- matrix(unlist(sapply(x, function(y) lapply(mcv_est_rr(y), function(z) 1 / z))), ncol = 1)
  c_est_vv <- matrix(unlist(sapply(x, mcv_est_vv)), ncol = 1)
  b_est_vv <- matrix(unlist(sapply(x, function(y) lapply(mcv_est_vv(y), function(z) 1 / z))), ncol = 1)
  c_est_vn <- matrix(unlist(sapply(x, mcv_est_vn)), ncol = 1)
  b_est_vn <- matrix(unlist(sapply(x, function(y) lapply(mcv_est_vn(y), function(z) 1 / z))), ncol = 1)
  c_est_az <- matrix(unlist(sapply(x, mcv_est_az)), ncol = 1)
  b_est_az <- matrix(unlist(sapply(x, function(y) lapply(mcv_est_az(y), function(z) 1 / z))), ncol = 1)
  
  sigma2_cb_est_mat <- sapply(x, sigma2_est_rr, nn = nn)
  sigma_rr_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_rr_b_est <- diag(sigma2_cb_est_mat[2, ])
  sigma2_cb_est_mat <- sapply(x, sigma2_est_vv, nn = nn)
  sigma_vv_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_vv_b_est <- diag(sigma2_cb_est_mat[2, ])
  sigma2_cb_est_mat <- sapply(x, sigma2_est_vn, nn = nn)
  sigma_vn_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_vn_b_est <- diag(sigma2_cb_est_mat[2, ])
  sigma2_cb_est_mat <- sapply(x, sigma2_est_az, nn = nn)
  sigma_az_c_est <- diag(sigma2_cb_est_mat[1, ])
  sigma_az_b_est <- diag(sigma2_cb_est_mat[2, ])
  
  # MCT all results in a data.frame
  
  i_r <- 1
  h_i <- t(h_mct[i_r, ])
  C_stat <- c(mct_test_stat_sep_rr(x, h_mct)[[1]][i_r], "", mct_test_stat_sep_vv(x, h_mct)[[1]][i_r], "",
              mct_test_stat_sep_vn(x, h_mct)[[1]][i_r], "", mct_test_stat_sep_az(x, h_mct)[[1]][i_r], "")
  C_cv <- c(temp_mct_q[1], temp_mct_q[9], temp_mct_q[3], temp_mct_q[11], temp_mct_q[5], temp_mct_q[13], 
            temp_mct_q[7], temp_mct_q[15])
  B_stat <- c(mct_test_stat_sep_rr(x, h_mct)[[2]][i_r], "", mct_test_stat_sep_vv(x, h_mct)[[2]][i_r], "",
              mct_test_stat_sep_vn(x, h_mct)[[2]][i_r], "", mct_test_stat_sep_az(x, h_mct)[[2]][i_r], "")
  B_cv <- c(temp_mct_q[2], temp_mct_q[10], temp_mct_q[4], temp_mct_q[12], temp_mct_q[6], temp_mct_q[14], 
            temp_mct_q[8], temp_mct_q[16])
  mct_res <- data.frame(Contrast = c(paste("(", paste(h_i, collapse = ", "), ")", 
                                           collapse = "", sep = ""), rep("", 7)),
                        MCV = c("RR", "", "VV", "", "VN", "", "AZ", ""),
                        Method = rep(c("asym", "boot"), 4),
                        hC_est = c(c(h_i %*% c_est_rr), "", c(h_i %*% c_est_vv), "", c(h_i %*% c_est_vn), "", c(h_i %*% c_est_az), ""),
                        hC_lwr = c(c(c(h_i %*% c_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[1]),
                                  c(c(h_i %*% c_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[9]),
                                  c(c(h_i %*% c_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[3]),
                                  c(c(h_i %*% c_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[11]),
                                  c(c(h_i %*% c_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[5]),
                                  c(c(h_i %*% c_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[13]),
                                  c(c(h_i %*% c_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[7]),
                                  c(c(h_i %*% c_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[15])),
                        hC_upr = c(c(c(h_i %*% c_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[1]),
                                  c(c(h_i %*% c_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[9]),
                                  c(c(h_i %*% c_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[3]),
                                  c(c(h_i %*% c_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[11]), 
                                  c(c(h_i %*% c_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[5]), 
                                  c(c(h_i %*% c_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[13]),
                                  c(c(h_i %*% c_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[7]),
                                  c(c(h_i %*% c_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[15])),
                        hC_stat = C_stat,
                        hC_cv = C_cv,
                        hC_decision = ifelse(rep(as.numeric(C_stat[seq(1, length(C_stat) - 1, by = 2)]), each = 2) > C_cv, "H1", "H0"),
                        hB_est = c(c(h_i %*% b_est_rr), "", c(h_i %*% b_est_vv), "", c(h_i %*% b_est_vn), "", c(h_i %*% b_est_az), ""),
                        hB_lwr = c(c(c(h_i %*% b_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[2]),
                                  c(c(h_i %*% b_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[10]),
                                  c(c(h_i %*% b_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[4]),
                                  c(c(h_i %*% b_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[12]),
                                  c(c(h_i %*% b_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[6]),
                                  c(c(h_i %*% b_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[14]),
                                  c(c(h_i %*% b_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[8]),
                                  c(c(h_i %*% b_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[16])),
                        hB_upr = c(c(c(h_i %*% b_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[2]),
                                  c(c(h_i %*% b_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[10]),
                                  c(c(h_i %*% b_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[4]),
                                  c(c(h_i %*% b_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[12]),
                                  c(c(h_i %*% b_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[6]),
                                  c(c(h_i %*% b_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[14]),
                                  c(c(h_i %*% b_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[8]),
                                  c(c(h_i %*% b_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[16])),
                        hB_stat = B_stat,
                        hB_cv = B_cv,
                        hB_decision = ifelse(rep(as.numeric(B_stat[seq(1, length(B_stat) - 1, by = 2)]), each = 2) > B_cv, "H1", "H0"))
  
  if (r > 1) {
    for (i_r in 2:r) {
      h_i <- t(h_mct[i_r, ])
      C_stat <- c(mct_test_stat_sep_rr(x, h_mct)[[1]][i_r], "", mct_test_stat_sep_vv(x, h_mct)[[1]][i_r], "",
                  mct_test_stat_sep_vn(x, h_mct)[[1]][i_r], "", mct_test_stat_sep_az(x, h_mct)[[1]][i_r], "")
      C_cv <- c(temp_mct_q[1], temp_mct_q[9], temp_mct_q[3], temp_mct_q[11], temp_mct_q[5], temp_mct_q[13], 
                temp_mct_q[7], temp_mct_q[15])
      B_stat <- c(mct_test_stat_sep_rr(x, h_mct)[[2]][i_r], "", mct_test_stat_sep_vv(x, h_mct)[[2]][i_r], "",
                  mct_test_stat_sep_vn(x, h_mct)[[2]][i_r], "", mct_test_stat_sep_az(x, h_mct)[[2]][i_r], "")
      B_cv <- c(temp_mct_q[2], temp_mct_q[10], temp_mct_q[4], temp_mct_q[12], temp_mct_q[6], temp_mct_q[14], 
                temp_mct_q[8], temp_mct_q[16])
      mct_res <- rbind(mct_res, data.frame(Contrast = c(paste("(", paste(h_i, collapse = ", "), ")", 
                                                              collapse = "", sep = ""), rep("", 7)),
                                           MCV = c("RR", "", "VV", "", "VN", "", "AZ", ""),
                                           Method = rep(c("asym", "boot"), 4),
                                           hC_est = c(c(h_i %*% c_est_rr), "", c(h_i %*% c_est_vv), "", c(h_i %*% c_est_vn), "", c(h_i %*% c_est_az), ""),
                                           hC_lwr = c(c(c(h_i %*% c_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[1]),
                                                      c(c(h_i %*% c_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[9]),
                                                      c(c(h_i %*% c_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[3]),
                                                      c(c(h_i %*% c_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[11]),
                                                      c(c(h_i %*% c_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[5]),
                                                      c(c(h_i %*% c_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[13]),
                                                      c(c(h_i %*% c_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[7]),
                                                      c(c(h_i %*% c_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[15])),
                                           hC_upr = c(c(c(h_i %*% c_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[1]),
                                                      c(c(h_i %*% c_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[9]),
                                                      c(c(h_i %*% c_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[3]),
                                                      c(c(h_i %*% c_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[11]), 
                                                      c(c(h_i %*% c_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[5]), 
                                                      c(c(h_i %*% c_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[13]),
                                                      c(c(h_i %*% c_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[7]),
                                                      c(c(h_i %*% c_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[15])),
                                           hC_stat = C_stat,
                                           hC_cv = C_cv,
                                           hC_decision = ifelse(rep(as.numeric(C_stat[seq(1, length(C_stat) - 1, by = 2)]), each = 2) > C_cv, "H1", "H0"),
                                           hB_est = c(c(h_i %*% b_est_rr), "", c(h_i %*% b_est_vv), "", c(h_i %*% b_est_vn), "", c(h_i %*% b_est_az), ""),
                                           hB_lwr = c(c(c(h_i %*% b_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[2]),
                                                     c(c(h_i %*% b_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[10]),
                                                     c(c(h_i %*% b_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[4]),
                                                     c(c(h_i %*% b_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[12]),
                                                     c(c(h_i %*% b_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[6]),
                                                     c(c(h_i %*% b_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[14]),
                                                     c(c(h_i %*% b_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[8]),
                                                     c(c(h_i %*% b_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[16])),
                                           hB_upr = c(c(c(h_i %*% b_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[2]),
                                                     c(c(h_i %*% b_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[10]),
                                                     c(c(h_i %*% b_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[4]),
                                                     c(c(h_i %*% b_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[12]),
                                                     c(c(h_i %*% b_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[6]),
                                                     c(c(h_i %*% b_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[14]),
                                                     c(c(h_i %*% b_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[8]),
                                                     c(c(h_i %*% b_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[16])),
                                           hB_stat = B_stat,
                                           hB_cv = B_cv,
                                           hB_decision = ifelse(rep(as.numeric(B_stat[seq(1, length(B_stat) - 1, by = 2)]), each = 2) > B_cv, "H1", "H0")))
      
      # asymptotic
      # c(c(h_i %*% c_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[1])
      # c(c(h_i %*% c_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[1])
      # c(c(h_i %*% b_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[2])
      # c(c(h_i %*% b_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[2])
      # 
      # c(c(h_i %*% c_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[3])
      # c(c(h_i %*% c_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[3])
      # c(c(h_i %*% b_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[4])
      # c(c(h_i %*% b_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[4])
      # 
      # c(c(h_i %*% c_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[5])
      # c(c(h_i %*% c_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[5])
      # c(c(h_i %*% b_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[6])
      # c(c(h_i %*% b_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[6])
      # 
      # c(c(h_i %*% c_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[7])
      # c(c(h_i %*% c_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[7])
      # c(c(h_i %*% b_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[8])
      # c(c(h_i %*% b_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[8])
      # 
      # bootstrap
      # c(c(h_i %*% c_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[9])
      # c(c(h_i %*% c_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_c_est %*% t(h_i)) * temp_mct_q[9])
      # c(c(h_i %*% b_est_rr) - nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[10])
      # c(c(h_i %*% b_est_rr) + nn_m05 * sqrt(h_i %*% sigma_rr_b_est %*% t(h_i)) * temp_mct_q[10])
      # 
      # c(c(h_i %*% c_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[11])
      # c(c(h_i %*% c_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_c_est %*% t(h_i)) * temp_mct_q[11])
      # c(c(h_i %*% b_est_vv) - nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[12])
      # c(c(h_i %*% b_est_vv) + nn_m05 * sqrt(h_i %*% sigma_vv_b_est %*% t(h_i)) * temp_mct_q[12])
      # 
      # c(c(h_i %*% c_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[13])
      # c(c(h_i %*% c_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_c_est %*% t(h_i)) * temp_mct_q[13])
      # c(c(h_i %*% b_est_vn) - nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[14])
      # c(c(h_i %*% b_est_vn) + nn_m05 * sqrt(h_i %*% sigma_vn_b_est %*% t(h_i)) * temp_mct_q[14])
      # 
      # c(c(h_i %*% c_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[15])
      # c(c(h_i %*% c_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_c_est %*% t(h_i)) * temp_mct_q[15])
      # c(c(h_i %*% b_est_az) - nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[16])
      # c(c(h_i %*% b_est_az) + nn_m05 * sqrt(h_i %*% sigma_az_b_est %*% t(h_i)) * temp_mct_q[16])
    }
  }
  
  if (dd > 1) {
    res <- list(overall_res = list(test_stat = overall_test_stat,
                                   p_values = overall_p_val),
                mct_res = mct_res,
                h_mct = h_mct,
                h_wald = h_wald,
                alpha = alpha)
  } else if (dd == 1) {
    overall_test_stat <- overall_test_stat[, 1:2]
    colnames(overall_test_stat) <- c("C", "B")
    overall_p_val <- overall_p_val[, 1:2]
    colnames(overall_p_val) <- c("C", "B")
    temp <- c()
    vec <- which(mct_res$MCV == "RR")
    for (i_t in seq_len(length(vec))) {
      temp <- c(temp, vec[i_t], vec[i_t] + 1)
    }
    mct_res <- mct_res[temp, ]
    mct_res$MCV <- NULL
    rownames(mct_res) <- seq_len(2 * r)
    res <- list(overall_res = list(test_stat = overall_test_stat,
                                   p_values = overall_p_val),
                mct_res = mct_res,
                h_mct = h_mct,
                h_wald = h_wald,
                alpha = alpha)
  }
  class(res) <- "gfdmcv"
  return(res)
}
