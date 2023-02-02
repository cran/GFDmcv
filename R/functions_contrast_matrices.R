#' Contrasts' matrices
#' 
#' The output are different contrast matrices, namely the centering matrix as well as the matrices of Dunnett's and Tukey's contrasts for given number of groups.
#' 
#' @param k an integer denoting a number of groups
#' @param type an character denoting type of contrasts. The possible values are \code{"center"} (default), 
#' \code{"Dunnett"}, \code{"Tukey"}.
#' 
#' @details The centering matrix is \eqn{\mathbf{P}_k = \mathbf{I}_k  - \mathbf{J}_k/k}, where 
#' \eqn{\mathbf{I}_k} is the unity matrix and \eqn{\mathbf{J}_k=\mathbf{1}\mathbf{1}^{\top} \in R^{k\times k}} consisting of \eqn{1}'s only.
#' The matrix of Dunnett's contrasts: \deqn{\left(\begin{array}{rrrrrr} 
#' -1 & 1 & 0 & \ldots & \ldots  & 0\\			
#' -1 & 0 & 1 & 0 & \ldots  & 0 \\
#' \vdots & \vdots & \vdots & \vdots & \vdots & \vdots  \\
#' -1 & 0 & \ldots & \ldots & 0  &  1\\
#' \end{array}\right)\in R^{(k-1)\times k}.} 
#' The matrix of Tukey's contrasts: \deqn{\left(\begin{array}{rrrrrrr}
#' -1 & 1 & 0 & \ldots & \ldots & 0 & 0\\		
#' -1 & 0 & 1 & 0 & \ldots & \ldots & 0 \\
#' \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
#' -1 & 0 & 0 & 0 & \ldots  & \ldots & 1 \\
#' 0 & -1 & 1 & 0 &\ldots & \ldots & 0 \\
#' \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
#' 0 & \ldots & \ldots & \ldots & 0 & -1 &  1\\
#' \end{array}\right)\in R^{{k \choose 2}\times k}.}
#'
#' @return The matrix of contrasts.
#' 
#' @examples
#' contr_mat(4, type = "center")
#' contr_mat(4, type = "Dunnett")
#' contr_mat(4, type = "Tukey")
#'
#' @references Ditzhaus M., Smaga L. (2022) Permutation test for the multivariate coefficient of variation in factorial designs. 
#' Journal of Multivariate Analysis 187, 104848.
#' 
#' Ditzhaus M., Smaga L. (2023) Inference for all variants of the multivariate coefficient of variation in factorial designs. Preprint https://arxiv.org/abs/2301.12009.
#' 
#' Dunnett C. (1955) A multiple comparison procedure for comparing several treatments with a control. Journal of the American Statistical Association 50, 1096-1121.
#' 
#' Tukey J.W. (1953) The problem of multiple comparisons. Princeton University.
#' 
#' @export

contr_mat <- function(k, type = c("center", "Dunnett", "Tukey")) {
  type <- match.arg(type)
  h_Dunnett <- function(k) {
    comp <- k - 1
    temp <- cbind(rep(-1, comp), diag(rep(1, comp)))
    return(temp)
  }
  if (type == "center") {
    return(diag(k) - (1 / k) * matrix(1, k, k))
  } else if (type == "Dunnett") {
    return(h_Dunnett(k))
  } else if (type == "Tukey") {
    temp <- h_Dunnett(k)
    if (k > 2) {
      comp <- choose(k, 2)
      j <- 1
      for (i in (k - 1):2) {
        temp <- rbind(temp, 
                      cbind(matrix(0, nrow = k - 1 - j, ncol = j), h_Dunnett(i)))
        j <- j + 1
      }
    }
    return(temp)
  }
}
