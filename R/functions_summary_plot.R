#' Print "gfdmcv" object
#' 
#' Prints the summary of the inference methods for CV and MCVs.
#' 
#' @param object an "gfdmcv" object.
#' @param ... integer indicating the number of decimal places to be used to present the numerical results, 
#' It can be named \code{digits} as in the \code{round()} function (see examples).
#' 
#' @details The function prints out the information about the significance level, constrast matrices, 
#' test statistics, \eqn{p}-values, critical values, simultaneous confidence intervals for contrasts
#' performed by the \code{GFDmcv()} function.
#' 
#' @return No return value, called for side effects.
#' 
#' @examples
#' # Some of the examples may run some time. 
#' # For more examples, see the documentation of the GFDmcv() function.
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
#' }
#' 
#' @export

summary.gfdmcv <- function(object, ...) {
  if (ncol(object$mct_res) == 15) {
    cat("#--- Inference for multivariate coefficients of variation and their reciprocals ---#", "\n", "\n")
    cat("- Significance level:", object$alpha)
    cat("\n", "\n")
    cat("#--- Constrasts for Wald-type tests -----------------------------------------------#", "\n")
    print(object$h_wald)
    cat("\n")
    cat("#--- Constrasts for multiple contrasts tests --------------------------------------#", "\n")
    print(object$h_mct)
    cat("\n")
    cat("#--- Overall test statistics ------------------------------------------------------#", "\n")
    print(round(object$overall_res$test_stat, ...))
    cat("\n")
    cat("#--- Overall p-values -------------------------------------------------------------#", "\n")
    print(round(object$overall_res$p_values, ...))
    cat("\n")
    cat("#--- Multiple contrast testing results --------------------------------------------#", "\n")
    temp_mct_res <- object$mct_res
    for (i in seq_len(nrow(object$mct_res))) {
      for (j in c(4:8, 10:14)) {
        if (object$mct_res[i, j] != "") {
          # temp_mct_res[i, j] <- round(as.numeric(stringr::str_extract_all(object$mct_res[i, j], "-*\\d+\\.*\\d*")[[1]]), ...)
          temp_mct_res[i, j] <- round(as.numeric(object$mct_res[i, j]), ...)
        }
      }
    }
    print(temp_mct_res)
    cat("#----------------------------------------------------------------------------------#", "\n")
  } else if (ncol(object$mct_res) == 14) {
    cat("#--- Inference for coefficient of variation and its reciprocal --------------------#", "\n", "\n")
    cat("- Significance level:", object$alpha)
    cat("\n", "\n")
    cat("#--- Constrasts for Wald-type tests -----------------------------------------------#", "\n")
    print(object$h_wald)
    cat("\n")
    cat("#--- Constrasts for multiple contrasts tests --------------------------------------#", "\n")
    print(object$h_mct)
    cat("\n")
    cat("#--- Overall test statistics ------------------------------------------------------#", "\n")
    print(round(object$overall_res$test_stat, ...))
    cat("\n")
    cat("#--- Overall p-values -------------------------------------------------------------#", "\n")
    print(round(object$overall_res$p_values, ...))
    cat("\n")
    cat("#--- Multiple contrast testing results --------------------------------------------#", "\n")
    temp_mct_res <- object$mct_res
    for (i in seq_len(nrow(object$mct_res))) {
      for (j in c(3:7, 9:13)) {
        if (object$mct_res[i, j] != "") {
          # temp_mct_res[i, j] <- round(as.numeric(stringr::str_extract_all(object$mct_res[i, j], "-*\\d+\\.*\\d*")[[1]]), ...)
          temp_mct_res[i, j] <- round(as.numeric(object$mct_res[i, j]), ...)
        }
      }
    }
    print(temp_mct_res)
    cat("#----------------------------------------------------------------------------------#", "\n")
  }
}

#' Plot simultaneous confidence intervals for contrasts
#' 
#' Simultaneous confidence intervals for contrasts for CV and MCVs and their reciprocals are plotted. 
#' 
#' @param x an "gfdmcv" object.
#' @param ... additional arguments not used.
#' 
#' @return No return value, called for side effects.
#' 
#' @examples
#' # Some of the examples may run some time. 
#' # For more examples, see the documentation of the GFDmcv() function.
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
#' oldpar <- par(mar = c(4, 5, 2, 0.3))
#' plot(res)
#' par(oldpar)
#' }
#' 
#' @importFrom graphics axis lines par segments
#' 
#' @export

plot.gfdmcv <- function(x, ...) {
  contr <- x$mct_res$Contrast[x$mct_res$Contrast != ""]
  contr_2 <- rep(contr, each = 2)
  r <- length(contr)
  temp_r <- 2 * r + 1
  
  if (ncol(x$mct_res) == 15) {
    mcvs <- c("RR", "VV", "VN", "AZ")
    x$mct_res$hC_est_2 <- rep(x$mct_res$hC_est[x$mct_res$hC_est != ""], each = 2)
    x$mct_res$hB_est_2 <- rep(x$mct_res$hB_est[x$mct_res$hB_est != ""], each = 2)
    x$mct_res$mcv_2 <- rep(c("RR", "VV", "VN", "AZ"), each = 2)
    for (i_r in seq_len(4)) {
      temp <- x$mct_res[x$mct_res$mcv_2 == mcvs[i_r], ]
      # C
      m1 <- min(temp$hC_lwr)
      m2 <- max(temp$hC_upr)
      x_lim <- c(m1, m2)
      if (m1 > 0) {
        x_lim[1] <- 0
      } else if (m2 < 0) {
        x_lim[2] <- 0
      }
      y_lim <- c(0.8, 2 * r + 0.2)
      plot(x = x$mct_res$hC_lwr, y = x$mct_res$hC_upr,
           xlim = x_lim, ylim = y_lim,
           main = paste(paste(mcvs[i_r], ",", sep = ""), 
                        paste(100 * (1 - x$alpha), "%", sep = ""), 
                        "simultaneous confidence intervals for hC"),
           xlab = "Contrast hC", ylab = "", type = "n", yaxt = "n")
      for (i_y in seq_len(2 * r)) {
        if ((2 * r - i_y + 1) %% 2 == 1) {
          axis(2, at = temp_r - i_y, las = 2,
               labels = paste(contr_2[i_y], "\n", "boot", sep = ""))
        } else {
          axis(2, at = temp_r - i_y, las = 2,
               labels = paste(contr_2[i_y], "\n", "asym", sep = ""))
        }
      }
      lines(x = rep(0, 100051), y = -50:100000, lwd = 2, lty = 2)
      for (i_p in 1:(2 * r)) {
        lines(x = seq(x_lim[1] - 100000, x_lim[2] + 100000, length.out = 1000), 
              y = rep(i_p, 1000), lwd = 0.5, lty = 3)
      }
      for (i_c in 1:(2 * r)) {
        segments(temp$hC_lwr[i_c], temp_r - i_c, temp$hC_upr[i_c], temp_r - i_c)
        segments(temp$hC_lwr[i_c], temp_r - i_c - 0.1, temp$hC_lwr[i_c], temp_r - i_c + 0.1)
        segments(temp$hC_upr[i_c], temp_r - i_c - 0.1, temp$hC_upr[i_c], temp_r - i_c + 0.1)
        segments(as.numeric(temp$hC_est_2[i_c]), temp_r - i_c - 0.1, as.numeric(temp$hC_est_2[i_c]), temp_r - i_c + 0.1)
      }
      # B
      m1 <- min(temp$hB_lwr)
      m2 <- max(temp$hB_upr)
      x_lim <- c(m1, m2)
      if (m1 > 0) {
        x_lim[1] <- 0
      } else if (m2 < 0) {
        x_lim[2] <- 0
      }
      y_lim <- c(0.8, 2 * r + 0.2)
      plot(x = x$mct_res$hB_lwr, y = x$mct_res$hB_upr,
           xlim = x_lim, ylim = y_lim,
           main = paste(paste(mcvs[i_r], ",", sep = ""), 
                        paste(100 * (1 - x$alpha), "%", sep = ""), 
                        "simultaneous confidence intervals for hB"),
           xlab = "Contrast hB", ylab = "", type = "n", yaxt = "n")
      for (i_y in seq_len(2 * r)) {
        if ((temp_r - i_y) %% 2 == 1) {
          axis(2, at = temp_r - i_y, las = 2,
               labels = paste(contr_2[i_y], "\n", "boot", sep = ""))
        } else {
          axis(2, at = temp_r - i_y, las = 2,
               labels = paste(contr_2[i_y], "\n", "asym", sep = ""))
        }
      }
      lines(x = rep(0, 100051), y = -50:100000, lwd = 2, lty = 2)
      for (i_p in 1:(2 * r)) {
        lines(x = seq(x_lim[1] - 100000, x_lim[2] + 100000, length.out = 1000), 
              y = rep(i_p, 1000), lwd = 0.5, lty = 3)
      }
      for (i_c in 1:(2 * r)) {
        segments(temp$hB_lwr[i_c], temp_r - i_c, temp$hB_upr[i_c], temp_r - i_c)
        segments(temp$hB_lwr[i_c], temp_r - i_c - 0.1, temp$hB_lwr[i_c], temp_r - i_c + 0.1)
        segments(temp$hB_upr[i_c], temp_r - i_c - 0.1, temp$hB_upr[i_c], temp_r - i_c + 0.1)
        segments(as.numeric(temp$hB_est_2[i_c]), temp_r - i_c - 0.1, as.numeric(temp$hB_est_2[i_c]), temp_r - i_c + 0.1)
      }
    }
  } else if (ncol(x$mct_res) == 14) {
    x$mct_res$hC_est_2 <- rep(x$mct_res$hC_est[x$mct_res$hC_est != ""], each = 2)
    x$mct_res$hB_est_2 <- rep(x$mct_res$hB_est[x$mct_res$hB_est != ""], each = 2)
    temp <- x$mct_res
    # C
    m1 <- min(temp$hC_lwr)
    m2 <- max(temp$hC_upr)
    x_lim <- c(m1, m2)
    if (m1 > 0) {
      x_lim[1] <- 0
    } else if (m2 < 0) {
      x_lim[2] <- 0
    }
    y_lim <- c(0.8, 2 * r + 0.2)
    plot(x = x$mct_res$hC_lwr, y = x$mct_res$hC_upr,
         xlim = x_lim, ylim = y_lim,
         main = paste("CV,", 
                      paste(100 * (1 - x$alpha), "%", sep = ""), 
                      "simultaneous confidence intervals for hC"),
         xlab = "Contrast hC", ylab = "", type = "n", yaxt = "n")
    for (i_y in seq_len(2 * r)) {
      if ((2 * r - i_y + 1) %% 2 == 1) {
        axis(2, at = temp_r - i_y, las = 2,
             labels = paste(contr_2[i_y], "\n", "boot", sep = ""))
      } else {
        axis(2, at = temp_r - i_y, las = 2,
             labels = paste(contr_2[i_y], "\n", "asym", sep = ""))
      }
      
    }
    lines(x = rep(0, 100051), y = -50:100000, lwd = 2, lty = 2)
    for (i_p in 1:(2 * r)) {
      lines(x = seq(x_lim[1] - 100000, x_lim[2] + 100000, length.out = 1000), 
            y = rep(i_p, 1000), lwd = 0.5, lty = 3)
    }
    for (i_c in 1:(2 * r)) {
      segments(temp$hC_lwr[i_c], temp_r - i_c, temp$hC_upr[i_c], temp_r - i_c)
      segments(temp$hC_lwr[i_c], temp_r - i_c - 0.1, temp$hC_lwr[i_c], temp_r - i_c + 0.1)
      segments(temp$hC_upr[i_c], temp_r - i_c - 0.1, temp$hC_upr[i_c], temp_r - i_c + 0.1)
      segments(as.numeric(temp$hC_est_2[i_c]), temp_r - i_c - 0.1, as.numeric(temp$hC_est_2[i_c]), temp_r - i_c + 0.1)
    }
    # B
    m1 <- min(temp$hB_lwr)
    m2 <- max(temp$hB_upr)
    x_lim <- c(m1, m2)
    if (m1 > 0) {
      x_lim[1] <- 0
    } else if (m2 < 0) {
      x_lim[2] <- 0
    }
    y_lim <- c(0.8, 2 * r + 0.2)
    plot(x = x$mct_res$hB_lwr, y = x$mct_res$hB_upr,
         xlim = x_lim, ylim = y_lim,
         main = paste("CV,", 
                      paste(100 * (1 - x$alpha), "%", sep = ""), 
                      "simultaneous confidence intervals for hB"),
         xlab = "Contrast hB", ylab = "", type = "n", yaxt = "n")
    for (i_y in seq_len(2 * r)) {
      if ((2 * r - i_y + 1) %% 2 == 1) {
        axis(2, at = temp_r - i_y, las = 2,
             labels = paste(contr_2[i_y], "\n", "boot", sep = ""))
      } else {
        axis(2, at = temp_r - i_y, las = 2,
             labels = paste(contr_2[i_y], "\n", "asym", sep = ""))
      }
      
    }
    lines(x = rep(0, 100051), y = -50:100000, lwd = 2, lty = 2)
    for (i_p in 1:(2 * r)) {
      lines(x = seq(x_lim[1] - 100000, x_lim[2] + 100000, length.out = 1000), 
            y = rep(i_p, 1000), lwd = 0.5, lty = 3)
    }
    for (i_c in 1:(2 * r)) {
      segments(temp$hB_lwr[i_c], temp_r - i_c, temp$hB_upr[i_c], temp_r - i_c)
      segments(temp$hB_lwr[i_c], temp_r - i_c - 0.1, temp$hB_lwr[i_c], temp_r - i_c + 0.1)
      segments(temp$hB_upr[i_c], temp_r - i_c - 0.1, temp$hB_upr[i_c], temp_r - i_c + 0.1)
      segments(as.numeric(temp$hB_est_2[i_c]), temp_r - i_c - 0.1, as.numeric(temp$hB_est_2[i_c]), temp_r - i_c + 0.1)
    }
  }
}
