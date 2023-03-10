% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_testing_all_mcv.R
\name{GFDmcv}
\alias{GFDmcv}
\title{Inference for four multivariate coefficients of variation}
\usage{
GFDmcv(
  x,
  h_mct,
  h_wald,
  alpha = 0.05,
  n_perm = 1000,
  n_boot = 1000,
  parallel = FALSE,
  n_cores = NULL
)
}
\arguments{
\item{x}{a list of length \eqn{k} with elements being \eqn{n_i\times d} matrices of data, \eqn{i=1,\dots,k}.}

\item{h_mct}{a \eqn{r\times k} contrast matrix \eqn{\mathbf{H}} of full row rank for multiple contrast tests. 
Remember to specify it correctly taking into account the order of elements of the list \code{x}.}

\item{h_wald}{a \eqn{q\times k} contrast matrix \eqn{\mathbf{H}} of full row rank for the Wald-type tests. 
Remember to specify it correctly taking into account the order of elements of the list \code{x}.}

\item{alpha}{a significance level (then \code{1-alpha} is the confidence level).}

\item{n_perm}{a number of permutation replicates.}

\item{n_boot}{a number of bootstrap replicates.}

\item{parallel}{a logical indicating whether to use parallelization.}

\item{n_cores}{if \code{parallel = TRUE}, a number of processes used in parallel computation. 
Its default value means that it will be equal to a number of cores of a computer used.}
}
\value{
A list of class \code{gfdmcv} containing the following components:
\item{overall_res}{a list of two elements representing the results for testing
the global null hypothesis. The first one is a matrix \code{test_stat} of the
test statistics, while the second is a matrix \code{p_values} of the respective \eqn{p}-values.}
\item{mct_res}{all results of MCT tests for particular hypothesis in \code{h_mct}, i.e., 
the estimators and simultaneous confidence intervals for \eqn{\mathbf{h_{\ell}^{\top}\pmb{C}}} and 
for \eqn{\mathbf{h_{\ell}^{\top}\pmb{B}}}, the test statistics and critical values as well as the decisions.}
\item{h_mct}{an argument \code{h_mct}.}
\item{h_wald}{an argument \code{h_wald}.}
\item{alpha}{an argument \code{alpha}.}
}
\description{
The function \code{GFDmcv()} calculates the Wald-type statistic for global null hypotheses 
and max-type statistics for multiple local null hypotheses, both in terms of the four variants 
of the multivariate coefficient of variation. Respective \eqn{p}-values
are obtained by a \eqn{\chi^2}-approximation, a pooled bootstrap strategy and a pooled permutation approach (only for the
Wald-type statistic), respectively.
}
\details{
The function \code{GFDmcv()} calculates the Wald-type statistic for global null hypotheses of the form
\deqn{ \mathcal H_0: \mathbf{H} (C_1,\ldots,C_k)^\top = \mathbf{0}\,\,\text{and}\,\,\mathcal H_0: \mathbf{H} (B_1,\ldots,B_k)^\top = \mathbf{0},}
where \eqn{\mathbf{H}} is a contrast matrix reflecting the research question of interest and
\eqn{C_i} (\eqn{B_i}) are the subgroup-specific MCVs (and their reciprocal) by Reyment (1960, RR), Van Valen (1974, VV), 
Voinov and Nikulin (1996, VN) or Albert and Zhang (2010, AZ), respectively. 
We refer to the function \code{e_mcv()} for the detailed definitions of the different
variants. The \eqn{p}-value of the Wald-type statistic relies on a \eqn{\chi^2}-approximation,
a (pooled) bootstrap or permutation approach.

Furthermore, the function \code{GFDmcv()} calculates a max-type test statistic for the multiple
comparison of \eqn{q} local null hypotheses:
 \deqn{ \mathcal H_{0,\ell}: \mathbf{h_\ell}^\top \mathbf{C} = \mathbf{0}\,\,
 \text{or}\,\,\mathcal H_{0,\ell}: \mathbf{h_\ell}^\top \mathbf{B} = \mathbf{0}, \,\,\ell=1,\ldots,q,}
where \eqn{\mathbf{C}=(C_1,\ldots,C_k)^\top} and  \eqn{\mathbf{B}=(B_1,\ldots,B_k)^\top}. The \eqn{p}-values are determined by a Gaussian approximation and a bootstrap approach, respectively.
In addition to the local test decisions, multiple adjusted confidence intervals for the contrasts \eqn{\mathbf{h_{\ell}^{\top}\pmb{C}}} and
\eqn{\mathbf{h_{\ell}^{\top}\pmb{B}}}, respectively, are calculated.

Please have a look on the plot and summary functions designed for this package. They can be used
to simplify the output of \code{GFDmcv()}.
}
\examples{
# Some of the examples may run some time.

# one-way analysis for MCV and CV
# d > 1 (MCV)
data_set <- lapply(list(iris[iris$Species == "setosa", 1:3],
                        iris[iris$Species == "versicolor", 1:3],
                        iris[iris$Species == "virginica", 1:3]),
                   as.matrix)
# estimators and confidence intervals of MCVs and their reciprocals
lapply(data_set, e_mcv)
# contrast matrices
k <- length(data_set)
# Tukey's contrast matrix
h_mct <- contr_mat(k, type = "Tukey")
# centering matrix P_k
h_wald <- contr_mat(k, type = "center")
\donttest{
# testing without parallel computing
res <- GFDmcv(data_set, h_mct, h_wald)
summary(res, digits = 3)
oldpar <- par(mar = c(4, 5, 2, 0.3))
plot(res)
par(oldpar)}
\donttest{
# testing with parallel computing
library(doParallel)
res <- GFDmcv(data_set, h_mct, h_wald, parallel = TRUE, n_cores = 2)
summary(res, digits = 3)
oldpar <- par(mar = c(4, 5, 2, 0.3))
plot(res)
par(oldpar)}
\dontshow{
data_set_2 <- lapply(data_set, function(x) x[1:5, ])
res <- GFDmcv(data_set, h_mct, h_wald, n_perm = 2, n_boot = 2)
summary(res, digits = 3)
oldpar <- par(mar = c(4, 5, 2, 0.3))
plot(res)
par(oldpar)}
# two-way analysis for CV (based on the example in Ditzhaus and Smaga, 2022)
library(HSAUR)
data_set <- lapply(list(BtheB$bdi.pre[BtheB$drug == "No" & BtheB$length == "<6m"],
                        BtheB$bdi.pre[BtheB$drug == "No" & BtheB$length == ">6m"],
                        BtheB$bdi.pre[BtheB$drug == "Yes" & BtheB$length == "<6m"],
                        BtheB$bdi.pre[BtheB$drug == "Yes" & BtheB$length == ">6m"]), 
                   as.matrix)
# estimators and confidence intervals of CV and its reciprocal
lapply(data_set, e_mcv)
\donttest{
# interaction
h_mct <- contr_mat(4, type = "Tukey")
h_wald <- kronecker(contr_mat(2, type = "center"), 
                    contr_mat(2, type = "center"))
res <- GFDmcv(data_set, h_mct, h_wald)
summary(res, digits = 3)
oldpar <- par(mar = c(4, 6, 2, 0.1))
plot(res)
par(oldpar)}
\donttest{
# main effect drug
h_mct <- matrix(c(1, 1, -1, -1), nrow = 1)
h_wald <- kronecker(contr_mat(2, type = "center"), 0.5 * matrix(1, 1, 2))
res <- GFDmcv(data_set, h_mct, h_wald)
summary(res, digits = 3)
oldpar <- par(mar = c(4, 6, 2, 0.1))
plot(res)
par(oldpar)}
\donttest{
# main effect length
h_mct <- matrix(c(1, -1, 1, -1), nrow = 1)
h_wald <- kronecker(0.5 * matrix(1, 1, 2), contr_mat(2, type = "center"))
res <- GFDmcv(data_set, h_mct, h_wald)
summary(res, digits = 3)
oldpar <- par(mar = c(4, 6, 2, 0.1))
plot(res)
par(oldpar)}

}
\references{
Albert A., Zhang L. (2010) A novel definition of the multivariate coefficient of variation. Biometrical Journal 52:667-675.

Ditzhaus M., Smaga L. (2022) Permutation test for the multivariate coefficient of variation in factorial designs. 
Journal of Multivariate Analysis 187, 104848.

Ditzhaus M., Smaga L. (2023) Inference for all variants of the multivariate coefficient of variation in factorial designs. Preprint https://arxiv.org/abs/2301.12009.

Reyment R.A. (1960) Studies on Nigerian Upper Cretaceous and Lower Tertiary Ostracoda: part 1. 
Senonian and Maastrichtian Ostracoda, Stockholm Contributions in Geology, vol 7.

Van Valen L. (1974) Multivariate structural statistics in natural history. Journal of Theoretical Biology 45:235-247.

Voinov V., Nikulin M. (1996) Unbiased Estimators and Their Applications, Vol. 2, Multivariate Case. Kluwer, Dordrecht.
}
