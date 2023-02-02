#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix D_tilde_rcpp(NumericVector x) {
  double d = x.length();
  NumericMatrix mat_temp(pow(d, 2), d);
  for (int r = 0; r < d; r++) {
    for (int a = 0; a < d; a++) {
      for (int s = 0; s < d; s++) {
        if ((s == a) & (a != r)) {
          mat_temp(r + a * d, s) = -x(r);
        }
        if ((s == a) & (a == r)) {
          mat_temp(r + a * d, s) = -2 * x(s);
        }
        if ((r == s) & (s != a)) {
          mat_temp(r + a * d, s) = -x(a);
        }
      }
    }
  }
  return mat_temp;
}

// [[Rcpp::export]]
List Psi_est_rcpp(NumericMatrix x) {
  double d = x.ncol();
  NumericMatrix Psi_3_est(pow(d, 2), d);
  NumericMatrix Psi_4_est(pow(d, 2), pow(d, 2));
  for (int a = 0; a < d; a++) {
    for (int r = 0; r < d; r++) {
      for (int s = 0; s < d; s++) {
        Psi_3_est(a * d + r, s) = mean(x(_, a) * x(_, r) * x(_, s)) - 
          mean(x(_, a) * x(_, r)) * mean(x(_, s));
        for (int b = 0; b < d; b++) {
          Psi_4_est(a * d + r, b * d + s) = mean(x(_, a) * x(_, r) * x(_, b) * x(_, s)) - 
            mean(x(_, a) * x(_, r)) * mean(x(_, b) * x(_, s));
        }
      }
    }
  }
  return List::create(Psi_3_est, Psi_4_est);
}

// [[Rcpp::export]]
NumericMatrix sigma_est_w_rccp(NumericMatrix x, NumericVector w, NumericVector mu_w) {
  double n_i = x.nrow();
  double d = x.ncol();
  NumericMatrix temp(d, d);
  for (int r = 0; r < d; r++) {
    for (int s = 0; s < d; s++) {
      for (int j = 0; j < n_i; j++) {
        temp(r, s) = temp(r, s) + w(j) * x(j, r) * x(j, s) - mu_w(r) * mu_w(s);
      }
    }
  }
  temp = temp / (n_i - 1);
  return temp;
}

// [[Rcpp::export]]
NumericMatrix sigma_est_ws_rccp(NumericMatrix x, double w, NumericVector mu_w) {
  double n_i = x.nrow();
  double d = x.ncol();
  NumericMatrix temp(d, d);
  for (int r = 0; r < d; r++) {
    for (int s = 0; s < d; s++) {
      for (int j = 0; j < n_i; j++) {
        temp(r, s) = temp(r, s) + w * x(j, r) * x(j, s) - mu_w(r) * mu_w(s);
      }
    }
  }
  temp = temp / (n_i - 1);
  return temp;
}

// [[Rcpp::export]]
List Psi_est_w_rcpp(NumericMatrix x, NumericVector w) {
  double d = x.ncol();
  NumericMatrix Psi_3_est(pow(d, 2), d);
  NumericMatrix Psi_4_est(pow(d, 2), pow(d, 2));
  for (int a = 0; a < d; a++) {
    for (int r = 0; r < d; r++) {
      for (int s = 0; s < d; s++) {
        Psi_3_est(a * d + r, s) = mean(w * x(_, a) * x(_, r) * x(_, s)) - 
          mean(w * x(_, a) * x(_, r)) * mean(w * x(_, s));
        for (int b = 0; b < d; b++) {
          Psi_4_est(a * d + r, b * d + s) = mean(w * x(_, a) * x(_, r) * x(_, b) * x(_, s)) - 
            mean(w * x(_, a) * x(_, r)) * mean(w * x(_, b) * x(_, s));
        }
      }
    }
  }
  return List::create(Psi_3_est, Psi_4_est);
}
