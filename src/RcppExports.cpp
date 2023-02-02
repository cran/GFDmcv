// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// D_tilde_rcpp
NumericMatrix D_tilde_rcpp(NumericVector x);
RcppExport SEXP _GFDmcv_D_tilde_rcpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(D_tilde_rcpp(x));
    return rcpp_result_gen;
END_RCPP
}
// Psi_est_rcpp
List Psi_est_rcpp(NumericMatrix x);
RcppExport SEXP _GFDmcv_Psi_est_rcpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Psi_est_rcpp(x));
    return rcpp_result_gen;
END_RCPP
}
// sigma_est_w_rccp
NumericMatrix sigma_est_w_rccp(NumericMatrix x, NumericVector w, NumericVector mu_w);
RcppExport SEXP _GFDmcv_sigma_est_w_rccp(SEXP xSEXP, SEXP wSEXP, SEXP mu_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu_w(mu_wSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_est_w_rccp(x, w, mu_w));
    return rcpp_result_gen;
END_RCPP
}
// sigma_est_ws_rccp
NumericMatrix sigma_est_ws_rccp(NumericMatrix x, double w, NumericVector mu_w);
RcppExport SEXP _GFDmcv_sigma_est_ws_rccp(SEXP xSEXP, SEXP wSEXP, SEXP mu_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu_w(mu_wSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_est_ws_rccp(x, w, mu_w));
    return rcpp_result_gen;
END_RCPP
}
// Psi_est_w_rcpp
List Psi_est_w_rcpp(NumericMatrix x, NumericVector w);
RcppExport SEXP _GFDmcv_Psi_est_w_rcpp(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(Psi_est_w_rcpp(x, w));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GFDmcv_D_tilde_rcpp", (DL_FUNC) &_GFDmcv_D_tilde_rcpp, 1},
    {"_GFDmcv_Psi_est_rcpp", (DL_FUNC) &_GFDmcv_Psi_est_rcpp, 1},
    {"_GFDmcv_sigma_est_w_rccp", (DL_FUNC) &_GFDmcv_sigma_est_w_rccp, 3},
    {"_GFDmcv_sigma_est_ws_rccp", (DL_FUNC) &_GFDmcv_sigma_est_ws_rccp, 3},
    {"_GFDmcv_Psi_est_w_rcpp", (DL_FUNC) &_GFDmcv_Psi_est_w_rcpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_GFDmcv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}