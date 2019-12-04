// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// knn_mi
double knn_mi(arma::mat datax, arma::mat datay, int k);
RcppExport SEXP _fastmit_knn_mi(SEXP dataxSEXP, SEXP dataySEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type datax(dataxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type datay(dataySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(knn_mi(datax, datay, k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastmit_knn_mi", (DL_FUNC) &_fastmit_knn_mi, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastmit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}