// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mfwht
SEXP mfwht(arma::mat data, int nrow, int ncol);
RcppExport SEXP _eSCAN_mfwht(SEXP dataSEXP, SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(mfwht(data, nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}
// big_mfwht
SEXP big_mfwht(arma::mat data);
RcppExport SEXP _eSCAN_big_mfwht(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(big_mfwht(data));
    return rcpp_result_gen;
END_RCPP
}
// compPval
SEXP compPval(double Q, arma::mat Covw_sub);
RcppExport SEXP _eSCAN_compPval(SEXP QSEXP, SEXP Covw_subSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Covw_sub(Covw_subSEXP);
    rcpp_result_gen = Rcpp::wrap(compPval(Q, Covw_sub));
    return rcpp_result_gen;
END_RCPP
}
// compCov
SEXP compCov(arma::mat G, arma::mat X, arma::vec working, double sigma, int fam);
RcppExport SEXP _eSCAN_compCov(SEXP GSEXP, SEXP XSEXP, SEXP workingSEXP, SEXP sigmaSEXP, SEXP famSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type working(workingSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type fam(famSEXP);
    rcpp_result_gen = Rcpp::wrap(compCov(G, X, working, sigma, fam));
    return rcpp_result_gen;
END_RCPP
}
// compx
SEXP compx(arma::mat G, arma::mat X, arma::vec working, double sigma, int fam, int times);
RcppExport SEXP _eSCAN_compx(SEXP GSEXP, SEXP XSEXP, SEXP workingSEXP, SEXP sigmaSEXP, SEXP famSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type working(workingSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type fam(famSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(compx(G, X, working, sigma, fam, times));
    return rcpp_result_gen;
END_RCPP
}
// compCovw
SEXP compCovw(int p, int w_num, arma::mat Cov, arma::mat weights);
RcppExport SEXP _eSCAN_compCovw(SEXP pSEXP, SEXP w_numSEXP, SEXP CovSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type w_num(w_numSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cov(CovSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(compCovw(p, w_num, Cov, weights));
    return rcpp_result_gen;
END_RCPP
}
// matrix_flip
SEXP matrix_flip(arma::mat G);
RcppExport SEXP _eSCAN_matrix_flip(SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_flip(G));
    return rcpp_result_gen;
END_RCPP
}
// CCT_pval
SEXP CCT_pval(arma::vec x, arma::vec weights);
RcppExport SEXP _eSCAN_CCT_pval(SEXP xSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(CCT_pval(x, weights));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_eSCAN_mfwht", (DL_FUNC) &_eSCAN_mfwht, 3},
    {"_eSCAN_big_mfwht", (DL_FUNC) &_eSCAN_big_mfwht, 1},
    {"_eSCAN_compPval", (DL_FUNC) &_eSCAN_compPval, 2},
    {"_eSCAN_compCov", (DL_FUNC) &_eSCAN_compCov, 5},
    {"_eSCAN_compx", (DL_FUNC) &_eSCAN_compx, 6},
    {"_eSCAN_compCovw", (DL_FUNC) &_eSCAN_compCovw, 4},
    {"_eSCAN_matrix_flip", (DL_FUNC) &_eSCAN_matrix_flip, 1},
    {"_eSCAN_CCT_pval", (DL_FUNC) &_eSCAN_CCT_pval, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_eSCAN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}