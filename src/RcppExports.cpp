// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// EM
Rcpp::List EM(const Rcpp::List oldpara, const Rcpp::List data, const arma::vec lambda);
RcppExport SEXP _pairedfda_EM(SEXP oldparaSEXP, SEXP dataSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type oldpara(oldparaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(EM(oldpara, data, lambda));
    return rcpp_result_gen;
END_RCPP
}
// minEM
const List minEM(const List data, const arma::vec lambda, const int ka, const int kb, const double tol, int maxiter);
RcppExport SEXP _pairedfda_minEM(SEXP dataSEXP, SEXP lambdaSEXP, SEXP kaSEXP, SEXP kbSEXP, SEXP tolSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int >::type ka(kaSEXP);
    Rcpp::traits::input_parameter< const int >::type kb(kbSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(minEM(data, lambda, ka, kb, tol, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// loglike
double loglike(const List data, const List para);
RcppExport SEXP _pairedfda_loglike(SEXP dataSEXP, SEXP paraSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const List >::type para(paraSEXP);
    rcpp_result_gen = Rcpp::wrap(loglike(data, para));
    return rcpp_result_gen;
END_RCPP
}
// orth_algo
const List orth_algo(arma::mat Th, arma::mat V);
RcppExport SEXP _pairedfda_orth_algo(SEXP ThSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Th(ThSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(orth_algo(Th, V));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pairedfda_EM", (DL_FUNC) &_pairedfda_EM, 3},
    {"_pairedfda_minEM", (DL_FUNC) &_pairedfda_minEM, 6},
    {"_pairedfda_loglike", (DL_FUNC) &_pairedfda_loglike, 2},
    {"_pairedfda_orth_algo", (DL_FUNC) &_pairedfda_orth_algo, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_pairedfda(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
