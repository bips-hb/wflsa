// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// genlassoRcpp
Rcpp::DoubleVector genlassoRcpp(Rcpp::DoubleVector y, const Rcpp::NumericMatrix W, const size_t m, const size_t c, const double eta1, const double eta2, double a, const double rho, const int max_iter, const double eps, const double truncate);
RcppExport SEXP _wfla_genlassoRcpp(SEXP ySEXP, SEXP WSEXP, SEXP mSEXP, SEXP cSEXP, SEXP eta1SEXP, SEXP eta2SEXP, SEXP aSEXP, SEXP rhoSEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP truncateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DoubleVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type W(WSEXP);
    Rcpp::traits::input_parameter< const size_t >::type m(mSEXP);
    Rcpp::traits::input_parameter< const size_t >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type eta1(eta1SEXP);
    Rcpp::traits::input_parameter< const double >::type eta2(eta2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const double >::type truncate(truncateSEXP);
    rcpp_result_gen = Rcpp::wrap(genlassoRcpp(y, W, m, c, eta1, eta2, a, rho, max_iter, eps, truncate));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_wfla_genlassoRcpp", (DL_FUNC) &_wfla_genlassoRcpp, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_wfla(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
