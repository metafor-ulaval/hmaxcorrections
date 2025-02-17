// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// H.cpp
double Hcpp(int n, cpp11::doubles P, cpp11::doubles h);
extern "C" SEXP _hmaxcorrections_Hcpp(SEXP n, SEXP P, SEXP h) {
  BEGIN_CPP11
    return cpp11::as_sexp(Hcpp(cpp11::as_cpp<cpp11::decay_t<int>>(n), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(P), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(h)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_hmaxcorrections_Hcpp", (DL_FUNC) &_hmaxcorrections_Hcpp, 3},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_hmaxcorrections(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
