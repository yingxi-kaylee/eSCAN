#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _eSCAN_big_mfwht(SEXP);
extern SEXP _eSCAN_CCT_pval(SEXP, SEXP);
extern SEXP _eSCAN_compCov(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _eSCAN_compCovw(SEXP, SEXP, SEXP, SEXP);
extern SEXP _eSCAN_compPval(SEXP, SEXP);
extern SEXP _eSCAN_compx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _eSCAN_matrix_flip(SEXP);
extern SEXP _eSCAN_mfwht(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_eSCAN_big_mfwht",   (DL_FUNC) &_eSCAN_big_mfwht,   1},
  {"_eSCAN_CCT_pval",    (DL_FUNC) &_eSCAN_CCT_pval,    2},
  {"_eSCAN_compCov",     (DL_FUNC) &_eSCAN_compCov,     5},
  {"_eSCAN_compCovw",    (DL_FUNC) &_eSCAN_compCovw,    4},
  {"_eSCAN_compPval",    (DL_FUNC) &_eSCAN_compPval,    2},
  {"_eSCAN_compx",       (DL_FUNC) &_eSCAN_compx,       6},
  {"_eSCAN_matrix_flip", (DL_FUNC) &_eSCAN_matrix_flip, 1},
  {"_eSCAN_mfwht",       (DL_FUNC) &_eSCAN_mfwht,       3},
  {NULL, NULL, 0}
};

void R_init_eSCAN(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}