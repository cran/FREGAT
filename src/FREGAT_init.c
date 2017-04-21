#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void Kernel_2wayIX(void *, void *, void *, void *);
extern void Kernel_IBS(void *, void *, void *, void *);
extern void Kernel_IBS_Weight(void *, void *, void *, void *, void *, void *);
extern void qfc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP asnum(SEXP);
extern SEXP readbed(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP subset(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"Kernel_2wayIX",     (DL_FUNC) &Kernel_2wayIX,      4},
    {"Kernel_IBS",        (DL_FUNC) &Kernel_IBS,         4},
    {"Kernel_IBS_Weight", (DL_FUNC) &Kernel_IBS_Weight,  6},
    {"qfc",               (DL_FUNC) &qfc,               11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"asnum",   (DL_FUNC) &asnum,   1},
    {"readbed", (DL_FUNC) &readbed, 5},
    {"subset",  (DL_FUNC) &subset,  3},
    {NULL, NULL, 0}
};

void R_init_FREGAT(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
