#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void F77_NAME(fitnoisef)(double *qlev, int *nquant, double *esim, int *numboot, int *nobs, int *nvars, double *x, double *y, double *pf, int *dfmax, int *pmax, int *nlam, double *flmin, double *ulam, double *eps, int *isd, int *intr, int *maxit, int *nalam, double *b0, double *beta, int *ibeta, int *nbeta, double *alam, int *npass, int *jerr, double *quant, double *hatlam);

static const R_FortranMethodDef FortranEntries[] = {
  {"fitnoisef",  (DL_FUNC) &F77_NAME(fitnoisef), 28},
  {NULL, NULL, 0}
};

void R_init_FAS(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

